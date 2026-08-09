# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Nextflow DSL2 pipeline (single `main.nf`) that takes paired-end Illumina reads from single-amplified genomes (SAGs) through to decontaminated assemblies: QC → trimming → complexity filtering → k-mer normalization → host/contaminant read removal (BWA) → assembly (SPAdes) → contig trimming/dedup → host/contaminant contig removal (BLAST) → per-sample assembly stats.

## Running the pipeline

```bash
nextflow run main.nf                       # normal run, reads from ./input/, writes to ./results/
nextflow run main.nf --dev                 # dev mode: only processes the first sample (CH_fastq.take(1))
nextflow run main.nf --indir <dir> --output <dir>
```

There is no test suite, linter, or build step — the pipeline is validated by running it against `input/`.

`nextflow.config` is gitignored (see `.gitignore`) — it does not exist in a fresh clone and must be created locally before running. It currently sets:
```groovy
process.container = 'quay.io/nextflow/bash'
docker { enabled = true; autoMounts = true }
```
Docker must be running locally; all processes execute inside per-process containers. `docker.autoMounts` triggers an "Unrecognized config option" warning on newer Nextflow versions (harmless, but it's a sign this config predates the Nextflow version currently in use — see the conda env below).

Nextflow itself lives in the `nextflow-env` conda environment (`/opt/anaconda3/envs/nextflow-env`), not on the base PATH — activate that env before running `nextflow` directly.

### Input naming convention

Files in `params.indir` must be paired FASTQ named `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz` (`.fastq`/`.fq`/`.fq.gz` also matched). The workflow derives the sample ID by stripping `_R1`/`_R2` from the filename and pairs mates via `groupTuple(size:2)` — mismatched or unpaired files will break grouping.

## Pipeline architecture

One process per stage in `main.nf`, wired together in the `workflow {}` block. Tool versions are pinned via container tags (some also encoded in the process name, e.g. `_v0_11_9`, though a few processes have drifted from their name — check the actual `container` line, not just the process name).

**Read processing:** `FASTQC_v0_11_9` → `TRIMMOMATIC_v0_32` → `COMPLEXITY_FILTER` (drops low-complexity pairs, `templates/complexity_filter.py`) → `KMERNORM_v1_0_0` (k-mer normalization) → `DEINTERLEAVE` (bbmap `reformat.sh`, splits the interleaved normalized file back into separate `r1_norm_<ID>`/`r2_norm_<ID>` files).

**Read decontamination:** `CONTAM_READ_FINDER` (`bwa aln`/`bwa sampe` against `params.contam_ref_fasta`) → `CONTAM_READ_REPORTER` (`samtools view -F0x0004`, keeps only reads that actually aligned to the contaminant reference) → `CONTAM_READ_REMOVER` (`templates/contam_read_remover.py`, drops any read pair whose name appears in the reporter's hit list from the *pre-alignment* interleaved `normalized_pe_<ID>.fastq.gz`).

**Assembly:** `SPADES_v3_15_2` (`--sc --careful`, single-cell mode) → `TRIM_CONTIGS` (`templates/trim_and_deduplicate_contigs.py`: length-filters by `params.assembly_minlength`, trims `params.assembly_lefttrim`/`assembly_righttrim` bp off each end, and drops contigs that are exact reverse-complement duplicates of another contig — a known SPAdes v3 artifact that gets contigs rejected by NCBI/GenBank).

**Contig decontamination:** `CONTAM_CONTIG_FINDER` (`blastn` of trimmed contigs against `params.contam_ref_fasta`, using `params.contam_min_length`/`contam_min_percid` as the relevant thresholds) → `CONTAM_CONTIG_REMOVER` (`templates/contam_contig_remover.py`).

**Stats:** `MEASURE_SAG` — per-sample max contig length, total assembly length, GC%.

### BWA_INDEX / BLAST_INDEX: build-once reference indices

`params.contam_ref_fasta` points at a large (multi-GB) external reference (human+mouse, currently `/Users/greggavelis/Desktop/SCGC_Refcontam/GRCh38_AG665_mm10.fa` — an absolute host path outside the repo; update it per-machine). Both `BWA_INDEX` and `BLAST_INDEX` are conditional: the workflow block checks for existing index files *next to* `contam_ref_fasta` before deciding whether to run them, and both `publishDir` their output back to that same directory (`mode: 'copy'`) so the index persists across runs and clean `work/` dirs instead of rebuilding every time.

- BWA: checks for `<contam_ref_fasta>.{amb,ann,bwt,pac,sa}`.
- BLAST: more involved, because `makeblastdb` splits large references into numbered volumes and only writes a top-level `.nal` alias file — the workflow parses `.nal`'s `DBLIST` line and verifies each listed volume's `.nhr` actually exists, rather than trusting the alias file's mere presence (this repo's actual reference had a stale `.nal` referencing volumes that didn't exist on disk).

Both index channels (`CH_bwa_fasta`/`CH_bwa_index`, `CH_blast_fasta`/`CH_blast_index`) are threaded as real process inputs into `CONTAM_READ_FINDER` and `CONTAM_CONTIG_FINDER` respectively — this is what forces Nextflow to wait for a rebuild to finish before the step that needs it runs, not just process ordering in the script.

`CONTAM_READ_FINDER` also has `memory '9.GB'` and `maxForks 1` — Docker Desktop's VM has a real memory ceiling (was raised to ~11.9GB from the default 7.75GB to make this reference usable at all) and `bwa aln` against this reference is the only step that comes close to it; `maxForks 1` stops two samples from stacking their memory demand past that ceiling even though each fits individually.

### Stepwise count logging

Per-sample `.count` files (one CSV row each) are emitted by most stages and rolled up with `channel.collectFile(seed: COUNT_HEADER, storeDir: CH_stepwise_counts)` directly in the workflow block — no separate `LOG_*` processes. Output is `${params.output}/sample_tracking/stepwise_counts/{1..8}_<stage>counts.csv`. `sort:` is intentionally omitted, so row order within a file reflects task completion order, not sample order. Step 8 (`8_clean_contigcounts.csv`, i.e. post-`CONTAM_CONTIG_REMOVER`) has no `collectFile()` call yet — nothing downstream of `TRIM_CONTIGS` currently emits a per-sample count in the right shape for it.

## Known gaps / quirks

- `CONTAM_CONTIG_REMOVER` has no active `publishDir` — its `publishDir` lines are all commented out (leftover from a shared/multi-project template, referencing an undefined `DIR_out` and `params.SPC`). Its outputs (`SCGC_<ID>_contigs.fasta`, the final decontaminated assembly) currently only exist in Nextflow's `work/` dir unless you add one.
- Several unrelated steps share the generic `brwnj/kmernorm:v1.0.0` image as a general-purpose Python/bash container (`COMPLEXITY_FILTER`, `KMERNORM_v1_0_0`, `CONTAM_READ_REMOVER`, `TRIM_CONTIGS`, `MEASURE_SAG`, `CONTAM_CONTIG_REMOVER`) — it's not just for k-mer normalization despite the name.
- `params.reference_threshold` (bwa `aln -n`) and `params.contam_min_length`/`contam_min_percid` (blast filtering thresholds) are declared but not all necessarily enforced in the same process that declares related params — check each `template_*.py`/shell block directly rather than assuming from the param's declared location.

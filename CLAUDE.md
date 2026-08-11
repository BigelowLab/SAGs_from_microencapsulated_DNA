# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository. `README.md` covers the human-facing quick start/usage/parameter reference — this file focuses on things a human README wouldn't (internal architecture, non-obvious bugs and their fixes, why things are built the way they are).

## What this is

A Nextflow DSL2 pipeline (single `main.nf`) that takes paired-end Illumina reads from single-amplified genomes (SAGs) through to decontaminated assemblies: QC → trimming → complexity filtering → k-mer normalization → host/contaminant read removal (BWA) → assembly (SPAdes) → contig trimming/dedup → host/contaminant contig removal (BLAST) → CheckM completeness → a single per-sample stats table.

## Running the pipeline

```bash
nextflow run main.nf                       # normal run, reads from ./input/, writes to ./results/
nextflow run main.nf --dev                 # dev mode: only processes the first sample (CH_fastq.take(1))
nextflow run main.nf --indir <dir> --output <dir>
nextflow run main.nf -stub-run --dev       # fast DAG/wiring smoke test — see "Stub-run support" below
```

There is no test suite or linter — `-stub-run` (see below) is the closest thing to one.

`nextflow.config` is gitignored (see `.gitignore`) — it does not exist in a fresh clone and must be created locally before running. It currently sets:
```groovy
process.container = 'quay.io/nextflow/bash'
docker { enabled = true; autoMounts = true }
```
Docker must be running locally; all processes execute inside per-process containers.

Nextflow itself lives in the `nextflow-env` conda environment (`/opt/anaconda3/envs/nextflow-env`), not on the base PATH — activate that env before running `nextflow` directly.

### Input naming convention

Files in `params.indir` must be paired FASTQ named `<sample>_R1.fastq.gz` / `<sample>_R2.fastq.gz` (`.fastq`/`.fq`/`.fq.gz` also matched). The workflow derives the sample ID by stripping `_R1`/`_R2` from the filename and pairs mates via `groupTuple(size:2)` — mismatched or unpaired files will break grouping.

## Pipeline architecture

One process per stage in `main.nf`, wired together in the `workflow {}` block. Tool versions are pinned via container tags (some also encoded in the process name, e.g. `_v0_11_9`, though a few processes have drifted from their name — check the actual `container` line, not just the process name).

**Read processing:** `FASTQC_v0_11_9` → `TRIMMOMATIC_v0_32` → `COMPLEXITY_FILTER` (drops low-complexity pairs, `templates/complexity_filter.py`) → `KMERNORM_v1_0_0` (k-mer normalization) → `DEINTERLEAVE` (bbmap `reformat.sh`, splits the interleaved normalized file back into separate `r1_norm_<ID>`/`r2_norm_<ID>` files). Per-sample readcounts (`Raw_readcount`, `Trimmed_readcount`, etc.) sum both mates (r1 + r2), not just r1.

**Read decontamination:** `CONTAM_READ_FINDER` (`bwa aln`/`bwa sampe` against `params.contam_ref_fasta`) → `CONTAM_READ_REPORTER` (`samtools view -F0x0004`, keeps only reads that actually aligned to the contaminant reference) → `CONTAM_READ_REMOVER` (`templates/contam_read_remover.py`, drops any read pair whose name appears in the reporter's hit list from the *pre-alignment* interleaved `normalized_pe_<ID>.fastq.gz`).

**Assembly:** `SPADES_v3_15_2` (`--sc --careful`, single-cell mode) → `TRIM_CONTIGS` (`templates/trim_and_deduplicate_contigs.py`: length-filters by `params.assembly_minlength`, trims `params.assembly_lefttrim`/`assembly_righttrim` bp off each end, and drops contigs that are exact reverse-complement duplicates of another contig — a known SPAdes v3 artifact that gets contigs rejected by NCBI/GenBank).

**Contig decontamination:** `CONTAM_CONTIG_FINDER` (`blastn` of trimmed contigs against `params.contam_ref_fasta`, using `params.contam_min_length`/`contam_min_percid` as the relevant thresholds) → `CONTAM_CONTIG_REMOVER` (`templates/contam_contig_remover.py`, publishes the final `SCGC_<ID>_contigs.fasta`).

**Stats:** `MEASURE_SAG` (max contig length, total assembly length, GC%) and `CHECKM_v1_1_9` (`checkm lineage_wf --reduced_tree`, genome completeness) run in parallel off `TRIM_CONTIGS.out.trimmed_contigs`, then `ASSEMBLY_STATS_TABULATOR` pivots every stage's counts into one wide `assembly_stats.csv` (one row per sample) — see "Stepwise count logging" below.

### BWA_INDEX / BLAST_INDEX: self-bootstrapping reference

`params.contam_ref_fasta` defaults to `./reference/GRCh38_AG665_mm10.fa` — relative and gitignored, so a fresh clone has nothing there. Both `BWA_INDEX` and `BLAST_INDEX` are conditional: the workflow block checks for existing files next to `contam_ref_fasta` before deciding whether to run them, and `log.info` prints a one-time-cost notice right before either actually triggers.

- **`BWA_INDEX`** no longer runs `bwa index` locally — it downloads the fasta *and* the prebuilt `.amb/.ann/.bwt/.pac/.sa` from Zenodo (DOI `10.5281/zenodo.21682938`, "GORG Dark - Reference Contaminant Dataset") and verifies each against Zenodo's published MD5 before extracting. Several non-obvious things about this process, all verified empirically (not documented anywhere) and all found by actually failing on real Linux CI, not just local macOS testing:
  - All 6 files served by that record are **zip archives**, each wrapping the real, already-decompressed file under its real name (fetching `<fasta>.ann.zip` returns a zip whose sole entry is literally `<fasta>.ann`) — extract with `unzip`, not `gunzip`.
  - The container is `curlimages/curl:8.21.0`, **not** a general-purpose image — it bakes in `ENTRYPOINT ["curl"]`, and BusyBox's `wget` (present in other biocontainers used elsewhere in this pipeline) cannot complete an HTTPS connection to zenodo.org at all. Because of the fixed entrypoint, this process needs `--entrypoint ""` (or Nextflow's normal `/bin/bash -ue .command.sh` invocation actually runs as `curl /bin/bash -ue .command.sh`, which curl silently misparses as `--user e` plus two bogus URLs — a confusing failure mode with no obvious link to the real cause).
  - This image also runs as a **non-root user** by default. Docker Desktop's bind-mount layer on macOS is lenient about that mismatch against the host-owned work dir, but real Linux Docker (GitHub Actions runners, HPC nodes) enforces it for real, so writing the downloaded/extracted files fails with `Permission denied` — fixed with `-u root` alongside the entrypoint override: `containerOptions '--entrypoint "" -u root'`.
  - This Alpine-based image has no `/bin/bash` at all, only `/bin/sh` (BusyBox ash) — `shell '/bin/sh', '-ue'` is required. Note the `shell` directive's syntax in this Nextflow version: bare comma-separated strings, not a `[...]` list literal — `shell ['/bin/sh', '-ue']` fails to parse ("Invalid process directive").
  - Because the shell is `/bin/sh`, not bash, **bash-only syntax silently breaks** rather than erroring — the `stub:` block originally used brace expansion (`touch ${fasta}{,.amb,.ann,...}`), which `sh` doesn't support. Under `sh` that whole expression is one literal filename, not six, so `touch` exits `0` while never producing the actually-expected output (`Missing output file(s)` from Nextflow, not a shell error). Fixed by listing each filename explicitly instead of relying on brace expansion.
- **`BLAST_INDEX`** still builds locally via `makeblastdb`, now depends on `CH_bwa_fasta` (not a raw `file()` variable) so it properly waits on `BWA_INDEX`'s download when the fasta doesn't exist yet, rather than relying on script ordering alone.
- BWA readiness check: `<contam_ref_fasta>` itself exists *and* `.{amb,ann,bwt,pac,sa}` all exist (the fasta is no longer assumed to pre-exist, so this can't use `checkIfExists: true`).
- BLAST readiness check: more involved, because `makeblastdb` splits large references into numbered volumes and only writes a top-level `.nal` alias file — the workflow parses `.nal`'s `DBLIST` line and resolves each listed volume against `contam_ref_fasta`'s own directory (the entries are bare filenames), verifying the volume's `.nhr` actually exists rather than trusting the alias file's mere presence.
- Both index channels (`CH_bwa_fasta`/`CH_bwa_index`, `CH_blast_fasta`/`CH_blast_index`) are threaded as real process inputs into `CONTAM_READ_FINDER` and `CONTAM_CONTIG_FINDER` — this is what forces Nextflow to wait for a download/rebuild to finish before the step that needs it runs, not just process ordering in the script.
- Both processes have `cache false` and `publishDir ..., enabled: !workflow.stubRun`. The `enabled` guard is load-bearing, not defensive boilerplate — without it, `-stub-run`'s placeholder output gets copied into the real, shared reference directory (this happened once for real: 0-byte `.nhr`/`.nin`/`.nsq` files leaked into the reference dir from a stub run before the guard existed). `cache false` closes a related risk — it's not confirmed whether Nextflow's task hash distinguishes a `stub:` execution from a `script:` execution of the same process, so without it a stub run's cached result could in principle be reused by a later real `-resume` run, or vice versa.
- `CONTAM_READ_FINDER` and `CHECKM_v1_1_9` are the two steps that come close to Docker Desktop's memory ceiling locally. Rather than a hardcoded `memory` value tuned to one machine (which broke portability to CI's smaller runners), both use `memory { N.GB * task.attempt }` with `errorStrategy { task.exitStatus in [137, 140] ? 'retry' : 'terminate' }` and `maxRetries 3` — retries with more memory on an OOM kill instead of assuming a fixed ceiling, adapting to whatever's actually available (laptop, CI, HPC node). `CONTAM_READ_FINDER` starts at `5.GB` (its real historical OOM point against this reference); `CHECKM_v1_1_9` starts at `6.GB` (its `--reduced_tree` flag only guarantees "<16GB", covered across two attempts). `maxForks 1` is kept on both regardless of environment, so retries on one sample can't stack against another sample's concurrent attempt.

### Stub-run support

10 of the pipeline's ~19 processes (`FASTQC_v0_11_9`, `TRIMMOMATIC_v0_32`, `COMPLEXITY_FILTER`, `KMERNORM_v1_0_0`, `BWA_INDEX`, `BLAST_INDEX`, `CONTAM_READ_FINDER`, `SPADES_v3_15_2`, `CONTAM_CONTIG_FINDER`, `CHECKM_v1_1_9`) have a `stub:` block, covering every step that's either genuinely slow (SPAdes, CheckM) or depends on the multi-GB reference (bwa/blast alignment and indexing). The rest are fast enough to just run for real even under `-stub-run` — Nextflow falls back to a process's normal `script:`/`shell:` when no `stub:` is defined, so partial coverage is safe.

Two things make the stubs more than placeholder `touch` calls:
- `TRIMMOMATIC_v0_32` → `COMPLEXITY_FILTER` → `KMERNORM_v1_0_0`'s stubs chain a **real subset** of actual read data (first 10 pairs) through each other, rather than empty files, because `DEINTERLEAVE` right after `KMERNORM_v1_0_0` is *not* stubbed and genuinely parses whatever it's given with `reformat.sh`.
- `CONTAM_CONTIG_FINDER`'s stub writes the same header line a real zero-hit run would produce, not an empty file, since `CONTAM_CONTIG_REMOVER` downstream isn't stubbed either and needs a properly-shaped TSV.

`ASSEMBLY_STATS_TABULATOR` prefixes every cell in its output with `STUB-RUN ` when `workflow.stubRun` is true, so stub output can never be mistaken for a real assembly's stats.

### Continuous integration

`.github/workflows/stub-run.yml` runs `nextflow run main.nf -stub-run --dev --indir .github/test_data -with-docker` on every push/PR, against tiny committed synthetic reads (`.github/test_data/CI_TEST_R{1,2}.fastq.gz` — `input/` has real sample data but is intentionally untracked, so CI needs its own). Two things this workflow has to work around that a local run doesn't:
- `nextflow.config` is gitignored, so there's nothing in a fresh CI checkout to turn Docker on — `-with-docker` does that directly on the CLI instead.
- The real reference doesn't exist in CI either, so `BWA_INDEX`/`BLAST_INDEX` run in stub mode too (placeholder files, no real ~8GB Zenodo download) — this is exactly the case `-stub-run` was built for.

Every bug fixed in the "self-bootstrapping reference" section above (`-u root`, the brace-expansion stub bug) was only discovered because this workflow failed on real Linux — none of it surfaced in local macOS testing. If `BWA_INDEX`/`BLAST_INDEX` change again, treat a green run here as the actual bar, not local success alone.

### Stepwise count logging

Per-sample `.count`/TSV outputs from each stage are rolled up with `channel.collectFile(seed: COUNT_HEADER, storeDir: CH_stepwise_counts, cache: false)` directly in the workflow block — no separate `LOG_*` processes. `cache: false` is required on all of these: `collectFile`'s own cache is tied to `-resume` and keyed on input data, so a pure change to the aggregation logic with unchanged upstream data can silently reuse a stale result otherwise (this happened for real once). Output is `${params.output}/sample_tracking/stepwise_counts/{1..9}_<stage>.csv`, one file per stage including `9_checkm_completeness.csv` (built via a `.map()` that pulls the `Completeness` column out of CheckM's TSV directly in Groovy — no separate parsing process/container needed).

All per-stage channels then `.concat()` (not `.mix()`) into one `all_stepwise_counts.csv`, each file's own header line stripped before stacking so the combined file has exactly one. `.concat()` guarantees stage order regardless of which stage's tasks actually finish first; `.mix()` would interleave by completion order instead. `collectFile`'s default `sort` (`'hash'`) reorders entries by content hash regardless of input order, which was silently undoing the `.concat()` ordering — `sort: false` ("append as produced") is required alongside it.

Finally, `ASSEMBLY_STATS_TABULATOR` reads `all_stepwise_counts.csv` and pivots it (long → wide) into `${params.output}/assembly_stats.csv`, one row per sample, columns in a fixed order (`LIST_col_order`). Handles the `Contam_filtered_readcount == "NO_CHANGE"` case (when `CONTAM_READ_REMOVER` found nothing to remove) by copying `Normalized_readcount` in as the numeric value.

## Known gaps / quirks

- Several unrelated steps share the generic `brwnj/kmernorm:v1.0.0` image as a general-purpose Python/bash container (`COMPLEXITY_FILTER`, `KMERNORM_v1_0_0`, `CONTAM_READ_REMOVER`, `TRIM_CONTIGS`, `MEASURE_SAG`, `CONTAM_CONTIG_REMOVER`, `ASSEMBLY_STATS_TABULATOR`) — it's not just for k-mer normalization despite the name. It's also a personal Docker Hub image with no `Dockerfile` in this repo, unlike everything else which pulls from quay.io/biocontainers.
- `params.reference_threshold` (bwa `aln -n`) and `params.contam_min_length`/`contam_min_percid` (blast filtering thresholds) are declared but not all necessarily enforced in the same process that declares related params — check each `template_*.py`/shell block directly rather than assuming from the param's declared location.
- `docker.autoMounts` in `nextflow.config` triggers an "Unrecognized config option" warning on the Nextflow version currently in use (harmless).

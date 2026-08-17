# SAGs from microencapsulated DNA

![Pipeline smoke test](https://github.com/BigelowLab/GORG-Dark-SAG-assembly/actions/workflows/stub-run.yml/badge.svg)

A Nextflow pipeline that takes paired-end Illumina reads from Atrandi combinatorial-barcoded single-amplified genome (SAG) libraries through to decontaminated assemblies: Atrandi demultiplexing → quality control → trimming → complexity filtering → k-mer normalization → host/contaminant read removal → assembly → contig trimming/deduplication → host/contaminant contig removal → genome completeness estimation → a single per-sample stats table.

Built for the GORG-Dark project (single-cell genomics of deep-ocean prokaryotes), but the decontamination steps are general-purpose against any BWA/BLAST-indexable reference.

## Requirements

- [Nextflow](https://www.nextflow.io/) (DSL2)
- [Docker](https://www.docker.com/), running locally — every process except `KMERNORM_v1_0_0` executes in its own container, no local tool installation needed
- Java (required by Nextflow itself)
- `kmernorm` on `PATH` — see [Installing kmernorm](#installing-kmernorm) below. Not bundled in a container; you provide your own build.

## Quick start

```bash
git clone https://github.com/BigelowLab/GORG-Dark-SAG-assembly.git
cd GORG-Dark-SAG-assembly
```

Create a local `nextflow.config` (gitignored — machine-specific, not part of the repo):

```groovy
process.container = 'quay.io/nextflow/bash'
docker.enabled = true
```

Drop paired, Atrandi-barcoded FASTQ files into `input/` (see naming convention below), then:

```bash
nextflow run main.nf
```

The first run downloads the decontamination reference (~8GB) automatically — see [Reference data](#reference-data) below.

### Try it without any data first

```bash
nextflow run main.nf -stub-run --dev --indir .github/test_data --output /tmp/gorg-test -with-docker
```

Runs the full pipeline wiring end-to-end against tiny bundled synthetic Atrandi-barcoded reads in under a minute. Atrandi demultiplexing runs for real (it's fast even on real data), while every slow step further downstream (assembly, alignment, indexing, CheckM) is swapped for a fast placeholder. This is also what runs in CI on every push/PR (`.github/workflows/stub-run.yml`) — good for confirming your Nextflow/Docker setup works, or for testing a pipeline change without waiting on a real run.

## Usage

```bash
nextflow run main.nf                               # normal run, reads from ./input/, writes to ./results/
nextflow run main.nf -resume                        # resume after a failure, reusing cached results
nextflow run main.nf --dev                          # only carry the first few capsules into assembly
nextflow run main.nf --dev --dev_num_capsules 5      # same, but carry 5 capsules instead of the default 3
nextflow run main.nf --indir <dir> --output <dir>   # override input/output locations
```

### Input naming convention

Files in `--indir` must be paired FASTQ named `<library>_R1.fastq.gz` / `<library>_R2.fastq.gz` (`.fastq`, `.fq`, and `.fq.gz` are also matched). The library ID is derived by stripping `_R1`/`_R2` from the filename; mismatched or unpaired files will break the pairing step.

Each pair is an **Atrandi combinatorial-barcode pool**, not a single SAG — many single-cell capsules multiplexed together via 4 independent 8bp barcode positions (D, C, B, A) embedded in the first 44bp of R2. The pipeline demultiplexes each pool before assembly; downstream stages then operate per-capsule, with each capsule's sample ID taking the form `<library>_<capsuleID>`.

### Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `--indir` | `./input/` | Directory of input FASTQ files |
| `--output` | `./results/` | Output directory |
| `--dev` | `false` | Only carry `--dev_num_capsules` capsules into assembly (every pool is still demultiplexed in full — see [Pipeline stages](#pipeline-stages)) |
| `--dev_num_capsules` | `3` | How many capsules `--dev` carries into assembly |
| `--publishmode` | `symlink` | How outputs are linked into `--output` (`symlink`, `copy`, etc.) |
| `--barcode_dir` | `./barcodes/` | Directory containing `bc{A,B,C,D}_24.txt`, the 4 Atrandi barcode lists |
| `--read_threshold` | `3` | Minimum reads a barcode combo needs to be treated as a real capsule, not noise |
| `--cell_threshold` | (uncapped) | Maximum number of capsules kept per pool, after ranking by read count |
| `--barcode_trim_length` | `45` | bp trimmed off the front of R2 (barcode + linker) after demultiplexing |
| `--sample_num_reads` | `1000000` | Reads subsampled per pool to estimate barcode frequencies before the real demux |
| `--sample_hamming_dist` / `--split_hamming_dist` | `1` / `1` | Barcode mismatch tolerance during frequency estimation / the real demux |
| `--contam_ref_fasta` | `./reference/GRCh38_AG665_mm10.fa` | Decontamination reference; auto-downloaded here if missing (see below) |
| `--complexity_threshold` | `0.05` | Low-complexity read filtering threshold |
| `--reference_threshold` | `0.05` | `bwa aln -n` mismatch threshold for read decontamination |
| `--contam_min_length` / `--contam_min_percid` | `100` / `95.0` | BLAST thresholds for contig decontamination |
| `--assembly_minlength` | `1000` | Minimum contig length kept after assembly |
| `--assembly_lefttrim` / `--assembly_righttrim` | `200` / `200` | bp trimmed off each contig end |
| `--kmernorm_opts` | `-k 21 -t 30 -c 3` | Passed directly to `kmernorm` |

## Pipeline stages

0. **Atrandi demultiplexing** — subsamples each pool to estimate barcode frequencies (Pheniqs), filters out low-count noise, assigns each real capsule an ID, then demultiplexes the full pool by combinatorial D/C/B/A barcode into one read pair per capsule, and trims the barcode/linker bases off R2
1. **QC & trimming** — FastQC, Trimmomatic
2. **Complexity filtering & normalization** — drops low-complexity pairs, then k-mer normalizes
3. **Read decontamination** — aligns to the reference with BWA, removes anything that hits it
4. **Assembly** — SPAdes (`--sc --careful`, single-cell mode)
5. **Contig trimming & dedup** — length filter, end-trimming, drops exact reverse-complement duplicate contigs (a known SPAdes artifact that gets genomes rejected by NCBI/GenBank)
6. **Contig decontamination** — BLASTs trimmed contigs against the reference, removes hits
7. **Stats** — assembly length/GC%/max contig length, plus CheckM genome completeness

Each stage's per-sample counts land in `results/sample_tracking/stepwise_counts/`, and everything gets combined into one final `results/assembly_stats.csv` — one row per sample.

## Reference data

Read and contig decontamination both need a reference to screen against (default: a combined human GRCh38 + mouse mm10 assembly used to catch common lab contaminants). It's not part of this repo — the `BWA_INDEX` and `BLAST_INDEX` processes each check whether it already exists at `--contam_ref_fasta`'s path, and download it from Zenodo automatically the first time it's needed:

> **GORG Dark - Reference Contaminant Dataset**
> DOI: [10.5281/zenodo.21682938](https://doi.org/10.5281/zenodo.21682938)

That download (fasta + prebuilt BWA index) is a one-time cost of a few GB; a matching BLAST database is then built locally once and also cached in place. Both are skipped automatically on every run after the first. Point `--contam_ref_fasta` at an already-populated path (e.g. a shared location on a cluster) to avoid downloading a fresh copy per clone.

## Continuous integration

Every push and pull request runs a fast wiring smoke test via GitHub Actions (`.github/workflows/stub-run.yml`) — see [Try it without any data first](#try-it-without-any-data-first) above for the equivalent local command.

## Installing kmernorm

`KMERNORM_v1_0_0` is the one step that does not run in a container. The underlying `kmernorm` tool (by Mingkun Li) has no license anywhere — not on its SourceForge page, not in its source archive — so this pipeline does not bundle it or depend on any third-party Docker image wrapping it. You need to build and provide your own copy on `PATH` before running the pipeline for real (`-stub-run` does not need it — see above).

```bash
curl -sL -o kmernorm.tar https://sourceforge.net/projects/kmernorm/files/latest/download
tar xf kmernorm.tar && make
cp kmernorm /opt/homebrew/bin/   # or anywhere else already on PATH
```

Installing and running this source is entirely at your own discretion — confirm license terms and satisfy yourself of the tool's suitability independently; this pipeline doesn't warrant or vouch for it. See `THIRD_PARTY_LICENSES.md` for more.

> **Warning — do not build kmernorm on macOS.** A macOS/ARM64 build of this exact source silently corrupts paired-end read ordering (verified: R1/R2 pairs get mismatched partway through real data, with no error at build or run time — it only surfaces later as a cryptic `bwa` "paired reads have different names" failure). The identical source built for Linux/x86_64 produces correctly paired output on the same data — this is specific to compiling on Mac. Build and run `kmernorm` on Linux until this is root-caused.

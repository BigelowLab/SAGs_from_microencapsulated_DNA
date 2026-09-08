# Third-Party Licenses

This pipeline is licensed under the GNU General Public License v3.0 or later
(see `LICENSE`). It wraps and invokes the following third-party
bioinformatics tools as separate external programs. Each tool is distributed
under, and remains subject to, its own license. Using this pipeline does not
relicense these tools; the GPL-3.0-or-later license applies only to the
pipeline's own code (Nextflow scripts, process definitions, and wrapper
logic).

If container images built from this repository bundle any of these tools'
binaries or source, the original license and copyright notices for each tool
must be preserved and included in the distributed image/artifact. For example, the viral annotators **DeepVirFinder** and **GeNomad** are not licensed for commercial use. For that reason, the --viral mode of this pipeline is set to 'false' by default, and should only be activated by non-commercial users (e.g. academic researchers). This pipeline does not bundle or distribute these third party tools, which must be installed at the discretion of the user. Likewise, reference *databases* (the contaminant reference, SILVA, GTDB, the Prokka SwissProt set) are user-supplied and not redistributed
by this pipeline — their own terms still apply and are noted where relevant.

## Read processing & demultiplexing

| Tool | Version | License | Notes |
|---|---|---|---|
| [seqtk](https://github.com/lh3/seqtk) | 1.2 | MIT | `SAMPLE_READS` — subsamples R2 to estimate Atrandi barcode frequencies. |
| [Pheniqs](https://github.com/biosails/pheniqs) | 2.1.0 | **Conflicting** | Individual source files carry AGPL-3.0-or-later headers, but the repository's top-level `LICENSE` states a separate, more restrictive NYU research license (internal, non-commercial use only; no redistribution/modification/sublicensing without a signed agreement). The intent (see the process comments and `README.md` "Installing Pheniqs") is for the user to install and provide their own `pheniqs` binary on `PATH`; installing/running it is entirely **at the user's discretion and risk**, and this pipeline does not warrant or vouch for it. |
| [matplotlib](https://matplotlib.org/) | (biocontainers mulled image) | Matplotlib License (BSD-style, PSF-derived) | `PHENIQS_PLOT_HIST` — the observed-barcode-distribution histogram. The image also bundles NumPy (BSD-3-Clause) and DejaVu fonts (permissive Bitstream Vera / public-domain additions). |
| [samtools](https://github.com/samtools/samtools) / [htslib](https://github.com/samtools/htslib) | 1.24 | MIT/Expat (samtools); MIT + modified-BSD (htslib) | `PHENIQS_COUNT_SORT_SAMPLE` (tally observed barcodes) and `CONTAM_READ_REPORTER` (`samtools view -F0x0004`). |
| [Trim Galore](https://github.com/FelixKrueger/TrimGalore) | 0.6.7 | GPL-3.0-or-later | `TRIM_BARCODE` — trims the Atrandi barcode + linker bases off the front of R2 after demux. Wraps **Cutadapt** (MIT-licensed) and FastQC internally. |
| [Cutadapt](https://github.com/marcelm/cutadapt) | (bundled in Trim Galore 0.6.7) | MIT | Invoked indirectly by Trim Galore. |
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.9 | GPL-3.0-or-later | `FASTQC_v0_11_9`. Bundles Picard SAM/BAM libraries and other third-party code under separate terms. |
| [Trimmomatic](https://github.com/usadellab/Trimmomatic) | 0.32 | GPL-3.0-or-later | `TRIMMOMATIC_v0_32`. Illumina adapter-sequence FASTA files included with the distribution are **not** GPL-licensed — owned by and used with permission of Illumina, Inc. Handle separately if redistributing. |
| [pysam](https://github.com/pysam-developers/pysam) (+ htslib) | 0.24.0 | MIT | Only real dependency of `COMPLEXITY_FILTER` (`templates/complexity_filter.py`), `CONTAM_READ_REMOVER` (`templates/contam_read_remover.py`), and `SSU_GET_GENE` (`pysam.faidx` / `FastaFile` to excise the SSU region from its hit contig). |
| [KMERNORM](https://sourceforge.net/projects/kmernorm/) | 1.0.0 | **Unstated / unverified** | No license file on the project page or in its source archive as of this writing. The user must build and provide their own `kmernorm` binary on `PATH`; installing/running either the tool or that image is entirely **at the user's discretion and risk**, and this pipeline does not warrant or vouch for it. See `README.md` for a known correctness issue with macOS-built copies. |
| [BBMap / BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/) | 38.90 | BSD-3-Clause-LBNL | `DEINTERLEAVE` (`reformat.sh`). Modified-BSD variant used by Lawrence Berkeley National Laboratory / JGI. |

## Reference download & decontamination

| Tool | Version | License | Notes |
|---|---|---|---|
| [Debian](https://www.debian.org/legal/licenses/) (`debian:bookworm-slim`) | bookworm | Composite (base system: GPL / LGPL / MIT / BSD / others per package) | Base image for `BWA_INDEX`, which `apt-get install`s `curl` + `unzip` to fetch the contaminant reference from Zenodo. No Debian source is redistributed by this repo; the image is pulled at runtime. |
| [BWA](https://github.com/lh3/bwa) (Burrows-Wheeler Aligner) | 0.7.17 | GPL-3.0-or-later + MIT | `CONTAM_READ_FINDER` (`bwa aln` / `bwa sampe`). GPLv3 governs the package as a whole (via BWT-SW-derived code); the sorting, hash-table, BWT, and IS libraries are separately MIT-licensed. |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | 2.11.0 | Public Domain (U.S. Government work) | Written by NCBI staff in the course of U.S. federal employment; not copyrightable in the U.S.; no warranty from NLM / U.S. Government. `CONTAM_CONTIG_FINDER`, `SSU_BLAST` (`-task megablast` vs SILVA), and `BLAST_INDEX` (`makeblastdb`). |
| GORG Dark Reference Contaminant Dataset | Zenodo DOI [10.5281/zenodo.21682938](https://doi.org/10.5281/zenodo.21682938) | (see the Zenodo record) | The combined human GRCh38 + mouse mm10 (+ `AG665`) contaminant reference + prebuilt BWA index, auto-downloaded by `BWA_INDEX`. Not part of this repo; terms are those of the Zenodo deposit. |

## Assembly, contig processing & completeness

| Tool | Version | License | Notes |
|---|---|---|---|
| [SPAdes](https://github.com/ablab/spades) | 3.15.2 | GPL-2.0-only | `SPADES_v3_15_2` (`--sc --careful`). |
| [Biopython](https://biopython.org/) | 1.84 (1.79 in the `PROKKA_GFF_2_TSV` mulled image) | Biopython License Agreement (permissive; some files dual-licensed BSD-3-Clause) | `TRIM_CONTIGS` / `MEASURE_SAG` / `CONTAM_CONTIG_REMOVER` templates; `templates/prokka_gff_2_tsv.py` (`Bio.SeqIO`, coding-density); `templates/ssu_classifier.py` (`Bio.Phylo` / `Bio.SeqIO`, SSU LCA classification). |
| [toolshed](https://github.com/brentp/toolshed) | 0.4.8 | MIT | `pip`-installed at task start inside `CONTAM_CONTIG_REMOVER` (`templates/contam_contig_remover.py`). |
| [InterLap](https://github.com/brentp/interlap) | 0.2.7 | MIT | `pip`-installed at task start inside `CONTAM_CONTIG_REMOVER` — interval-overlap logic. |
| [CheckM](https://github.com/Ecogenomics/CheckM) | 1.1.9 | GPL-3.0-or-later | `CHECKM_v1_1_9` (`lineage_wf --reduced_tree`). CheckM v1 is unmaintained upstream (superseded by CheckM2); this does not affect its license terms. |

## Annotation & taxonomy

| Tool | Version | License | Notes |
|---|---|---|---|
| [Prokka](https://github.com/tseemann/prokka) | 1.14.6 | GPL-3.0-only | `PROKKA_v1_14_6`. At runtime invokes several separately-licensed components (BLAST+, HMMER, Aragorn, Barrnap, Infernal, MinCED, Prodigal, tbl2asn) as external binaries within its container — see Prokka's own `doc/LICENSE.*` files. `--proteins` points at a user-supplied SwissProt FASTA (UniProt; CC BY 4.0), not redistributed here. |
| [pandas](https://pandas.pydata.org/) | 2.2.1 (1.5.2 in the `PROKKA_GFF_2_TSV` mulled image) | BSD-3-Clause | `ASSEMBLY_STATS_TABULATOR`, the Pheniqs helper scripts, `PROKKA_GFF_2_TSV`, `PARSE_CLASSIFIER`, `PARSE_GTDBTK`, `EGGNOG_HITS_TO_CELL_OR_VIRUS`. |
| [NumPy](https://numpy.org/) | — | BSD-3-Clause | Bundled with pandas / matplotlib in the images above. |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | 2.0.0 | GPL-3.0-or-later | `GTDBTK_v2_0_0` (`classify_wf`). Pinned to 2.0.0 (with GTDB reference data **release 207**) for fidelity to the published results, not because it is current. GTDB reference data is released under CC BY-SA 4.0 and is user-supplied, not redistributed here. |
| [SILVA rRNA database](https://www.arb-silva.de/) / CREST | `silvamod` (v128-era) | SILVA: CC BY 4.0. CREST (Lanzén et al.): GPL-3.0 | User-supplied SILVA BLAST DB + `.map` / `.tree`, not redistributed by this pipeline. The lowest-common-ancestor logic in `templates/ssu_classifier.py` is derived from CREST (Lanzén A. *et al.*, 2012). |

## Viral classification (`--viral`, default off)

| Tool | Version | License | Notes |
|---|---|---|---|
| [geNomad](https://github.com/apcamargo/genomad) | 1.11.1 | **Berkeley Lab Academic / Non-Commercial License** | `GENOMAD_v1_11_1` (`end-to-end`). **Restricts *use* itself** — academic, internal research & development, non-commercial only; commercial rights reserved by Lawrence Berkeley National Laboratory. Confirm eligibility before any commercial use of this pipeline or its results. `params.DB_genomad_v1_11_1` is a user-supplied prebuilt database. |
| [VirSorter2](https://github.com/jiarong/VirSorter2) (`docker://jiarong/virsorter:2.2.3`) | 2.2.3 | GPL-2.0-or-later | `VIRSORTER_v2_2_3` calls `virsorter run` — the VirSorter2 CLI, not the original VirSorter (`simroux/VirSorter`, a different codebase). Pinned to the `2.2.3` tag (was `:latest`); confirmed identical content to `:latest` via `virsorter --version` (2026-09-04) — `:latest` hasn't moved since both were published together on 2021-12-27, but `:2.2.3` is immutable going forward. |
| [CheckV](https://bitbucket.org/berkeleylab/checkv/) | 1.0.1 | BSD-3-Clause-LBNL | `CHECKV_v1_0_1` (`end_to_end`). Modified-BSD variant used by Lawrence Berkeley National Laboratory. `params`-referenced CheckV DB is user-supplied. |
| [DeepVirFinder](https://github.com/jessieren/DeepVirFinder) | pinned by digest (`replikation/deepvirfinder@sha256:cc9666...`) | **USC-RL v1.0 — academic / non-commercial only** | `DEEPVIRFINDER` (`dvf.py`). Commercial use requires a separate paid license from the University of Southern California. Upstream ships no official container and has no numbered releases since ~2019 (the image content traces to upstream commit `475d883`, 2019-01-04). `replikation/deepvirfinder` has only ever published one tag (`:latest`, since 2019-05-13), so this pipeline pins its actual content digest instead of the mutable tag name — immutable even if that tag is ever reused. Image Dockerfile traces to the peer-reviewed "What the Phage" pipeline (`replikation/What_the_Phage`). |
| [HMMER](https://github.com/EddyRivasLab/hmmer) | 3.4 | BSD-3-Clause | `PROTEINS_VS_EGGNOG_4dot5` (`hmmsearch -E 0.00001` of each capsule's Prokka-predicted proteins against a user-supplied eggNOG HMM database). Distinct from the HMMER 3.3.2 that Prokka bundles internally, above. |

## Summary of obligations

- **Using these tools as external executables** (subprocess / CLI calls, passing
  data via files or stdout/stdin) is treated as "mere aggregation," not linking
  or derivative-work creation, so it does not require this pipeline's own code
  to adopt terms from tools under a different license.
- **Redistributing any of these tools' source or binaries** (e.g. in a
  Docker / Singularity image) requires preserving their original license and
  copyright notices; you may not relicense their code.
- **Modifying and redistributing** any GPL-licensed tool's source requires
  releasing those modifications under the same GPL version (or, where the tool
  is "-or-later," a compatible later version).
- **Non-commercial-only tools** — **geNomad** and **DeepVirFinder** restrict
  more than redistribution: they restrict *use* itself to academic /
  non-commercial contexts. Confirm eligibility before running this pipeline, or
  distributing its results, in any commercial setting. (Both are in the
  `--viral` block, which now **defaults to `false`** — only set `--viral true`
  if your use qualifies.)
- **Pheniqs** carries a repository `LICENSE` (restrictive NYU research license)
  that conflicts with its per-file AGPL-3.0-or-later headers, so this pipeline
  does not bundle or reference any Pheniqs container — downloading, installing,
  and running it is left entirely to the user's own discretion and risk (same
  treatment as KMERNORM, below).
- **KMERNORM** has no specified license, as of this writing — installation and use is at the user's own discretion.
- Reference **databases** (the Zenodo contaminant set, SILVA, GTDB, the Prokka
  SwissProt FASTA, and the geNomad / CheckV databases) are user-supplied and not
  redistributed by this repository; their own terms (SILVA CC BY 4.0, GTDB
  CC BY-SA 4.0, UniProt CC BY 4.0, etc.) still apply.

Last reviewed: September 2026.

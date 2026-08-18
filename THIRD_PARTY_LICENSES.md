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
must be preserved and included in the distributed image/artifact.

| Tool | Version | License | Notes |
|---|---|---|---|
| [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) | 0.11.9 | GPL-3.0-or-later | Bundles Picard SAM/BAM libraries and other third-party code under separate terms. |
| [Trimmomatic](https://github.com/usadellab/Trimmomatic) | 0.32 | GPL-3.0-or-later | Illumina adapter sequence FASTA files included with the distribution are **not** GPL-licensed — owned by and used with permission of Illumina, Inc. Handle separately if redistributing. |
| [KMERNORM](https://sourceforge.net/projects/kmernorm/) | 1.0.0 | **Unstated / unverified** | No license file found on the project page or in its source archive, as of this writing. This pipeline runs it via a user-provided binary rather than depending on any third-party container; installing and running it is entirely **at the discretion of the user**, who should independently confirm license terms and satisfy themselves of the tool's suitability before doing so — this pipeline does not warrant or vouch for it. See `README.md` for a known correctness issue with macOS-built copies. |
| [BWA](https://github.com/lh3/bwa) (Burrows-Wheeler Aligner) | — | GPL-3.0-or-later + MIT | GPLv3 governs the package as a whole (via BWT-SW-derived code); sorting, hash table, BWT, and IS libraries are separately MIT-licensed. |
| [SPAdes](https://github.com/ablab/spades) | 3.15.2 | GPL-2.0-only | |
| [CheckM](https://github.com/Ecogenomics/CheckM) | 1.1.9 | GPL-3.0-or-later | CheckM v1 is unmaintained upstream (superseded by CheckM2); this does not affect its license terms. |
| [Prokka](https://github.com/tseemann/prokka) | 1.14.6 | GPL-3.0-or-later | At runtime, invokes several other external annotation tools (BLAST+, HMMER, Aragorn, Infernal, minced, Prodigal, tbl2asn) as separate binaries within its container, each under its own license — not covered individually here. |
| [Biopython](https://biopython.org/) | 1.79 | Biopython License Agreement (permissive; some files dual-licensed BSD-3-Clause) | Used by `templates/prokka_gff_2_tsv.py` (`Bio.SeqIO`) for coding-density calculation. |
| [pandas](https://pandas.pydata.org/) | — | BSD-3-Clause | Version varies by container across processes (1.5.2 in the mulled image backing `PROKKA_GFF_2_TSV`, 2.2.1 elsewhere, e.g. `ASSEMBLY_STATS_TABULATOR`). |
| [NumPy](https://numpy.org/) | — | BSD-3-Clause | Bundled with pandas in the containers above; not directly imported by this pipeline's own scripts. |

## Summary of obligations

- **Using these tools as external executables** (subprocess/CLI calls,
  passing data via files or stdout/stdin) is treated as "mere aggregation,"
  not linking or derivative work creation, so it does not require this
  pipeline's own code to adopt terms from tools under a different license.
- **Redistributing any of these tools' source or binaries** (e.g., in a
  Docker/Singularity image) requires preserving their original license and
  copyright notices; you may not relicense their code as GPL-3.0-or-later
  or any other license.
- **Modifying and redistributing** any GPL-licensed tool's source requires
  releasing those modifications under the same GPL version (or, where the
  tool is "-or-later," a compatible later version).

Last reviewed: August 2026.

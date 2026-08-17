#!/usr/bin/env python
"""
Generates the tiny synthetic Atrandi-barcoded fastq pair CI runs against
(ATRANDI_TEST_R{1,2}.fastq.gz). Deterministic and committed alongside its
output so the fixture is reproducible/inspectable rather than a black box.

R2 layout mirrors the Pheniqs token windows the pipeline's demux configs use
(see PHENIQS_MAKE_SAMPLE_CONFIG/PHENIQS_NAME_CAPSULES in main.nf):
    [8bp D][4bp linker][8bp C][4bp linker][8bp B][4bp linker][8bp A][genomic filler]
Linker content (between the 8bp barcode windows) is arbitrary filler — Pheniqs
only reads the fixed token windows, so this isn't meant to be authentic
Atrandi linker sequence, just structurally correct for exercising the demux.

Three barcode combos are emitted, using real entries from barcodes/bc*_24.txt:
  - combo 1 (index 0 of each list): 10 read pairs -> clears read_threshold (3)
  - combo 2 (index 5 of each list): 5 read pairs  -> clears threshold, proves
    a second, distinct capsule resolves correctly
  - combo 3 (index 10, "noise"): 1 read pair -> below threshold, must be
    filtered out and land in the undetermined bucket
"""
import gzip
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
BARCODE_DIR = REPO_ROOT / "barcodes"
OUT_DIR = pathlib.Path(__file__).resolve().parent

READ_LENGTH = 100
LINKER = "GGGG"
GENOMIC_FILLER_R1 = ("ACGTACGTAGCTAGCTAGCATCGATCGATCGTAGCTAGCTAGCTAGCTAGCT" * 3)[:READ_LENGTH]


def load_barcodes(letter):
    lines = (BARCODE_DIR / f"bc{letter}_24.txt").read_text().splitlines()
    return [line.strip() for line in lines if line.strip()]


def build_r2_prefix(bc_d, bc_c, bc_b, bc_a):
    return bc_d + LINKER + bc_c + LINKER + bc_b + LINKER + bc_a


def main():
    bc_d, bc_c, bc_b, bc_a = (load_barcodes(letter) for letter in "DCBA")

    combos = [
        {"idx": 0, "n_pairs": 10, "name": "combo1"},
        {"idx": 5, "n_pairs": 5, "name": "combo2"},
        {"idx": 10, "n_pairs": 1, "name": "combo3_noise"},
    ]

    r1_records = []
    r2_records = []
    for combo in combos:
        i = combo["idx"]
        r2_prefix = build_r2_prefix(bc_d[i], bc_c[i], bc_b[i], bc_a[i])
        r2_seq = (r2_prefix + GENOMIC_FILLER_R1)[:READ_LENGTH]
        qual = "I" * READ_LENGTH
        for read_num in range(combo["n_pairs"]):
            read_id = f"ATRANDI_TEST:1:{combo['name']}:{read_num:04d}"
            r1_records.append((f"{read_id}/1", GENOMIC_FILLER_R1, qual))
            r2_records.append((f"{read_id}/2", r2_seq, qual))

    write_fastq(OUT_DIR / "ATRANDI_TEST_R1.fastq.gz", r1_records)
    write_fastq(OUT_DIR / "ATRANDI_TEST_R2.fastq.gz", r2_records)


def write_fastq(path, records):
    with gzip.open(path, "wt") as handle:
        for name, seq, qual in records:
            handle.write(f"@{name}\n{seq}\n+\n{qual}\n")


if __name__ == "__main__":
    main()

```bash
#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Filtering of unmapped reads
#
# Description:
#   Filters unmapped FASTQ reads based on sequence quality, nucleotide
#   composition, ambiguous bases, and repetitive sequence content.
#
# Usage:
#   bash 02_filter_unmapped.sh SAMPLE OUTDIR
#
# Input:
#   OUTDIR/02_unmapped_fastq/SAMPLE_unmapped_R1.fastq.gz
#   OUTDIR/02_unmapped_fastq/SAMPLE_unmapped_R2.fastq.gz
#
# Output:
#   OUTDIR/03_filtered_unmapped/SAMPLE_R1.filtered.fastq.gz
#   OUTDIR/03_filtered_unmapped/SAMPLE_R2.filtered.fastq.gz
#
# Requirements:
#   Python 3
# -----------------------------------------------------------------------------

SAMPLE="$1"
OUTDIR="$2"

IN_DIR="$OUTDIR/02_unmapped_fastq"
FILTERED_DIR="$OUTDIR/03_filtered_unmapped"

R1_IN="$IN_DIR/${SAMPLE}_unmapped_R1.fastq.gz"
R2_IN="$IN_DIR/${SAMPLE}_unmapped_R2.fastq.gz"

R1_OUT="$FILTERED_DIR/${SAMPLE}_R1.filtered.fastq.gz"
R2_OUT="$FILTERED_DIR/${SAMPLE}_R2.filtered.fastq.gz"

mkdir -p "$FILTERED_DIR"

filter_fastq() {

  INPUT="$1"
  OUTPUT="$2"

  python3 - "$INPUT" "$OUTPUT" <<'PY'
import gzip
import re
import sys
from itertools import islice

input_file = sys.argv[1]
output_file = sys.argv[2]

mono_run_re = re.compile(r"(A{21,}|C{21,}|G{21,}|T{21,})")

di_motifs = [
    a + b
    for a in "ACGT"
    for b in "ACGT"
    if a != b
] + [base * 2 for base in "ACGT"]

tri_motifs = [
    a + b + c
    for a in "ACGT"
    for b in "ACGT"
    for c in "ACGT"
]


def too_repetitive(sequence):

    if mono_run_re.search(sequence):
        return True

    for motif in di_motifs:
        if motif * 11 in sequence:
            return True

    for motif in tri_motifs:
        if motif * 11 in sequence:
            return True

    return False


def keep(record):

    header, sequence, separator, quality = record
    sequence = sequence.rstrip()
    quality = quality.rstrip()

    length = len(sequence)

    if length == 0:
        return False

    if sequence.count("N") / length > 0.10:
        return False

    quality_scores = [ord(character) - 33 for character in quality]

    if not quality_scores:
        return False

    trim_count = max(1, int(0.10 * len(quality_scores)))
    trimmed_scores = sorted(quality_scores)[trim_count:]

    if sum(trimmed_scores) / len(trimmed_scores) < 25:
        return False

    nucleotide_counts = [
        sequence.count("A"),
        sequence.count("C"),
        sequence.count("G"),
        sequence.count("T"),
    ]

    if max(nucleotide_counts) / length > 0.60:
        return False

    if min(nucleotide_counts) / length < 0.10:
        return False

    if too_repetitive(sequence):
        return False

    return True


with gzip.open(input_file, "rt") as input_handle, \
     gzip.open(output_file, "wt") as output_handle:

    records = iter(input_handle)

    while True:

        record = list(islice(records, 4))

        if not record:
            break

        if len(record) != 4:
            raise ValueError("Incomplete FASTQ record encountered")

        if keep(record):
            output_handle.writelines(record)
PY
}

echo "[$(date)] Filtering $SAMPLE R1"
filter_fastq "$R1_IN" "$R1_OUT"

echo "[$(date)] Filtering $SAMPLE R2"
filter_fastq "$R2_IN" "$R2_OUT"

echo "[$(date)] Done"
```


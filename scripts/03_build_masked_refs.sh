#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Build mismatch-masked reference genomes
#
# Description:
#   Generates one reference genome for each substitution class by replacing all
#   occurrences of a source nucleotide with a destination nucleotide, then
#   indexes each modified reference using BWA.
#
# Usage:
#   bash 03_build_masked_refs.sh REF.fa classes.tsv OUTDIR
#
# Input:
#   REF.fa
#   classes.tsv
#
# Output:
#   OUTDIR/REF_AG/
#   OUTDIR/REF_AC/
#   ...
#
# Requirements:
#   Python 3
#   BWA 0.7.18
# -----------------------------------------------------------------------------

REF="$1"
CLASSES="$2"
OUTDIR="$3"

mkdir -p "$OUTDIR"

while read LABEL SRC DST; do

    echo "[$(date)] Building ${LABEL} masked reference"

    REFDIR="$OUTDIR/REF_${LABEL}"
    mkdir -p "$REFDIR"

    python3 - "$REF" "$REFDIR/refined.${LABEL}.fa" "$SRC" "$DST" <<'PY'
import sys

inp, outp, src, dst = sys.argv[1:5]

translation = str.maketrans({
    src.upper(): dst.upper(),
    src.lower(): dst.lower()
})

with open(inp) as fin, open(outp, "w") as fout:
    for line in fin:
        if line.startswith(">"):
            fout.write(line)
        else:
            fout.write(line.translate(translation))
PY

    bwa index "$REFDIR/refined.${LABEL}.fa"

done < "$CLASSES"

echo "[$(date)] Done"

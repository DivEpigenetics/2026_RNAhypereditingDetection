```bash
#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Alignment against a mismatch-masked reference
#
# Description:
#   Converts the reads according to one substitution class and aligns them
#   against the corresponding mismatch-masked reference using BWA aln.
#   R1 and R2 are processed independently as single-end reads.
#
# Usage:
#   bash 04_masked_align.sh \
#     SAMPLE CLASS OUTDIR MASKED_REF_DIR classes.tsv
#
# Example:
#   bash 04_masked_align.sh \
#     SAMPLE01 AG results/SAMPLE01 masked_references \
#     config/classes.tsv
#
# Input:
#   OUTDIR/03_filtered_unmapped/SAMPLE_R1.filtered.fastq.gz
#   OUTDIR/03_filtered_unmapped/SAMPLE_R2.filtered.fastq.gz
#
# Requirements:
#   BWA 0.7.18
#   SAMtools 1.20
# -----------------------------------------------------------------------------

SAMPLE="$1"
CLASS="$2"
OUTDIR="$3"
MASKED_REF_DIR="$4"
CLASSES="$5"

THREADS=8

IN_DIR="$OUTDIR/03_filtered_unmapped"
ALIGN_DIR="$OUTDIR/04_masked_align/$CLASS"

REF="$MASKED_REF_DIR/REF_${CLASS}/refined.${CLASS}.fa"

SRC=$(awk -v class="$CLASS" '$1 == class {print $2}' "$CLASSES")
DST=$(awk -v class="$CLASS" '$1 == class {print $3}' "$CLASSES")

mkdir -p \
  "$ALIGN_DIR/sai" \
  "$ALIGN_DIR/bam"

echo "[$(date)] Masked alignment: $SAMPLE, class $CLASS"

for END in R1 R2; do

  FASTQ="$IN_DIR/${SAMPLE}_${END}.filtered.fastq.gz"
  SAI="$ALIGN_DIR/sai/${SAMPLE}_${END}.${CLASS}.sai"
  BAM="$ALIGN_DIR/bam/${SAMPLE}_${END}.${CLASS}.bam"

  echo "[$(date)] Processing $END"

  zcat "$FASTQ" \
    | awk -v SRC="$SRC" -v DST="$DST" '
        BEGIN {
          OFS = "\n"
        }

        /^@/ {
          header = $0
          getline sequence
          getline separator
          getline quality

          gsub(SRC, DST, sequence)
          gsub(tolower(SRC), tolower(DST), sequence)

          print header, sequence, "+", quality
        }
      ' \
    | bwa aln \
        -t "$THREADS" \
        -n 2 \
        -o 0 \
        -l 32 \
        "$REF" \
        - \
    > "$SAI"

  bwa samse \
      "$REF" \
      "$SAI" \
      <(zcat "$FASTQ") \
    | samtools view -bh - \
    > "$BAM"

done

echo "[$(date)] Done"
```


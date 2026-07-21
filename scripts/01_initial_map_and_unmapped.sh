# -----------------------------------------------------------------------------
# Initial alignment and extraction of unmapped read pairs
#
# Description:
#   Aligns paired-end RNA-seq reads to the reference genome using BWA-MEM,
#   sorts and indexes the resulting BAM file, and extracts read pairs for which
#   both mates are unmapped.
# Usage:
# bash 01_initial_map_and_unmapped.sh SAMPLE R1.fastq.gz R2.fastq.gz REF.fa OUTDIR

# Requirements:
#   BWA 0.7.18
#   SAMtools 1.20
# -----------------------------------------------------------------------------

SAMPLE="$1"
R1="$2"
R2="$3"
REF="$4"
OUTDIR="$5"

THREADS=8

mkdir -p \
  "$OUTDIR/01_initial_mapping" \
  "$OUTDIR/02_unmapped_fastq"

echo "[$(date)] Initial mapping: $SAMPLE"

bwa mem -t "$THREADS" "$REF" "$R1" "$R2" \
  | samtools view -bh - \
  | samtools sort -@ 7 -o "$OUTDIR/01_initial_mapping/${SAMPLE}.bam"

samtools index \
  "$OUTDIR/01_initial_mapping/${SAMPLE}.bam"

echo "[$(date)] Extracting unmapped read pairs"

samtools fastq \
  -f 12 \
  -F 256 \
  -@ 7 \
  -1 "$OUTDIR/02_unmapped_fastq/${SAMPLE}_unmapped_R1.fastq.gz" \
  -2 "$OUTDIR/02_unmapped_fastq/${SAMPLE}_unmapped_R2.fastq.gz" \
  -0 /dev/null \
  -s /dev/null \
  "$OUTDIR/01_initial_mapping/${SAMPLE}.bam"

echo "[$(date)] Done"

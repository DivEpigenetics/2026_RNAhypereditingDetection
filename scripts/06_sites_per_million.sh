```bash
#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Normalization of hyperediting sites
#
# Description:
#   Counts hyperediting sites for each substitution class and normalizes them
#   by the number of primary mapped reads and mapped query bases obtained from
#   the initial alignment.
#
# Usage:
#   bash 06_sites_per_million.sh \
#     SAMPLE OUTDIR classes.tsv
#
# Input:
#   OUTDIR/01_initial_mapping/SAMPLE.bam
#   OUTDIR/05_calls/CLASS/SAMPLE.CLASS.sites.bed
#
# Output:
#   OUTDIR/06_summary/SAMPLE.sites_per_million_reads_and_bases.tsv
#
# Requirements:
#   SAMtools 1.20
# -----------------------------------------------------------------------------

SAMPLE="$1"
OUTDIR="$2"
CLASSES="$3"

BAM="$OUTDIR/01_initial_mapping/${SAMPLE}.bam"
SUMMARY_DIR="$OUTDIR/06_summary"

mkdir -p "$SUMMARY_DIR"

OUTPUT_TSV="$SUMMARY_DIR/${SAMPLE}.sites_per_million_reads_and_bases.tsv"

echo -e \
  "sample\tclass\tsites\tmapped_primary_reads\tmapped_bases\tSPMM\tSPMB" \
  > "$OUTPUT_TSV"


# Count primary mapped reads.
# Excludes unmapped reads (0x4) and secondary alignments (0x100).
count_mapped_primary_reads() {

  samtools view -c -F 260 "$1"

}


# Sum mapped query bases from primary alignments.
# Includes M, =, X, and I CIGAR operations.
sum_mapped_bases() {

  samtools view -F 260 "$1" \
    | awk '
        function cigar_query_bases(cigar, total, token, length, operation) {

          total = 0

          while (match(cigar, /[0-9]+[MIDNSHP=X]/)) {

            token = substr(cigar, RSTART, RLENGTH)
            length = substr(token, 1, RLENGTH - 1) + 0
            operation = substr(token, RLENGTH, 1)

            if (
              operation == "M" ||
              operation == "=" ||
              operation == "X" ||
              operation == "I"
            ) {
              total += length
            }

            cigar = substr(cigar, RSTART + RLENGTH)
          }

          return total
        }

        {
          sum += cigar_query_bases($6)
        }

        END {
          printf "%.0f\n", sum
        }
      '

}


echo "[$(date)] Calculating normalization denominators: $SAMPLE"

MAPPED_READS=$(count_mapped_primary_reads "$BAM")
MAPPED_BASES=$(sum_mapped_bases "$BAM")


while read -r CLASS SRC DST; do

  [[ -z "$CLASS" || "$CLASS" == \#* ]] && continue

  SITES_BED="$OUTDIR/05_calls/$CLASS/${SAMPLE}.${CLASS}.sites.bed"

  if [[ -s "$SITES_BED" ]]; then
    SITES=$(wc -l < "$SITES_BED")
  else
    SITES=0
  fi

  SPMM=$(awk \
    -v sites="$SITES" \
    -v mapped_reads="$MAPPED_READS" \
    'BEGIN {
       printf "%.6f", sites / (mapped_reads / 1000000)
     }'
  )

  SPMB=$(awk \
    -v sites="$SITES" \
    -v mapped_bases="$MAPPED_BASES" \
    'BEGIN {
       printf "%.6f", sites / (mapped_bases / 1000000)
     }'
  )

  echo -e \
    "${SAMPLE}\t${CLASS}\t${SITES}\t${MAPPED_READS}\t${MAPPED_BASES}\t${SPMM}\t${SPMB}" \
    >> "$OUTPUT_TSV"

done < "$CLASSES"

echo "[$(date)] Done"
```


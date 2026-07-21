```bash
#!/usr/bin/env bash
set -euo pipefail

# -----------------------------------------------------------------------------
# Calling of hyperedited sites and clusters
#
# Description:
#   Calls hyperedited sites and clusters from the masked-reference alignments.
#   All substitution classes listed in classes.tsv are processed.
#
# Usage:
#   bash 05_call_and_cluster.sh \
#     SAMPLE OUTDIR classes.tsv call_hyperediting_sites.py
#
# Input:
#   OUTDIR/04_masked_align/CLASS/bam/SAMPLE_R1.CLASS.bam
#   OUTDIR/04_masked_align/CLASS/bam/SAMPLE_R2.CLASS.bam
#
# Requirements:
#   Python 3
#   SAMtools 1.20
# -----------------------------------------------------------------------------

SAMPLE="$1"
OUTDIR="$2"
CLASSES="$3"
CALL_SCRIPT="$4"

mkdir -p "$OUTDIR/05_calls"

while read -r CLASS SRC DST; do

  [[ -z "$CLASS" || "$CLASS" == \#* ]] && continue

  echo "[$(date)] Calling sites and clusters: $SAMPLE, class $CLASS"

  mkdir -p "$OUTDIR/05_calls/$CLASS"

  python3 "$CALL_SCRIPT" \
    "$OUTDIR" \
    "$SAMPLE" \
    "$CLASS"

done < "$CLASSES"

echo "[$(date)] Done"
```


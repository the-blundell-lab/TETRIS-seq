#!/bin/bash
###############################################################################
# fetch_reference.sh
#
# Downloads the b37 / GRCh37 reference genome the pipeline expects, and builds
# the BWA index if it is not supplied alongside it.
#
#   scripts/fetch_reference.sh [destination_directory]
#
# Default destination: $EXTERNAL_TOOLS (or ./Pipeline_tools if unset).
#
# Downloads ~3 GB. Building the BWA index takes around an hour and needs ~5 GB
# more; it is skipped if an index is already present.
###############################################################################

set -euo pipefail

DEST="${1:-${EXTERNAL_TOOLS:-./Pipeline_tools}}"
BASE="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0"
FASTA="Homo_sapiens_assembly19.fasta"

mkdir -p "$DEST"
cd "$DEST"

for f in "$FASTA" "$FASTA.fai" "Homo_sapiens_assembly19.dict"; do
    if [ -s "$f" ]; then
        echo "have    $f"
    else
        echo "fetch   $f"
        wget -c -O "$f" "$BASE/$f"
    fi
done

# sanity check: the FASTA should be ~3 GB and start with a FASTA header
[ "$(stat -c%s "$FASTA" 2>/dev/null || stat -f%z "$FASTA")" -gt 2000000000 ] || {
    echo "ERROR: $FASTA looks truncated - delete it and re-run." >&2; exit 1; }
head -c 1 "$FASTA" | grep -q '>' || {
    echo "ERROR: $FASTA does not start with a FASTA header." >&2; exit 1; }

if [ -s "$FASTA.bwt" ]; then
    echo "have    BWA index"
else
    command -v bwa >/dev/null 2>&1 || {
        echo "bwa not found - activate the conda environment first:" >&2
        echo "    conda activate tetris-seq-pipeline" >&2; exit 1; }
    echo "index   building the BWA index (about an hour)"
    bwa index "$FASTA"
fi

echo
echo "Reference ready in $DEST:"
ls -lh "$FASTA"* "Homo_sapiens_assembly19.dict" 2>/dev/null | awk '{printf "  %6s  %s\n", $5, $9}'
echo
echo "Set EXTERNAL_TOOLS to this directory in config/config.sh, then run scripts/check_setup.sh"

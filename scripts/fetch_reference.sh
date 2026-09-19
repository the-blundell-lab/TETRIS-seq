#!/bin/bash
###############################################################################
# fetch_reference.sh
#
# Obtains the b37 / GRCh37 reference genome the pipeline expects.
#
#   scripts/fetch_reference.sh [destination_directory]
#
# Default destination: $EXTERNAL_TOOLS (or ./Pipeline_tools if unset).
#
# By default it downloads a prebuilt copy (FASTA + .fai + .dict + BWA index,
# 4 GB compressed) from Zenodo: https://doi.org/10.5281/zenodo.22846473
#
# With --from-broad it instead downloads the FASTA from the Broad's public
# bucket and builds the BWA index locally (~3 GB download, ~1 hour of indexing).
###############################################################################

set -euo pipefail

SOURCE="zenodo"
DEST=""
for arg in "$@"; do
    case "$arg" in
        --from-broad) SOURCE="broad" ;;
        -h|--help) sed -n '3,20p' "$0"; exit 0 ;;
        *) DEST="$arg" ;;
    esac
done
DEST="${DEST:-${EXTERNAL_TOOLS:-./Pipeline_tools}}"

ZENODO_RECORD=22846473
TARBALL=Homo_sapiens_assembly19_b37_with_bwa_index.tar.gz
TARBALL_MD5=a6363075cf6905816e4a1d42c6c2eedf
BASE="https://storage.googleapis.com/gcp-public-data--broad-references/hg19/v0"
FASTA="Homo_sapiens_assembly19.fasta"

mkdir -p "$DEST"
cd "$DEST"

if [ -s "$FASTA.bwt" ] && [ -s "$FASTA" ]; then
    echo "have    an indexed reference already in $DEST - nothing to do"
elif [ "$SOURCE" = "zenodo" ]; then
    echo "fetch   $TARBALL from Zenodo record $ZENODO_RECORD (4 GB)"
    wget -c -O "$TARBALL" "https://zenodo.org/records/$ZENODO_RECORD/files/$TARBALL?download=1"

    echo "verify  md5"
    if command -v md5sum >/dev/null 2>&1; then got=$(md5sum "$TARBALL" | cut -d" " -f1)
    else got=$(md5 -q "$TARBALL"); fi
    [ "$got" = "$TARBALL_MD5" ] || { echo "ERROR: md5 mismatch ($got) - delete $TARBALL and retry." >&2; exit 1; }

    echo "unpack  $TARBALL"
    tar -xzf "$TARBALL"
    rm -f "$TARBALL"
else
    for f in "$FASTA" "$FASTA.fai" "Homo_sapiens_assembly19.dict"; do
        [ -s "$f" ] && echo "have    $f" || { echo "fetch   $f"; wget -c -O "$f" "$BASE/$f"; }
    done
    size=$(stat -c%s "$FASTA" 2>/dev/null || stat -f%z "$FASTA")
    [ "$size" -gt 2000000000 ] || { echo "ERROR: $FASTA looks truncated." >&2; exit 1; }
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

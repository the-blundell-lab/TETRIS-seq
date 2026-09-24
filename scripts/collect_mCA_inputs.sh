#!/bin/bash
###############################################################################
# collect_mCA_inputs.sh
#
# Gathers the two files per sample that the mCA callers need into one directory.
#
#   collect_mCA_inputs.sh <search_root> <output_directory> [--copy]
#
# The callers (mCA_caller/watson_code_mCA_caller_unphased.ipynb and
# _phased.ipynb) read everything from a single flat directory, set as
# CNV_DEPOSIT_DIR at the top of each notebook. The pipeline, by contrast, leaves
# its output nested per sample:
#
#   <sample>_<UDI>/..._variant_calling_only_SNPs_annovar_annotated.txt   (CNV panel)
#   <library>/<sample>/PON_normalised_log2ratios_Feb2026/
#                          <sample>_PON_normalised_read_depths_and_LRR.txt  (PON step)
#
# This script searches <search_root> for both kinds and symlinks them into
# <output_directory>, which you then point CNV_DEPOSIT_DIR at. Use --copy for
# real copies instead of links.
###############################################################################

set -euo pipefail

[ $# -ge 2 ] || { sed -n '3,22p' "$0"; exit 1; }

ROOT="$(cd "$1" && pwd)"
OUT="$2"
MODE="link"
[ "${3:-}" = "--copy" ] && MODE="copy"

mkdir -p "$OUT"
OUT="$(cd "$OUT" && pwd)"

LRR_SUFFIX="_PON_normalised_read_depths_and_LRR.txt"
SNP_SUFFIX="variant_calling_only_SNPs_annovar_annotated.txt"

n_lrr=0; n_snp=0; skipped=0
while IFS= read -r f; do
    base="$(basename "$f")"
    target="$OUT/$base"
    if [ -e "$target" ]; then skipped=$((skipped+1)); continue; fi
    if [ "$MODE" = "copy" ]; then cp "$f" "$target"; else ln -s "$f" "$target"; fi
    case "$base" in
        *"$LRR_SUFFIX") n_lrr=$((n_lrr+1)) ;;
        *"$SNP_SUFFIX") n_snp=$((n_snp+1)) ;;
    esac
done < <(find "$ROOT" \( -name "*$LRR_SUFFIX" -o -name "*$SNP_SUFFIX" \) -type f | sort)

echo "searched   $ROOT"
echo "output     $OUT  ($MODE)"
echo "LRR files  $n_lrr"
echo "SNP files  $n_snp"
[ "$skipped" -gt 0 ] && echo "skipped    $skipped already present"

if [ "$n_lrr" -eq 0 ] || [ "$n_snp" -eq 0 ]; then
    echo
    echo "Warning: the callers need BOTH a *${LRR_SUFFIX} and a"
    echo "         *${SNP_SUFFIX} for each sample."
    [ "$n_lrr" -eq 0 ] && echo "         No LRR files found - has the PON notebook been run (step 1)?"
    [ "$n_snp" -eq 0 ] && echo "         No SNP files found - these come from the CNV panel."
fi

echo
echo "Now set this at the top of the caller notebooks:"
echo "    CNV_DEPOSIT_DIR = '$OUT'"

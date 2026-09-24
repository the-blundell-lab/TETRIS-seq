#!/bin/bash
###############################################################################
# rearrangement_SSCS_recall.sh
#
# Second pass of chromosomal rearrangement calling: once pass 1 has found a
# rearrangement in the raw mapped merged BAM, re-call it on consensus (SSCS)
# reads over the breakpoint regions only, to get an accurate VAF.
#
# It chains the three steps that would otherwise be run by hand:
#   1. SSCS consensus calling restricted to the breakpoint regions
#      (reads the SAME raw mapped merged BAM as pass 1 - it regroups by UMI)
#   2. samtools sort + index of that SSCS BAM
#   3. the translocation caller again, this time on the sorted SSCS BAM
#
# Usage:
#   rearrangement_SSCS_recall.sh -i <mapped_merged.bam> -s <sample_name> \
#       -o <output_directory> -r "<chrom start end [chrom start end ...]>"
#
#   -i  raw mapped merged BAM from the CNV / rearrangement panel (as pass 1)
#   -s  sample name; outputs are prefixed <sample_name>_SSCS
#   -o  output directory (created if needed)
#   -r  breakpoint regions, quoted: chrom start end triples, e.g.
#         -r "21 36210000 36215000 8 93078000 93080000"
#       Take these from the LEFT COORDINATE and RIGHT COORDINATE of the call in
#       pass 1's *_grouped_and_filtered.csv, with a window around each.
#   -l  read length (default: 146)
#   -q  minimum mapping quality (default: 20)
#
# Consensus settings are the region-restricted caller's defaults, as used in the
# manuscript. Reference genome, panel BED and ideogram come from the pipeline
# configuration, as for every other script: set EXTERNAL_TOOLS (or REF) in
# config/config.sh, and ANNOVAR_HOME or EXTERNAL_TOOLS for the annotation step.
###############################################################################

set -euo pipefail

P="$(cd "$(dirname "$0")" && pwd)"
PIPELINE_TOOLS="${PIPELINE_TOOLS:-$P/../pipeline_tools}"
CALLER_DIR="$P/../chromosomal_rearrangement_caller"

[ -f "$P/../config/config.sh" ] && . "$P/../config/config.sh"

: "${EXTERNAL_TOOLS:=$HOME/Pipeline_tools}"
: "${REF:=$EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta}"
: "${ANNOVAR_HOME:=$EXTERNAL_TOOLS/annovar}"
export ANNOVAR_HOME EXTERNAL_TOOLS

CNV_BED="$PIPELINE_TOOLS/TWIST_CNV_panel_TE-95031423_h19.bed"
TARGET_BED="$CALLER_DIR/Translocation_regions_of_interest.bed"
IDEOGRAM="$PIPELINE_TOOLS/chromosome_ideogram_hg19.txt"

SSCS_SCRIPT="$CALLER_DIR/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1_specific_regions.py"
CALLER="$CALLER_DIR/Watson_code_translocation_calling_all_types_2025_v1_targeted.py"

# refuse to start under a Python whose shelve would fall back to dbm.dumb
. "$P/require_python_backend.sh"

INFILE=""; SAMPLE=""; OUTDIR=""; REGIONS=""; READLEN=146; MAPQ=20

usage() { sed -n '3,32p' "$0"; }

while getopts "i:s:o:r:l:q:h" opt; do
  case "$opt" in
    i) INFILE="$OPTARG" ;;
    s) SAMPLE="$OPTARG" ;;
    o) OUTDIR="$OPTARG" ;;
    r) REGIONS="$OPTARG" ;;
    l) READLEN="$OPTARG" ;;
    q) MAPQ="$OPTARG" ;;
    h) usage; exit 0 ;;
    *) usage; exit 1 ;;
  esac
done

for pair in "INFILE:i" "SAMPLE:s" "OUTDIR:o" "REGIONS:r"; do
    v="${pair%%:*}"; flag="${pair##*:}"
    if [ -z "${!v}" ]; then echo "Error: -$flag is required"; echo; usage; exit 1; fi
done
[ -f "$INFILE" ] || { echo "Error: BAM not found: $INFILE"; exit 1; }
[ -f "$REF" ]    || { echo "Error: reference not found: $REF (set EXTERNAL_TOOLS or REF)"; exit 1; }

# regions must be whole triples
set -- $REGIONS
if [ $(( $# % 3 )) -ne 0 ] || [ $# -eq 0 ]; then
    echo "Error: -r must be chrom start end triples, e.g. \"21 36210000 36215000 8 93078000 93080000\""
    exit 1
fi
echo "regions: $(( $# / 3 ))"

# the consensus caller writes its read-distribution metrics into this subdirectory
mkdir -p "$OUTDIR/Metrics_and_images" TEMP

SSCS_BAM="$OUTDIR/${SAMPLE}_SSCS_specific_regions.bam"
SSCS_UNPAIRED="$OUTDIR/${SAMPLE}_SSCS_specific_regions_unpaired.bam"
SSCS_SORTED="$OUTDIR/${SAMPLE}_SSCS_specific_regions_sorted.bam"

echo
echo ">> 1. SSCS consensus calling over the breakpoint regions"
# shellcheck disable=SC2086
python "$SSCS_SCRIPT" \
    --infile "$INFILE" \
    --sample-name "${SAMPLE}_SSCS" \
    --min-mapping-quality "$MAPQ" \
    --regions $REGIONS \
    --outbam "$SSCS_BAM" \
    --unpaired-outbam "$SSCS_UNPAIRED" \
    --out-directory "$OUTDIR"

echo
echo ">> 2. Sorting and indexing the SSCS BAM"
samtools sort -o "$SSCS_SORTED" "$SSCS_BAM"
samtools index "$SSCS_SORTED"

echo
echo ">> 3. Re-calling rearrangements on the SSCS BAM"
python "$CALLER" \
    --infile "$SSCS_SORTED" \
    --sample-name "${SAMPLE}_SSCS" \
    --min-mapping-quality "$MAPQ" \
    --read-length "$READLEN" \
    --bed "$CNV_BED" \
    --targeted_bed "$TARGET_BED" \
    --min-reads 5 \
    --min-softclip-length 10 \
    --chromosomal_ideogram "$IDEOGRAM" \
    --ref "$REF" \
    --out-directory "$OUTDIR"

rm -f TEMP/${SAMPLE}_SSCS_*

echo
echo "Done. VAFs for the rearrangement are in:"
echo "  $OUTDIR/${SAMPLE}_SSCS_translocations_found_just_those_specifically_targeted_both_sides_panel_grouped_and_filtered.csv"

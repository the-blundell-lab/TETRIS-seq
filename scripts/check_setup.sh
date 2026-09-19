#!/bin/bash
###############################################################################
# check_setup.sh
#
# Preflight check: reports which tools, JARs, reference files and panel
# resources the pipeline can find, using the same resolution rules as
# Watson_code_environment_setup_v2.6.sh (config/config.sh, then defaults).
#
# Run this before processing any samples:
#
#   scripts/check_setup.sh
#
# Exit status is 0 if everything required was found, 1 otherwise.
###############################################################################

P="$(cd "$(dirname "$0")" && pwd)"
CONFIG_FILE="${TETRIS_CONFIG:-$P/../config/config.sh}"
[ -f "$CONFIG_FILE" ] && source "$CONFIG_FILE" && echo "config: $CONFIG_FILE" \
                                               || echo "config: none found (using defaults)"

: "${EXTERNAL_TOOLS:=$HOME/Pipeline_tools}"
: "${REF:=$EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta}"
: "${VARDICT_HOME:=$HOME/VarDictJava}"
: "${SNPEFF_HOME:=$HOME/snpEff}"
: "${ANNOVAR_HOME:=$EXTERNAL_TOOLS/annovar}"
: "${PINDEL_DIR:=$EXTERNAL_TOOLS/pindel}"
: "${PIPELINE_TOOLS:=$P/../pipeline_tools}"

missing=0
ok()   { printf '  [ ok ]      %-34s %s\n' "$1" "$2"; }
bad()  { printf '  [MISSING]   %-34s %s\n' "$1" "$2"; missing=$((missing+1)); }
note() { printf '  [ note ]    %-34s %s\n' "$1" "$2"; }

check_cmd()  { if command -v "$1" >/dev/null 2>&1; then ok "$1" "$(command -v "$1")"; else bad "$1" "not on PATH"; fi; }
check_file() { if [ -f "$2" ]; then ok "$1" "$2"; else bad "$1" "$2"; fi; }
check_dir()  { if [ -d "$2" ]; then ok "$1" "$2"; else bad "$1" "$2"; fi; }

echo
echo "Command-line tools"
for c in bwa samtools tabix; do check_cmd "$c"; done

if command -v python >/dev/null 2>&1; then ok "python" "$(command -v python) ($(python --version 2>&1))"
elif command -v python3 >/dev/null 2>&1; then ok "python3" "$(command -v python3) ($(python3 --version 2>&1))"
else bad "python" "not on PATH"; fi

if command -v java >/dev/null 2>&1; then
    jv=$(java -version 2>&1 | head -1)
    case "$jv" in
        *\"1.8*|*\"8.*) ok "java" "$jv" ;;
        *) note "java" "$jv - Picard 2.18, GATK 3.8 and fgbio 1.3 expect Java 8" ;;
    esac
else bad "java" "not on PATH"; fi

echo
echo "Java archives (\$EXTERNAL_TOOLS = $EXTERNAL_TOOLS)"
check_file "picard.jar"            "$EXTERNAL_TOOLS/picard.jar"
check_file "fgbio-1.3.0.jar"       "$EXTERNAL_TOOLS/fgbio-1.3.0.jar"
check_file "GenomeAnalysisTK.jar"  "$EXTERNAL_TOOLS/GenomeAnalysisTK.jar"

echo
echo "Other third-party tools"
check_file "VarDict"               "$VARDICT_HOME/build/install/VarDict/bin/VarDict"
check_file "snpEff.jar"            "$SNPEFF_HOME/snpEff.jar"
check_file "SnpSift.jar"           "$SNPEFF_HOME/SnpSift.jar"
check_file "ANNOVAR table_annovar" "$ANNOVAR_HOME/table_annovar.pl"
check_dir  "ANNOVAR humandb/"      "$ANNOVAR_HOME/humandb"
check_file "pindel"                "$PINDEL_DIR/pindel"
check_file "pindel2vcf"            "$PINDEL_DIR/pindel2vcf"

echo
echo "Reference genome"
check_file "reference FASTA"       "$REF"
check_file "  .fai index"          "${REF}.fai"
check_file "  .dict"               "${REF%.fasta}.dict"
check_file "  BWA index (.bwt)"    "${REF}.bwt"
check_dir  "dbSNP intervals"       "$EXTERNAL_TOOLS/dbSNP"

echo
echo "Panel resources and custom scripts (\$PIPELINE_TOOLS = $PIPELINE_TOOLS)"
for f in Watson_code_SSCS_calling_2.1.py \
         Watson_code_DCS_calling_1.4.py \
         Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py \
         TWIST_SNV_panel-TE-92048328_hg19.bed \
         TWIST_CNV_panel_TE-95031423_h19.bed \
         TWIST_v2_transcript_details.csv \
         chromosome_ideogram_hg19.txt; do
    check_file "$f" "$PIPELINE_TOOLS/$f"
done

echo
if [ "$missing" -eq 0 ]; then
    echo "All checks passed."
else
    echo "$missing item(s) missing - see docs/INSTALL.md."
fi
exit $(( missing > 0 ))

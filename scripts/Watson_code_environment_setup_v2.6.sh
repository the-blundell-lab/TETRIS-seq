# Initializes the pipeline environment so that scripts can run as expected
function initialize() {
    IFS=$'\n\t'
    detect_platform

    # --------------------------------------------------------------------------
    # User-configurable install locations.
    # Copy config/config.sh.example to config/config.sh and edit it, OR export
    # these variables in your shell, to point the pipeline at your own installs.
    # If nothing is set, the original defaults ($HOME/...) are used unchanged.
    # --------------------------------------------------------------------------
    CONFIG_FILE="${TETRIS_CONFIG:-$P/../config/config.sh}"
    [ -f "$CONFIG_FILE" ] && source "$CONFIG_FILE"
    : "${PIPELINE_TOOLS:=$P/../pipeline_tools}"                 # ships with the repo: custom scripts + panel BED/CSV resources
    : "${EXTERNAL_TOOLS:=$HOME/Pipeline_tools}"                 # you install: jars, reference genome, ANNOVAR, Pindel, dbSNP
    : "${REF:=$EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta}"  # Broad b37 / GRCh37 reference FASTA
    : "${ANNOVAR_HOME:=$EXTERNAL_TOOLS/annovar}"               # ANNOVAR install (annotate_variation.pl etc.)
    : "${PINDEL_DIR:=$EXTERNAL_TOOLS/pindel}"                  # Pindel install

    # bwa/tabix: use the copies on PATH (e.g. from the conda environment) if present,
    # otherwise fall back to binaries bundled under scripts/bin/<platform>/.
    bwa=${BWA:-$(command -v bwa || echo "$P/bin/$PLATFORM/bwa")}
    tabix=${TABIX:-$(command -v tabix || echo "$P/bin/$PLATFORM/tabix")}
    # JARs: an explicit setting wins, then $EXTERNAL_TOOLS, then the active conda
    # environment (conda installs them under $CONDA_PREFIX/share/<pkg>/).
    find_jar() {   # $1 = file name, $2 = conda share glob
        if   [ -f "$EXTERNAL_TOOLS/$1" ]; then echo "$EXTERNAL_TOOLS/$1"
        elif [ -n "${CONDA_PREFIX:-}" ] && [ -f "$(echo $CONDA_PREFIX/share/$2 | cut -d" " -f1)" ]; then
             echo $CONDA_PREFIX/share/$2 | cut -d" " -f1
        else echo "$EXTERNAL_TOOLS/$1"
        fi
    }
    picard=${PICARD_JAR:-$(find_jar picard.jar "picard*/picard.jar")}
    fgbio=${FGBIO_JAR:-$(find_jar fgbio-1.3.0.jar "fgbio*/fgbio.jar")}
    GATK=${GATK_JAR:-$EXTERNAL_TOOLS/GenomeAnalysisTK.jar}
    chromosome_ideogram=$PIPELINE_TOOLS/chromosome_ideogram_hg19.txt
    dbSNP_directory=$EXTERNAL_TOOLS/dbSNP
    watson_call_SSCS=$PIPELINE_TOOLS/Watson_code_SSCS_calling_2.1.py
    watson_call_DCS=$PIPELINE_TOOLS/Watson_code_DCS_calling_1.4.py
    watson_VCF_text=$PIPELINE_TOOLS/Watson_code_VCF_SNP_calling_v1.3.py
    watson_filter_SSCS=$PIPELINE_TOOLS/Watson_code_filter_SSCS_1.1.py
    watson_VarDict_to_text=$PIPELINE_TOOLS/Watson_code_VarDict_to_TXT.py
    watson_VarDict_annotate=$PIPELINE_TOOLS/Watson_code_VarDict_annotation_v1.2.py
    watson_code_sample=$PIPELINE_TOOLS/Watson_code_sample_name_from_UDI_index_v1.py
    watson_code_sample_UDI=$PIPELINE_TOOLS/Watson_code_sample_name_from_UDI_index_with_UDI_v1.py
    watson_code_library_name=$PIPELINE_TOOLS/Watson_code_library_name_v1.py
    REF=${REF:-$EXTERNAL_TOOLS/Homo_sapiens_assembly19.fasta}

    watson_code_unmap_bam=$PIPELINE_TOOLS/Watson_code_convert_mapped_to_unmapped_BAM_v1.py

    #specific to SNV panel
    watson_read_depths=$PIPELINE_TOOLS/Watson_code_SNV_panel_read_coverage.py
    watson_VCF_text_SNV_panel=$PIPELINE_TOOLS/Watson_code_VCF_SNP_calling_SNV_panel_v1.2.py #contains quality scores
    watson_VCF_text_SNV_panel_duplex=$PIPELINE_TOOLS/Watson_code_duplex_VCF_SNP_calling_SNV_panel_v1.2.py #replaces quality scores with total UMI family sizes
    ANNOVAR=$ANNOVAR_HOME/annotate_variation.pl
    ANNOVAR_annotate=$ANNOVAR_HOME/table_annovar.pl
    ANNOVAR_humandb=$ANNOVAR_HOME/humandb/
    annotate_DCS=$PIPELINE_TOOLS/Watson_code_DCS_annotation_v1.2.py
    annotate_SSCS=$PIPELINE_TOOLS/Watson_code_SSCS_annotation_v1.2.py
    TWIST_transcript_details=$PIPELINE_TOOLS/TWIST_v2_transcript_details.csv
    SNV_bed=$PIPELINE_TOOLS/TWIST_SNV_panel-TE-92048328_hg19.bed
    SNV_bed_GATK=$PIPELINE_TOOLS/TWIST_SNV_panel_TE-92048328_GATK_compatible.bed

    #FLT3 calling (SNV panel)
    watson_call_DCS_FLT3=$PIPELINE_TOOLS/Watson_code_DCS_calling_1.4_for_FLT3_calling.py
    pindel_config_file=$PIPELINE_TOOLS/Watson_code_FLT3_create_pindel_config_file.py
    # pindel: $PINDEL_DIR if you installed it yourself, otherwise the conda environment
    pindel=$( [ -x "$PINDEL_DIR/pindel" ] && echo "$PINDEL_DIR/pindel" || command -v pindel || echo "$PINDEL_DIR/pindel" )
    pindel2vcf=$( [ -x "$PINDEL_DIR/pindel2vcf" ] && echo "$PINDEL_DIR/pindel2vcf" || command -v pindel2vcf || echo "$PINDEL_DIR/pindel2vcf" )
    watson_call_FLT3_ITD=$PIPELINE_TOOLS/Watson_code_FLT3_ITD_calling_from_pindel_vcf_v1.py

    #specific to CNV panel
    watson_SSCS_depths=$PIPELINE_TOOLS/Watson_code_sample_normalised_read_depths_v1.8.py
    KMT2A_coordinates=$PIPELINE_TOOLS/KMT2A_coordinates_hg19.csv
    watson_KMT2A_PTD=$PIPELINE_TOOLS/Watson_code_KMT2A-PTD_calling_v1.9.py
    # PON-normalised LRR for mCA calling is now produced by
    # mCA_caller/watson_code_create_PON.ipynb (supersedes the old
    # Watson_code_PON_normalised_read_depths_v1.1.py, which is no longer used).
    watson_translocation_breakpoint_depths=$PIPELINE_TOOLS/Watson_code_translocation_breakpoint_coverage_v1.2.py
    watson_SNP_plots=$PIPELINE_TOOLS/Watson_code_CNV_SNP_plotting_v1.2.py
    targeted_SNPs=$PIPELINE_TOOLS/TWIST_CNV_panel_targeted_SNPs_positions.csv
    watson_phased_SNP_plots=$PIPELINE_TOOLS/Watson_code_CNV_phased_SNP_plotting_v1.1.py
    CNV_bed=$PIPELINE_TOOLS/TWIST_CNV_panel_TE-95031423_h19.bed
    CNV_bed_GATK=$PIPELINE_TOOLS/TWIST_CNV_panel_TE-95031423_GATK_compatible.bed

    #Translocation_calling
    # One consensus caller serves both the FLT3-ITD arm and the rearrangement arm:
    # it is the standard SSCS caller with the SAM flag filter removed. The
    # rearrangement arm uses the --regions variant in chromosomal_rearrangement_caller/.
    watson_call_SSCS_translocations_or_FLT3=$PIPELINE_TOOLS/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py

    check_dependencies
}

# Function to check and ensure that the dependencies necessary to run the pipeline
# are available
function check_dependencies() {
    # Check that an appropriate version of Java is available
    if ! type -p java > /dev/null ; then
        fail "Could not find java executble. Please ensure Java 1.8 or above is on the path." \
             "To download the latest version go to: https://java.com/en/download/manual.jsp"
    fi

    java_version=$(java -version 2>&1 | fgrep version | cut -d\  -f3 | tr -d '"' | cut -d. -f1-2)
    major=$(echo $java_version | cut -d. -f1)
    minor=$(echo $java_version | cut -d. -f2)

    if [ $major -lt 2 ] && [ $minor -lt 8 ]; then
        fail "Detected java version $java_version on the path. Java version 1.8+ is required" \
             "to run this pipeline. To download the latest version go to: " \
             "https://java.com/en/download/manual.jsp"
    fi

    # TODO: CHECK THAT R is available
}

# Function to detect that platform being executed on and set PLATFORM
function detect_platform() {
  if   $(uname -a | fgrep -i darwin > /dev/null); then PLATFORM="mac"
  elif $(uname -a | fgrep -i linux  > /dev/null); then PLATFORM="linux"
  else fail "Could not detect supported operating system."
  fi
}

#
# Function used to exit after printing a large error message
function fail() {
    banner $*
    exit 1
}

function execute() {
    log "--------------------------------------------------------------------------------"
    log "- Executing: "
    for i in `seq 1 $#`; do
        if [ $i -eq $# ]; then lineterm=""; else lineterm=" \\"; fi
        if [ $i -eq 1 ];  then prefix=""; else prefix="    "; fi
        log "- ${prefix}${!i}${lineterm}"
    done
    log "--------------------------------------------------------------------------------"

    OLD=$IFS
    IFS=" "
    command="$@"
    eval $command
    IFS=$OLD
}

# Logs a message to the console along with the date an time
function log() {
    echo [`date`] $*
}

# Simple function to take arguments and return them as a file-system-safe string
function make_fss() {
    echo $* | tr '!$#()[]' '-'
}

# A short function for echoing a string to the screen in a banner
function banner() {
    echo
    echo "################################################################################"
    for line in "$@"; do
        echo "# $line"
    done
    echo "################################################################################"
    echo
}

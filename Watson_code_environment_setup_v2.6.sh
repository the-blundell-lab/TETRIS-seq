# Initializes the pipeline environment so that scripts can run as expected
function initialize() {
    IFS=$'\n\t'
    detect_platform
    bwa=$P/bin/$PLATFORM/bwa
    tabix=$P/bin/$PLATFORM/tabix
    picard=$HOME/Pipeline_tools/picard.jar
    fgbio=$HOME/Pipeline_tools/fgbio-1.3.0.jar
    GATK=$HOME/Pipeline_tools/GenomeAnalysisTK.jar
    VarDict=$HOME/VarDictJava/build/install/VarDict/bin/VarDict
    Teststrandbias=$HOME/VarDictJava/build/install/VarDict/bin/teststrandbias.R
    Var2VCF=$HOME/VarDictJava/build/install/VarDict/bin/var2vcf_valid.pl
    SNPEff=$HOME/snpEff/snpEff.jar
    SNPSift_one_line=$HOME/snpEff/scripts/vcfEffOnePerLine.pl
    SNPSift=$HOME/snpEff/SnpSift.jar
    chromosome_ideogram=$HOME/Pipeline_tools/chromosome_ideogram_hg19.txt
    dbSNP_directory=$HOME/Pipeline_tools/dbSNP
    watson_call_SSCS=$HOME/Pipeline_tools/Watson_code_SSCS_calling_2.1.py
    watson_call_DCS=$HOME/Pipeline_tools/Watson_code_DCS_calling_1.4.py
    watson_VCF_text=$HOME/Pipeline_tools/Watson_code_VCF_SNP_calling_v1.3.py
    watson_filter_SSCS=$HOME/Pipeline_tools/Watson_code_filter_SSCS_1.1.py
    watson_VarDict_to_text=$HOME/Pipeline_tools/Watson_code_VarDict_to_TXT.py
    watson_VarDict_annotate=$HOME/Pipeline_tools/Watson_code_VarDict_annotation_v1.2.py
    watson_code_sample=$HOME/Pipeline_tools/Watson_code_sample_name_from_UDI_index_v1.py
    watson_code_sample_UDI=$HOME/Pipeline_tools/Watson_code_sample_name_from_UDI_index_with_UDI_v1.py
    watson_code_library_name=$HOME/Pipeline_tools/Watson_code_library_name_v1.py
    REF=$HOME/Pipeline_tools/Homo_sapiens_assembly19.fasta

    watson_code_unmap_bam=$HOME/Pipeline_tools/Watson_code_convert_mapped_to_unmapped_BAM_v1.py

    #specific to SNV panel
    watson_read_depths=$HOME/Pipeline_tools/Watson_code_SNV_panel_read_coverage.py
    watson_VCF_text_SNV_panel=$HOME/Pipeline_tools/Watson_code_VCF_SNP_calling_SNV_panel_v1.2.py #contains quality scores
    watson_VCF_text_SNV_panel_duplex=$HOME/Pipeline_tools/Watson_code_duplex_VCF_SNP_calling_SNV_panel_v1.2.py #replaces quality scores with total UMI family sizes
    ANNOVAR=$HOME/Pipeline_tools/annovar/annotate_variation.pl
    ANNOVAR_annotate=$HOME/Pipeline_tools/annovar/table_annovar.pl
    ANNOVAR_humandb=$HOME/Pipeline_tools/annovar/humandb/
    annotate_DCS=$HOME/Pipeline_tools/Watson_code_DCS_annotation_v1.2.py
    annotate_SSCS=$HOME/Pipeline_tools/Watson_code_SSCS_annotation_v1.2.py
    TWIST_transcript_details=$HOME/Pipeline_tools/TWIST_v2_transcript_details.csv
    SNV_bed=$HOME/Pipeline_tools/TWIST_SNV_panel-TE-92048328_hg19.bed
    SNV_bed_GATK=$HOME/Pipeline_tools/TWIST_SNV_panel_TE-92048328_GATK_compatible.bed

    #FLT3 calling (SNV panel)
    watson_call_DCS_FLT3=$HOME/Pipeline_tools/Watson_code_DCS_calling_1.4_for_FLT3_calling.py
    pindel_config_file=$HOME/Pipeline_tools/Watson_code_FLT3_create_pindel_config_file.py
    pindel=$HOME/Pipeline_tools/pindel/pindel
    pindel2vcf=$HOME/Pipeline_tools/pindel/pindel2vcf
    watson_call_FLT3_ITD=$HOME/Pipeline_tools/Watson_code_FLT3_ITD_calling_from_pindel_vcf_v1.py

    #specific to CNV panel
    watson_SSCS_depths=$HOME/Pipeline_tools/Watson_code_sample_normalised_read_depths_v1.8.py
    KMT2A_coordinates=$HOME/Pipeline_tools/KMT2A_coordinates_hg19.csv
    watson_KMT2A_PTD=$HOME/Pipeline_tools/Watson_code_KMT2A-PTD_calling_v1.9.py
    watson_PON_normalised_LRR=$HOME/Pipeline_tools/Watson_code_PON_normalised_read_depths_v1.1.py
    watson_translocation_breakpoint_depths=$HOME/Pipeline_tools/Watson_code_translocation_breakpoint_coverage_v1.2.py
    watson_SNP_plots=$HOME/Pipeline_tools/Watson_code_CNV_SNP_plotting_v1.2.py
    targeted_SNPs=$HOME/Pipeline_tools/TWIST_CNV_panel_targeted_SNPs_positions.csv
    watson_phased_SNP_plots=$HOME/Pipeline_tools/Watson_code_CNV_phased_SNP_plotting_v1.2.py
    CNV_bed=$HOME/Pipeline_tools/TWIST_CNV_panel_TE-95031423_h19.bed
    CNV_bed_GATK=$HOME/Pipeline_tools/TWIST_CNV_panel_TE-95031423_GATK_compatible.bed

    #Translocation_calling
    watson_call_SSCS_translocations=$HOME/Pipeline_tools/Watson_code_SSCS_calling_for_translocations_1.1.py
    watson_call_SSCS_translocations_or_FLT3=$HOME/Pipeline_tools/Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py

    MANTA=$HOME/Pipeline_tools/manta-1.6.0.centos6_x86_64/bin/configManta.py
    FILTER_MANTA=$HOME/Pipeline_tools/Filter_Manta_output.py

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

#!/bin/bash
##################################################################################

set -euo pipefail #this command makes sure the script will be stopped when the first error is encountered.

################################################################################
# Argument processing and setup
################################################################################

#DEFAULTS
SAMPLE=''
sample_name_UDI=''
THREADS=4
CALLMIN=3
INSERTSIZE=300
SISYPHUS_OUTFOLDER=''
MINDEPTH=500

usage() {
  cat << EOF
Usage: $0 [options]

Analyzes mapped merged BAM file.
OPTIONS:
   -s <SAMPLE> Sample name
   -u <sample_name_UDI> Sample name with UDI
   -i <insert_size> The insert size for the sequencing output (default: $INSERTSIZE)
   -m <SSCS-min-reads> The minimum number of reads to form a SSCS read (default: $CALLMIN)
   -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
   -o <sisyphus-outfolder> directory on sisyphus in which to write the files
   -d <MINDEPTH> minimum total depth at positions to look at (default:$MINDEPTH)
EOF
}

while getopts “s:u:i:m:t:o:d:” OPTION; do
  case $OPTION in
    s) SAMPLE=$OPTARG;;
    u) sample_name_UDI=$OPTARG;;
    i) INSERTSIZE=$OPTARG;;
    m) CALLMIN=$OPTARG;;
    t) THREADS=$OPTARG;;
    o) SISYPHUS_OUTFOLDER=$OPTARG;;
    d) MINDEPTH=$OPTARG;;
    h) echo "Unknown option ${OPTION}"; usage; exit;;
    [?]) usage; exit;;
  esac
done

# FQ1=${@:$OPTIND:1} #read1 fastq file (demultipled)
# FQ2=${@:$OPTIND:2} #read2 fastq file (demultiplexed)

if [[ -z "$SAMPLE" ]]; then usage; echo; echo "Error: sample name is required."; exit 1; fi
if [[ -z "$sample_name_UDI" ]]; then usage; echo; echo "Error: sample name with UDI is required"; exit 1; fi
if [[ -z "$SISYPHUS_OUTFOLDER" ]]; then usage; echo; echo "Error: Sisyphus folder to transfer data to must be supplied"; exit 1; fi

################################################################################
# Environment setup
################################################################################
P=`dirname $0`
source $P/Watson_code_environment_setup_v2.6.sh
initialize

################################################################################
# Paths that are created/used in the pipeline
################################################################################
#FILES PRODUCED/ USED IN THE PROCESS:
SSCS_pindel_config_file="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_config_file_SSCS.txt"
DCS_pindel_config_file="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_config_file_DCS.txt"

SSCS_bam="$sample_name_UDI/FLT3_calling/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_mapped_merged.bam"
DCS_bam="$sample_name_UDI/FLT3_calling/DCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_mapped_merged.bam"

pindel_SSCS_output="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_SSCS"
pindel_DCS_output="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_DCS"

pindel_tandem_duplications_SSCS="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_SSCS_TD"
pindel_tandem_duplications_DCS="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_DCS_TD"

pindel_vcf_SSCS="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_SSCS"
pindel_vcf_DCS="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_pindel_DCS"

FLT3_ITD_calls_SSCS="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_FLT3_ITD_calls_SSCS.csv"
FLT3_ITD_calls_DCS="$sample_name_UDI/FLT3_calling/Pindel/${SAMPLE}_FLT3_ITD_calls_DCS.csv"

################################################################################
# Run the pipeline
################################################################################

execute "mkdir -p TEMP"
execute "mkdir -p $sample_name_UDI/FLT3_calling/Pindel"

################################################################################
banner "Create text configuration files for Pindel"

execute "python $pindel_config_file --infile $SSCS_bam --outfile $SSCS_pindel_config_file --insert_size 300 --sample_name ${SAMPLE}_SSCS"
execute "python $pindel_config_file --infile $DCS_bam --outfile $DCS_pindel_config_file --insert_size 300 --sample_name ${SAMPLE}_DCS"

################################################################################
## SSCS ######
################################################################################

banner "Run Pindel to look for FLT3 ITD..."

execute "$pindel -f $REF -i $SSCS_pindel_config_file -c 13 -o $pindel_SSCS_output -T $THREADS -e 0.000001"

banner "Convert Pindel output file to VCF..."

execute "$pindel2vcf -r $REF -R hg19 -d 00000000 -p $pindel_tandem_duplications_SSCS -v $pindel_vcf_SSCS -c 13"

banner "Convert Pindel VCF output to text file..."

execute "python $watson_call_FLT3_ITD --pindel_VCF $pindel_vcf_SSCS --outfile $FLT3_ITD_calls_SSCS --sample_name ${SAMPLE}_SSCS --min_depth $MINDEPTH --DCS_or_SSCS SSCS"

################################################################################
## SSCS ######
################################################################################

banner "Run Pindel to look for FLT3 ITD..."

execute "$pindel -f $REF -i $DCS_pindel_config_file -c 13 -o $pindel_DCS_output -T $THREADS -e 0.000001"

banner "Convert Pindel output file to VCF..."

execute "$pindel2vcf -r $REF -R hg19 -d 00000000 -p $pindel_tandem_duplications_DCS -v $pindel_vcf_DCS -c 13"

banner "Convert Pindel VCF output to text file..."

execute "python $watson_call_FLT3_ITD --pindel_VCF $pindel_vcf_DCS --outfile $FLT3_ITD_calls_DCS --sample_name ${SAMPLE}_DCS --min_depth $MINDEPTH --DCS_or_SSCS DCS"


# # ################################################################################

execute "echo -ne '\007'" #make a beep noise

banner "Completed."

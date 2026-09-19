#!/bin/bash
##################################################################################
# Prepares the consensus BAMs used by the FLT3-ITD (Pindel) caller.
#
# The FLT3-ITD arm cannot re-use the SNV panel's SSCS/DCS BAMs. Those are called
# with Watson_code_SSCS_calling_2.1.py, which filters reads on their SAM flags;
# the reads that carry an ITD map discordantly or with long soft-clips, so that
# filtering removes the evidence. SSCS and DCS are therefore called a second
# time from the same mapped-merged BAM, using the consensus callers that do no
# SAM-flag filtering, and written where Watson_code_Pindel_FLT3_caller.sh looks:
#
#   <sample_name_UDI>/FLT3_calling/SSCS_files_MUFS<m>/<SAMPLE>_SNV_watson_code_SSCS_mapped_merged.bam
#   <sample_name_UDI>/FLT3_calling/DCS_files_MUFS<m>/<SAMPLE>_SNV_watson_code_DCS_mapped_merged.bam
#
# Order: Watson_code_SNV_panel_v1.7.sh -> this script -> Watson_code_Pindel_FLT3_caller.sh
#
# No ClipBam and no GATK indel realignment here, unlike the SNV panel: hard-clipping
# the read ends and the overlapping portion of a pair removes bases Pindel needs to
# anchor the duplicated segment, and Pindel does its own split-read analysis. This
# matches what Watson_code_Pindel_FLT3_caller.sh expects -- its inputs are the
# *_mapped_merged.bam files, i.e. the MergeBamAlignment output, upstream of both.
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
CALLMINDUPLEX=3
QUALMIN=20
BASEQUALMIN=20
PROPMIN=0.9
MAXN=1.0
Duplex_MAXN=1.0
OUTFOLDER=''
INPUT_BAM=''

usage() {
  cat << EOF
Usage: $0 [options] -s <SAMPLE> -u <sample_name_UDI>

Calls SSCS and DCS consensus reads for FLT3-ITD detection, from the mapped
merged BAM produced by Watson_code_SNV_panel_v1.7.sh. Run this before
Watson_code_Pindel_FLT3_caller.sh.

OPTIONS:
   -s <SAMPLE> Sample name
   -u <sample_name_UDI> Sample name with UDI
   -B <bam> mapped merged BAM to start from
            (default: <sample_name_UDI>/<SAMPLE>_SNV_mapped_merged_bam.bam)
   -m <SSCS-min-reads> The minimum number of reads to form a SSCS read (default: $CALLMIN)
   -a <SSCS-min-reads-for-DCS-calling> The minimum SSCS family size to be included in DCS read calling (default: $CALLMINDUPLEX)
   -q <min-mapping-quality> The minimum mapping quality to include (default: $QUALMIN)
   -b <min-base-quality> The minimum base quality to include in consensus calling (default: $BASEQUALMIN)
   -p <SSCS-threshold> The minimum proportion of nucleotides at a position in a read that must be the same in order for a SSCS to be called at that position (default: $PROPMIN)
   -n <SSCS-max-N> Maximum fraction of Ns permitted in a SSCS consensus (default: $MAXN)
   -d <DCS-max-N> Maximum fraction of Ns permitted in a DCS consensus (default: $Duplex_MAXN)
   -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
   -o <output-folder> library folder name (accepted for consistency with the other panel scripts)
EOF
}

while getopts "s:u:B:m:a:q:b:p:n:d:t:o:" OPTION; do
  case $OPTION in
    s) SAMPLE=$OPTARG;;
    u) sample_name_UDI=$OPTARG;;
    B) INPUT_BAM=$OPTARG;;
    m) CALLMIN=$OPTARG;;
    a) CALLMINDUPLEX=$OPTARG;;
    q) QUALMIN=$OPTARG;;
    b) BASEQUALMIN=$OPTARG;;
    p) PROPMIN=$OPTARG;;
    n) MAXN=$OPTARG;;
    d) Duplex_MAXN=$OPTARG;;
    t) THREADS=$OPTARG;;
    o) OUTFOLDER=$OPTARG;;
    h) echo "Unknown option ${OPTION}"; usage; exit;;
    [?]) usage; exit;;
  esac
done

if [[ -z "$SAMPLE" ]]; then usage; echo; echo "Error: sample name is required."; exit 1; fi
if [[ -z "$sample_name_UDI" ]]; then usage; echo; echo "Error: sample name with UDI is required"; exit 1; fi

################################################################################
# Environment setup
################################################################################
P="$(dirname "$0")"
source "$P/Watson_code_environment_setup_v2.6.sh"
initialize

################################################################################
# Paths that are created/used in the pipeline
################################################################################
# Input: the mapped merged BAM from the SNV panel (raw reads, pre-consensus)
if [[ -z "$INPUT_BAM" ]]; then
  INPUT_BAM="$sample_name_UDI/${SAMPLE}_SNV_mapped_merged_bam.bam"
fi

if [[ ! -f "$INPUT_BAM" ]]; then
  echo "Error: mapped merged BAM not found: $INPUT_BAM"
  echo "Run Watson_code_SNV_panel_v1.7.sh for this sample first, or pass the BAM with -B."
  exit 1
fi

# The Pindel caller builds these directory names from its own -m value only, so
# both consensus levels live in MUFS${CALLMIN} folders. Keep the two scripts'
# -m arguments identical.
FLT3_directory="$sample_name_UDI/FLT3_calling"
SSCS_directory="$FLT3_directory/SSCS_files_MUFS${CALLMIN}"
DCS_directory="$FLT3_directory/DCS_files_MUFS${CALLMIN}"

#FILES PRODUCED IN THE PROCESS:
unpaired_SSCS_bam="$SSCS_directory/${SAMPLE}_SNV_watson_code_unpaired_SSCS_bam.bam"
SSCS_bam="$SSCS_directory/${SAMPLE}_SNV_watson_code_SSCS_bam.bam"
SSCS_bam_unmapped="$SSCS_directory/${SAMPLE}_SNV_watson_code_SSCS_unmapped.bam"
SSCS_bam_unmapped_sorted="$SSCS_directory/${SAMPLE}_SNV_watson_code_SSCS_unmapped_sorted.bam"
SSCS_adapter_marked_bam="$SSCS_directory/${SAMPLE}_SNV_watson_code_SSCS_adapter_marked.bam"
SSCS_mapped_merged="$SSCS_directory/${SAMPLE}_SNV_watson_code_SSCS_mapped_merged.bam"

DCS_bam="$DCS_directory/${SAMPLE}_SNV_watson_code_DCS_bam.bam"
DCS_bam_unmapped="$DCS_directory/${SAMPLE}_SNV_watson_code_DCS_unmapped.bam"
DCS_bam_unmapped_sorted="$DCS_directory/${SAMPLE}_SNV_watson_code_DCS_unmapped_sorted.bam"
DCS_adapter_marked_bam="$DCS_directory/${SAMPLE}_SNV_watson_code_DCS_adapter_marked.bam"
DCS_mapped_merged="$DCS_directory/${SAMPLE}_SNV_watson_code_DCS_mapped_merged.bam"

#METRICS FILES:
SSCS_adapter_metrics="$FLT3_directory/Metrics_and_images/${SAMPLE}_SNV_FLT3_SSCS_MarkIlluminaAdapter_metrics.txt"
DCS_adapter_metrics="$FLT3_directory/Metrics_and_images/${SAMPLE}_SNV_FLT3_DCS_MarkIlluminaAdapter_metrics.txt"

################################################################################
# Run the pipeline
################################################################################

execute "mkdir -p TEMP"
execute "mkdir -p $FLT3_directory"
execute "mkdir -p $FLT3_directory/Metrics_and_images"
execute "mkdir -p $SSCS_directory"
execute "mkdir -p $DCS_directory"

############### SSCS CALLING (no SAM flag filtering) ###########################
################################################################################
banner "Grouping reads by UMI and calling SSCS for FLT3 calling (watson code)..."

execute "python $watson_call_SSCS_translocations_or_FLT3 --infile $INPUT_BAM --sample-name $SAMPLE --min-family-size $CALLMIN" \
        " --threshold $PROPMIN --max_N $MAXN --outbam $SSCS_bam --min-mapping-quality $QUALMIN"  \
        " --min-base-quality $BASEQUALMIN --unpaired-outbam $unpaired_SSCS_bam --out-directory $FLT3_directory"

execute "rm TEMP/${SAMPLE}_grouped_reads_dict*"

################################################################################

banner "Create unmapped version of SSCS bam file..."

execute "python $watson_code_unmap_bam --infile $SSCS_bam --outbam $SSCS_bam_unmapped --sample-name $SAMPLE"

banner "Sort unmapped SSCS bam file..."

execute "java -jar $picard SortSam I=$SSCS_bam_unmapped O=$SSCS_bam_unmapped_sorted SORT_ORDER=queryname"

banner "Mark Illumina adapters in unmapped SSCS bam file..."

execute "java -Xmx32G -jar $picard MarkIlluminaAdapters I=$SSCS_bam_unmapped_sorted" \
        " O=$SSCS_adapter_marked_bam METRICS=$SSCS_adapter_metrics TMP_DIR=TEMP"

banner "SamToFastq, BWA and MergeBamAlignment..."

execute "java -Xmx32G -jar $picard SamToFastq I=$SSCS_adapter_marked_bam F=/dev/stdout INTERLEAVE=true" \
        " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X TMP_DIR=TEMP | bwa mem -p -t $THREADS $REF /dev/stdin |" \
        " java -Xmx32G -jar $picard MergeBamAlignment UNMAPPED=$SSCS_adapter_marked_bam ALIGNED=/dev/stdin" \
        " O=$SSCS_mapped_merged R=$REF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1" \
        " ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=TEMP "

# Pindel expects a <file>.bam.bai index; picard writes <file>.bai
execute "samtools index $SSCS_mapped_merged"

############### DCS CALLING ####################################################
################################################################################
banner "Calling DCS for FLT3 calling (watson code)..."

execute "python $watson_call_DCS_FLT3 --infile $SSCS_bam --sample-name $SAMPLE --min-family-size-SSCS $CALLMINDUPLEX" \
        " --max_N $Duplex_MAXN --outbam $DCS_bam --min-base-quality $BASEQUALMIN"  \
        " --out-directory $FLT3_directory"

################################################################################

banner "Create unmapped version of DCS bam file..."

execute "python $watson_code_unmap_bam --infile $DCS_bam --outbam $DCS_bam_unmapped --sample-name $SAMPLE"

banner "Sort unmapped DCS bam file..."

execute "java -jar $picard SortSam I=$DCS_bam_unmapped O=$DCS_bam_unmapped_sorted SORT_ORDER=queryname"

banner "Mark Illumina adapters in unmapped DCS bam file..."

execute "java -Xmx32G -jar $picard MarkIlluminaAdapters I=$DCS_bam_unmapped_sorted" \
        " O=$DCS_adapter_marked_bam METRICS=$DCS_adapter_metrics TMP_DIR=TEMP"

banner "SamToFastq, BWA and MergeBamAlignment..."

execute "java -Xmx32G -jar $picard SamToFastq I=$DCS_adapter_marked_bam F=/dev/stdout INTERLEAVE=true" \
        " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X TMP_DIR=TEMP | bwa mem -p -t $THREADS $REF /dev/stdin |" \
        " java -Xmx32G -jar $picard MergeBamAlignment UNMAPPED=$DCS_adapter_marked_bam ALIGNED=/dev/stdin" \
        " O=$DCS_mapped_merged R=$REF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1" \
        " ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=TEMP "

execute "samtools index $DCS_mapped_merged"

################################################################################

execute "rm $SSCS_bam_unmapped" #have the sorted version so don't need this as well
execute "rm $DCS_bam_unmapped" #have the sorted version so don't need this as well

################################################################################

banner "Completed. Now run Watson_code_Pindel_FLT3_caller.sh for this sample."

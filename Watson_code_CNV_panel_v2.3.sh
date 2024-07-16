#!/bin/bash
##################################################################################

set -euo pipefail #this command makes sure the script will be stopped when the first error is encountered.

################################################################################
# Argument processing and setup
################################################################################

#DEFAULTS
CSV=""
FQ1=""
FQ2=""
THREADS=4
CALLMIN=1
QUALMIN=20
BASEQUALMIN=20
PROPMIN=0.9
MAXN=1.0
CLIPBASES=3

usage() {
  cat << EOF
Usage: $0 [options] -s <sample_indexes.csv> -o <library folder name> -f <read1.fq> -g <read2.fq>

Analyzes paired end fastq files (fastq read1 file and fastq read2 file).
OPTIONS:
   -s <sample_indexes.csv> The csv file containing the link between the sample name and the UDI index name
   -m <SSCS-min-reads> The minimum number of reads to form a SSCS read (default: $CALLMIN)
   -q <min-mapping-quality> The minimum mapping quality to include (default: $QUALMIN)
   -b <min-base-quality> The minimum base quality to include in consensus calling (default: $BASEQUALMIN)
   -p <SSCS-threshold> The minimum proportion of nucleotides at a position in a read that must be the same in order for a SSCS to be called at that position [default = $PROPMIN]
   -n <SSCS-max-N> Maximum fraction of Ns permitted in a SSCS consensus [default = $MAXN]
   -c <clip-bases> Number of bases to clip from the 5' and 3' end of each read [default = $CLIPBASES]
   -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
   -o <library folder name> The name of the library folder on Sisyphus to copy the files too (e.g. 'SLX_20125')
   -f <read1.fq> read1 fastq file
   -g <read2.fq> read2 fastq file
EOF
}

while getopts “f:g:s:m:p:q:b:n:c:t:o:” OPTION; do
  case $OPTION in
    f) FQ1=$OPTARG;;
    g) FQ2=$OPTARG;;
    s) CSV=$OPTARG;;
    m) CALLMIN=$OPTARG;;
    p) PROPMIN=$OPTARG;;
    q) QUALMIN=$OPTARG;;
    b) BASEQUALMIN=$OPTARG;;
    n) MAXN=$OPTARG;;
    c) CLIPBASES=$OPTARG;;
    t) THREADS=$OPTARG;;
    o) SISYPHUS_OUTFOLDER=$OPTARG;;
    h) echo "Unknown option ${OPTION}"; usage; exit;;
    [?]) usage; exit;;
  esac
done

# FQ1=${@:$OPTIND:1} #read1 fastq file (demultipled)
# FQ2=${@:$OPTIND:2} #read2 fastq file (demultiplexed)

if [[ -z "$CSV" ]]; then usage; echo; echo "Error: csv file containing link between sample name and UDI index is required."; exit 1; fi
if [[ -z "$FQ1" ]]; then usage; echo; echo "Error: fastq read1 file must be supplied"; exit 1; fi
if [[ -z "$FQ2" ]]; then usage; echo; echo "Error: fastq read2 file must be supplied"; exit 1; fi
if [[ -z "$SISYPHUS_OUTFOLDER" ]]; then usage; echo; echo "Error: Sisyphus folder to transfer data to must be supplied"; exit 1; fi

################################################################################
# Environment setup
################################################################################
P=`dirname $0`
source $P/Watson_code_environment_setup_v2.6.sh
initialize

################################################################################
# RETRIEVE SAMPLE NAME FOR FILE NAMING
################################################################################
banner "Extract sample name from UDI index name"

SAMPLE=$(python $watson_code_sample --samples_csv_file $CSV --fastq $FQ1 2>&1)
sample_name_UDI=$(python $watson_code_sample_UDI --samples_csv_file $CSV --fastq $FQ1 2>&1)
library_name=$(python $watson_code_library_name --fastq $FQ1 2>&1)

################################################################################
# Paths that are created/used in the pipeline
################################################################################
#FILES PRODUCED IN THE PROCESS:
unmapped_bam="$sample_name_UDI/${SAMPLE}_unmapped_bam.bam"
unmapped_bam_with_UMI="$sample_name_UDI/${SAMPLE}_unmapped_bam_with_RX-UMI.bam"
adapter_marked_bam="$sample_name_UDI/${SAMPLE}_adapter_marked_bam.bam"
mapped_merged_BAM="$sample_name_UDI/${SAMPLE}_mapped_merged_bam.bam"
unpaired_SSCS_bam="$sample_name_UDI/${SAMPLE}_watson_code_unpaired_SSCS_bam.bam"
SSCS_bam="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_bam.bam"
SSCS_bam_unmapped="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_unmapped.bam"
SSCS_bam_unmapped_sorted="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_unmapped_sorted.bam"
SSCS_adapter_marked_bam="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_adapter_marked.bam"
SSCS_mapped_merged="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_mapped_merged.bam"
SSCS_clipped="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_mapped_clipped.bam"
GATK_intervals="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_mapped_clipped_GATK.intervals"
GATK_bam="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_mapped_clipped_GATK_realigned.bam"
annovar_annotated_SSCS="$sample_name_UDI/${SAMPLE}_watson_code_SSCS_variant_calling_annovar"

#METRICS FILES:
adapter_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_MarkIlluminaAdapter_metrics.txt"
insert_size_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_insert_size_metrics.txt"
insert_size_histogram="$sample_name_UDI/Metrics_and_images/${SAMPLE}_insert_size_histogram.pdf"
SSCS_overlap_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SSCS_clip_overlap_metrics.txt"
SSCS_adapter_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SSCS_MarkIlluminaAdapter_metrics.txt"

################################################################################
# Run the pipeline
################################################################################

execute "mkdir -p $sample_name_UDI"
execute "mkdir -p TEMP"
execute "mkdir -p $sample_name_UDI/Metrics_and_images"
execute "mkdir -p $sample_name_UDI/CNV_read_depths"
execute "mkdir -p $sample_name_UDI/SNPs_plots"
execute "mkdir -p $sample_name_UDI/KMT2A-PTD"
execute "mkdir -p $sample_name_UDI/Translocation_breakpoint_read_depths"

banner "Convert Fastq to Unmapped BAM..."

execute "java -Xmx32G -jar $picard FastqToSam F1=$FQ1 F2=$FQ2 O=$unmapped_bam SM=$SAMPLE LB=$library_name "\
        "PL=Illumina DS=$SAMPLE RG=$SAMPLE TMP_DIR=TEMP"

banner "Extract inline UMIs and store in RX tag of of BAM..."

execute "java -Xmx32G -jar $fgbio --tmp-dir ./TEMP ExtractUmisFromBam --input=$unmapped_bam --output=$unmapped_bam_with_UMI "\
        "--read-structure=3M2S146T 3M2S146T --molecular-index-tags=ZA ZB --single-tag=RX"

banner "MarkIlluminaAdapters..."

execute "java -Xmx32G -jar $picard MarkIlluminaAdapters I=$unmapped_bam_with_UMI O=$adapter_marked_bam METRICS=$adapter_metrics TMP_DIR=TEMP"

banner "SamToFastq, BWA and MergeBamAlignment..."

execute "java -Xmx32G -jar $picard SamToFastq I=$adapter_marked_bam F=/dev/stdout INTERLEAVE=true" \
        " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X TMP_DIR=TEMP | bwa mem -p -t 8 $REF /dev/stdin |" \
        " java -Xmx32G -jar $picard MergeBamAlignment UNMAPPED=$adapter_marked_bam ALIGNED=/dev/stdin" \
        " O=$mapped_merged_BAM R=$REF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1" \
        " ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=TEMP"

banner "CollectInsertSizeMetrics..."

execute "java -jar $picard CollectInsertSizeMetrics I=$mapped_merged_BAM O=$insert_size_metrics H=$insert_size_histogram"

################################################################################
banner "Grouping reads by UMI and calling SSCS (watson code)..."

execute "python $watson_call_SSCS --infile $mapped_merged_BAM --sample-name $SAMPLE --min-family-size $CALLMIN" \
        " --threshold $PROPMIN --max_N $MAXN --outbam $SSCS_bam --min-mapping-quality $QUALMIN"  \
        " --min-base-quality $BASEQUALMIN --unpaired-outbam $unpaired_SSCS_bam --out-directory $sample_name_UDI"

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
         " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X TMP_DIR=TEMP | bwa mem -p -t 8 $REF /dev/stdin |" \
         " java -Xmx32G -jar $picard MergeBamAlignment UNMAPPED=$SSCS_adapter_marked_bam ALIGNED=/dev/stdin" \
         " O=$SSCS_mapped_merged R=$REF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1" \
         " ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=TEMP "

banner "Clip Overlapping Reads and 3 bases from end of each read (ClipBam)"

execute "java -Xmx32G -jar $fgbio --tmp-dir TEMP ClipBam --input $SSCS_mapped_merged"\
        " --output $SSCS_clipped --ref $REF --clipping-mode Hard --clip-overlapping-reads true" \
        " --read-one-five-prime $CLIPBASES --read-one-three-prime $CLIPBASES --read-two-five-prime $CLIPBASES"\
        " --read-two-three-prime $CLIPBASES --metrics $SSCS_overlap_metrics"

################################################################################

banner "Perform realignment around indels"

execute "java -Xmx32G -jar $GATK -T RealignerTargetCreator -R $REF"\
        " -I $SSCS_clipped -o $GATK_intervals"

execute "java -Xmx32G -jar $GATK -T IndelRealigner -R $REF"\
        " -I $SSCS_clipped -targetIntervals $GATK_intervals -o $GATK_bam"

################################################################################
banner "Calculating and Plotting Read Depths"

execute "python $watson_SSCS_depths --infile $GATK_bam --sample-name $SAMPLE --panel_bed $CNV_bed"\
        " --chromosome-ideogram $chromosome_ideogram --out-directory $sample_name_UDI"

################################################################################
banner "Calculating and Plotting KMT2A Read Depths"

execute "python $watson_KMT2A_PTD --infile $GATK_bam --sample-name $SAMPLE --panel-bed $CNV_bed"\
        " --KMT2A-coordinates $KMT2A_coordinates --out-directory $sample_name_UDI"

################################################################################
banner "Plotting translocation breakpoint read coverages"

execute "python $watson_translocation_breakpoint_depths --infile $GATK_bam --sample-name $SAMPLE"\
        " --panel-bed $CNV_bed --out-directory $sample_name_UDI"

###############################################################################
banner "Create samtools mpileup and custom watson code VCF file and txt file"

# execute "export PATH=$HOME/Pipeline_tools/samtools-1.11:$PATH"

execute "samtools mpileup -BOa -d1000000 -f$REF -l$CNV_bed_GATK -Q0"\
        " --output-MQ --output-QNAME --reverse-del $GATK_bam | $watson_VCF_text --infile /dev/stdin"\
        " --sample-name $SAMPLE --path-to-reference-genome $REF --ECS_type SSCS"\
        " --dbSNP_directory $dbSNP_directory --mean_base_qual_filter 22.5 --mapq_filter 0"\
        " --read_pos_filter 8 --bias_filter 100 --min_family_size $CALLMIN --out-directory $sample_name_UDI"

execute "rm $sample_name_UDI/dbSNP_dictionary*"

###############################################################################
banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR -out $annovar_annotated_SSCS -build hg19 -hgvs $sample_name_UDI/${SAMPLE}_watson_code_SSCS_variant_calling_variants_and_SNPs_MUFS${CALLMIN}.txt"\
        " $ANNOVAR_humandb"

###############################################################################
banner "Plot SNP BAFs"

execute "python $watson_SNP_plots --targeted_SNPs $targeted_SNPs --annovar_exonic ${annovar_annotated_SSCS}.exonic_variant_function"\
        " --annovar_variant ${annovar_annotated_SSCS}.variant_function --annovar_invalid ${annovar_annotated_SSCS}.invalid_input"\
        " --sample-name $SAMPLE --panel_bed $CNV_bed --chromosome-ideogram $chromosome_ideogram --out-directory $sample_name_UDI"

###############################################################################
banner "Completed."

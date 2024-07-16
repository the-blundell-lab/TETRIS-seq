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
CALLMIN=3
CALLMINDUPLEX=3
QUALMIN=20
BASEQUALMIN=20
PROPMIN=0.9
MAXN=1.0
Duplex_MAXN=1.0
CLIPBASES=3

usage() {
  cat << EOF
Usage: $0 [options] -s <sample_indexes.csv> -o <library folder name> -f <read1.fq> -g <read2.fq>

Analyzes paired end fastq files (fastq read1 file and fastq read2 file).
OPTIONS:
    -s <sample_indexes.csv> The csv file containing the link between the sample name and the UDI index name
    -m <SSCS-min-reads> The minimum number of reads to form a SSCS read (default: $CALLMIN)
    -a <SSCS-min-reads-for-DCS-calling> The minimum SSCS family size to be included in DCS read calling (default: $CALLMINDUPLEX)
    -q <min-mapping-quality> The minimum mapping quality to include (default: $QUALMIN)
    -b <min-base-quality> The minimum base quality to include in consensus calling (default: $BASEQUALMIN)
    -p <SSCS-threshold> The minimum proportion of nucleotides at a position in a read that must be the same in order for a SSCS to be called at that position [default = $PROPMIN]
    -n <SSCS-max-N> Maximum fraction of Ns permitted in a SSCS consensus [default = $MAXN]
    -d <DCS-max-N> Maximum fraction of Ns permitted in a DCS consensus [default = $Duplex_MAXN]
    -c <clip-bases> Number of bases to clip from the 5' and 3' end of each read [default = $CLIPBASES]
    -t <cpus> The maximum number of threads/CPUs to use during analysis (default: $THREADS)
    -o <library folder name> The name of the library folder to copy the files too (e.g. 'SLX_20125')
    -f <read1.fq> read1 fastq file
    -g <read2.fq> read2 fastq file
EOF
}

while getopts “f:g:s:m:a:q:b:p:n:d:c:t:o:” OPTION; do
  case $OPTION in
    f) FQ1=$OPTARG;;
    g) FQ2=$OPTARG;;
    s) CSV=$OPTARG;;
    m) CALLMIN=$OPTARG;;
    a) CALLMINDUPLEX=$OPTARG;;
    q) QUALMIN=$OPTARG;;
    b) BASEQUALMIN=$OPTARG;;
    p) PROPMIN=$OPTARG;;
    n) MAXN=$OPTARG;;
    d) Duplex_MAXN=$OPTARG;;
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
adapter_marked_bam="$sample_name_UDI/${SAMPLE}_SNV_adapter_marked_bam.bam"
mapped_merged_BAM="$sample_name_UDI/${SAMPLE}_SNV_mapped_merged_bam.bam"

unpaired_SSCS_bam="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_unpaired_SSCS_bam.bam"
SSCS_bam="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_bam.bam"
SSCS_bam_unmapped="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_unmapped.bam"
SSCS_bam_unmapped_sorted="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_unmapped_sorted.bam"
SSCS_adapter_marked_bam="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_adapter_marked.bam"
SSCS_mapped_merged="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_mapped_merged.bam"
SSCS_clipped="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_mapped_clipped.bam"
SSCS_GATK_intervals="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_mapped_clipped_GATK.intervals"
SSCS_GATK_bam="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_mapped_clipped_GATK_realigned.bam"
SSCS_PILEUP="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_samtools_mpileup.txt"
SSCS_variants_txt="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_variants_MUFs_${CALLMIN}.txt"
SSCS_variants_VCF="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_variants_MUFs_${CALLMIN}_all_positions.vcf"
annovar_annotated_SSCS="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_variants_MUFs_${CALLMIN}_annovar"
annotated_SSCS="$sample_name_UDI/SSCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_variants_MUFs_${CALLMIN}_annotated.txt"
VARDICT_SSCS_VCF="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava.vcf"
VARDICT_SSCS_TXT="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava.txt"
annovar_annotated_VarDict_SSCS="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava_annovar"
annotated_VarDict_SSCS="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava_annotated.txt"

VARDICT_SSCS_VCF_UNALIGNED="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava_unaligned.vcf"
VARDICT_SSCS_TXT_UNALIGNED="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava_unaligned.txt"
annovar_annotated_VarDict_SSCS_unaligned="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava_unaligned_annovar"
annotated_VarDict_SSCS_unaligned="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_SSCS_VarDictJava_unaligned_annotated.txt"

DCS_bam="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_bam.bam"
DCS_bam_unmapped="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_unmapped.bam"
DCS_bam_unmapped_sorted="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_unmapped_sorted.bam"
DCS_adapter_marked_bam="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_adapter_marked.bam"
DCS_mapped_merged="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_mapped_merged.bam"
DCS_clipped="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_mapped_clipped.bam"
DCS_GATK_intervals="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_mapped_clipped_GATK.intervals"
DCS_GATK_bam="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_mapped_clipped_GATK_realigned.bam"
DCS_PILEUP="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_samtools_mpileup.txt"
DCS_variants_txt="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_variants_MUFs_${CALLMINDUPLEX}.txt"
DCS_variants_VCF="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_variants_MUFs_${CALLMINDUPLEX}_all_positions.vcf"
annovar_annotated_DCS="$sample_name_UDI/DCS_files_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_variants_MUFs_${CALLMIN}_annovar"
annotated_DCS="$sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}/${SAMPLE}_SNV_watson_code_DCS_variants_MUFs_${CALLMINDUPLEX}_annotated.txt"
VARDICT_DCS_VCF="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava.vcf"
VARDICT_DCS_TXT="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava.txt"
annovar_annotated_VarDict_DCS="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava_annovar"
annotated_VarDict_DCS="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava_annotated.txt"

VARDICT_DCS_VCF_UNALIGNED="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava_unaligned.vcf"
VARDICT_DCS_TXT_UNALIGNED="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava_unaligned.txt"
annovar_annotated_VarDict_DCS_unaligned="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava_unaligned_annovar"
annotated_VarDict_DCS_unaligned="$sample_name_UDI/VarDictJava_MUFS${CALLMIN}/${SAMPLE}_SNV_watson_code_DCS_VarDictJava_unaligned_annotated.txt"

#METRICS FILES:
adapter_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_MarkIlluminaAdapter_metrics.txt"
insert_size_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_insert_size_metrics.txt"
insert_size_histogram="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_insert_size_histogram.pdf"
SSCS_adapter_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_SSCS_MarkIlluminaAdapter_metrics.txt"
DCS_adapter_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_DCS_MarkIlluminaAdapter_metrics.txt"
SSCS_overlap_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_SSCS_clip_overlap_metrics_MUFS${CALLMIN}.txt"
DCS_overlap_metrics="$sample_name_UDI/Metrics_and_images/${SAMPLE}_SNV_DCS_clip_overlap_metrics_MUFS${CALLMINDUPLEX}.txt"

################################################################################
# Run the pipeline
################################################################################

execute "mkdir -p $sample_name_UDI"
execute "mkdir -p TEMP"
execute "mkdir -p $sample_name_UDI/Metrics_and_images"
execute "mkdir -p $sample_name_UDI/Gene_read_depths"
execute "mkdir -p $sample_name_UDI/SSCS_files_MUFS${CALLMIN}"
execute "mkdir -p $sample_name_UDI/DCS_files_MUFS${CALLMINDUPLEX}"
execute "mkdir -p $sample_name_UDI/VarDictJava_MUFS${CALLMIN}"

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
        " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X TMP_DIR=./TEMP | bwa mem -p -t 8 $REF /dev/stdin |" \
        " java -Xmx32G -jar $picard MergeBamAlignment UNMAPPED=$adapter_marked_bam ALIGNED=/dev/stdin" \
        " O=$mapped_merged_BAM R=$REF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1" \
        " ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=./TEMP"

banner "CollectInsertSizeMetrics..."

execute "java -jar $picard CollectInsertSizeMetrics I=$mapped_merged_BAM O=$insert_size_metrics H=$insert_size_histogram"

################################################################################
banner "Plot read coverage for raw bam"

execute "python $watson_read_depths --infile $mapped_merged_BAM --sample-name $SAMPLE --panel_bed $SNV_bed "\
        "--BAM_type 'raw bam' --out-directory $sample_name_UDI --min-family-size-SSCS $CALLMIN"

############### SSCS CALLING ###############################################
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

execute "java -Xmx32G -jar $fgbio --tmp-dir ./TEMP ClipBam --input $SSCS_mapped_merged"\
       " --output $SSCS_clipped --ref $REF --clipping-mode Hard --clip-overlapping-reads true" \
       " --read-one-five-prime $CLIPBASES --read-one-three-prime $CLIPBASES --read-two-five-prime $CLIPBASES"\
       " --read-two-three-prime $CLIPBASES --metrics $SSCS_overlap_metrics"

################################################################################

banner "Perform realignment around indels"

execute "java -Xmx32G -jar $GATK -T RealignerTargetCreator -R $REF"\
       " -I $SSCS_clipped -o $SSCS_GATK_intervals"

execute "java -Xmx32G -jar $GATK -T IndelRealigner -R $REF"\
       " -I $SSCS_clipped -targetIntervals $SSCS_GATK_intervals -o $SSCS_GATK_bam"

################################################################################
banner "Plot read coverage for SSCS bam"

execute "python $watson_read_depths --infile $SSCS_GATK_bam --sample-name $SAMPLE --panel_bed $SNV_bed "\
       "--BAM_type 'SSCS' --out-directory $sample_name_UDI --min-family-size-SSCS $CALLMIN"

###############################################################################
banner "Create samtools mpileup and custom watson code VCF file and txt file"

execute "samtools mpileup -BOa -d1000000 -f$REF -l$SNV_bed_GATK -Q0"\
       " --output-MQ --output-QNAME --reverse-del $SSCS_GATK_bam | $watson_VCF_text_SNV_panel --infile /dev/stdin"\
       " --sample-name $SAMPLE --path-to-reference-genome $REF --ECS_type SSCS"\
       " --dbSNP_directory $dbSNP_directory --mean_base_qual_filter 22.5 --mapq_filter 0"\
       " --read_pos_filter 8 --bias_filter 100 --min_family_size $CALLMIN --out-directory $sample_name_UDI"\
       " --variants_txt $SSCS_variants_txt --variants_vcf $SSCS_variants_VCF"

execute "rm $sample_name_UDI/dbSNP_dictionary*"

###############################################################################
banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR_annotate $SSCS_variants_txt $ANNOVAR_humandb/ -buildver hg19"\
       " -out $annovar_annotated_SSCS"\
       " -remove -protocol refGene,cosmic92_coding,cosmic92_noncoding,exac03,gnomad_genome,clinvar_20200316 -operation g,f,f,f,f,f -nastring . -csvout"

banner "Create annotated SSCS txt file"

execute "python $annotate_SSCS --variants_infile $SSCS_variants_txt --annotated_infile ${annovar_annotated_SSCS}.hg19_multianno.csv"\
       " --min_family_size $CALLMIN --sample-name $SAMPLE --gene_transcripts $TWIST_transcript_details --outfile $annotated_SSCS"

############## SSCS Vardict Java CALLING ###############################################
###############################################################################

execute "vardict-java -G $REF -f 0.0000001 -b $SSCS_GATK_bam -z 0 -c 1 -S 2 -E 3 -g 4 -r 1 -q 0 -o 0 -u -y -th $THREADS $SNV_bed |"\
       " python $watson_VarDict_to_text --infile /dev/stdin --sample-name $SAMPLE --path-to-reference-genome $REF"\
       " --variants_txt $VARDICT_SSCS_TXT"

banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR_annotate $VARDICT_SSCS_TXT $ANNOVAR_humandb/ -buildver hg19"\
       " -out $annovar_annotated_VarDict_SSCS"\
       " -remove -protocol refGene,cosmic92_coding,cosmic92_noncoding,exac03,gnomad_genome,clinvar_20200316 -operation g,f,f,f,f,f -nastring . -csvout"

banner "Create annotated SSCS txt file"

execute "python $watson_VarDict_annotate --variants_infile $VARDICT_SSCS_TXT --annotated_infile ${annovar_annotated_VarDict_SSCS}.hg19_multianno.csv"\
       " --min_family_size $CALLMIN --sample-name $SAMPLE --gene_transcripts $TWIST_transcript_details --outfile $annotated_VarDict_SSCS"

execute "rm $VARDICT_SSCS_TXT"
execute "rm ${annovar_annotated_VarDict_SSCS}.hg19_multianno.csv"
execute "rm ${annovar_annotated_VarDict_SSCS}.refGene.invalid_input"
execute "rm ${annovar_annotated_VarDict_SSCS}.invalid_input"

############### SSCS Vardict Java CALLING - without VARDICT local indel realignment ###############################################
################################################################################

execute "vardict-java -G $REF -f 0.0000001 -b $SSCS_GATK_bam -z 0 -c 1 -S 2 -E 3 -g 4 -r 1 -q 0 -o 0 -k 0 -u -y -th $THREADS $SNV_bed |"\
       " python $watson_VarDict_to_text --infile /dev/stdin --sample-name $SAMPLE --path-to-reference-genome $REF"\
       " --variants_txt $VARDICT_SSCS_TXT_UNALIGNED"

banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR_annotate $VARDICT_SSCS_TXT_UNALIGNED $ANNOVAR_humandb/ -buildver hg19"\
       " -out $annovar_annotated_VarDict_SSCS_unaligned"\
       " -remove -protocol refGene,cosmic92_coding,cosmic92_noncoding,exac03,gnomad_genome,clinvar_20200316 -operation g,f,f,f,f,f -nastring . -csvout"

banner "Create annotated SSCS txt file"

execute "python $watson_VarDict_annotate --variants_infile $VARDICT_SSCS_TXT_UNALIGNED --annotated_infile ${annovar_annotated_VarDict_SSCS_unaligned}.hg19_multianno.csv"\
       " --min_family_size $CALLMIN --sample-name $SAMPLE --gene_transcripts $TWIST_transcript_details --outfile $annotated_VarDict_SSCS_unaligned"

execute "rm $VARDICT_SSCS_TXT_UNALIGNED"
execute "rm ${annovar_annotated_VarDict_SSCS_unaligned}.hg19_multianno.csv"
execute "rm ${annovar_annotated_VarDict_SSCS_unaligned}.refGene.invalid_input"
execute "rm ${annovar_annotated_VarDict_SSCS_unaligned}.invalid_input"


############## DCS CALLING ###############################################
################################################################################
banner "Calling DCS (watson code)..."

execute "python $watson_call_DCS --infile $SSCS_bam --sample-name $SAMPLE --min-family-size-SSCS $CALLMINDUPLEX" \
       " --max_N $Duplex_MAXN --outbam $DCS_bam --min-base-quality $BASEQUALMIN"  \
       " --out-directory $sample_name_UDI"

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
        " CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=X TMP_DIR=TEMP | bwa mem -p -t 8 $REF /dev/stdin |" \
        " java -Xmx32G -jar $picard MergeBamAlignment UNMAPPED=$DCS_adapter_marked_bam ALIGNED=/dev/stdin" \
        " O=$DCS_mapped_merged R=$REF SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1" \
        " ORIENTATIONS=FR VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true TMP_DIR=TEMP "

banner "Clip Overlapping Reads and 3 bases from end of each read (ClipBam)"

execute "java -Xmx32G -jar $fgbio --tmp-dir ./TEMP ClipBam --input $DCS_mapped_merged"\
       " --output $DCS_clipped --ref $REF --clipping-mode Hard --clip-overlapping-reads true" \
       " --read-one-five-prime $CLIPBASES --read-one-three-prime $CLIPBASES --read-two-five-prime $CLIPBASES"\
       " --read-two-three-prime $CLIPBASES --metrics $DCS_overlap_metrics"

################################################################################

banner "Perform realignment around indels"

execute "java -Xmx32G -jar $GATK -T RealignerTargetCreator -R $REF"\
       " -I $DCS_clipped -o $DCS_GATK_intervals"

execute "java -Xmx32G -jar $GATK -T IndelRealigner -R $REF"\
       " -I $DCS_clipped -targetIntervals $DCS_GATK_intervals -o $DCS_GATK_bam"

################################################################################
banner "Plot read coverage for DCS bam"

execute "python $watson_read_depths --infile $DCS_GATK_bam --sample-name $SAMPLE --panel_bed $SNV_bed "\
       "--BAM_type 'duplex' --out-directory $sample_name_UDI --min-family-size-SSCS $CALLMINDUPLEX"

###############################################################################
banner "Create samtools mpileup and custom watson code VCF file and txt file (from DCS)"

execute "samtools mpileup -BOa -d1000000 -f$REF -l$SNV_bed_GATK -Q0"\
       " --output-MQ --output-QNAME --reverse-del $DCS_GATK_bam | $watson_VCF_text_SNV_panel_duplex --infile /dev/stdin"\
       " --sample-name $SAMPLE --path-to-reference-genome $REF --ECS_type DCS"\
       " --dbSNP_directory $dbSNP_directory --mapq_filter 0"\
       " --read_pos_filter 8 --bias_filter 100 --min_family_size $CALLMINDUPLEX --out-directory $sample_name_UDI"\
       " --variants_txt $DCS_variants_txt --variants_vcf $DCS_variants_VCF"

execute "rm $sample_name_UDI/dbSNP_dictionary*"

###############################################################################
banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR_annotate $DCS_variants_txt $ANNOVAR_humandb/ -buildver hg19"\
       " -out $annovar_annotated_DCS"\
       " -remove -protocol refGene,cosmic92_coding,cosmic92_noncoding,exac03,gnomad_genome,clinvar_20200316 -operation g,f,f,f,f,f -nastring . -csvout"

banner "Create annotated DCS txt file"

execute "python $annotate_DCS --variants_infile $DCS_variants_txt --annotated_infile ${annovar_annotated_DCS}.hg19_multianno.csv"\
       " --min_family_size $CALLMIN --sample-name $SAMPLE --gene_transcripts $TWIST_transcript_details --outfile $annotated_DCS"

############### DCS Vardict Java CALLING ###############################################
################################################################################

execute "vardict-java -G $REF -f 0.0000001 -b $DCS_GATK_bam -z 0 -c 1 -S 2 -E 3 -g 4 -r 1 -q 0 -o 0 -u -y -th $THREADS $SNV_bed |"\
       " python $watson_VarDict_to_text --infile /dev/stdin --sample-name $SAMPLE --path-to-reference-genome $REF"\
       " --variants_txt $VARDICT_DCS_TXT"

banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR_annotate $VARDICT_DCS_TXT $ANNOVAR_humandb/ -buildver hg19"\
       " -out $annovar_annotated_VarDict_DCS"\
       " -remove -protocol refGene,cosmic92_coding,cosmic92_noncoding,exac03,gnomad_genome,clinvar_20200316 -operation g,f,f,f,f,f -nastring . -csvout"

banner "Create annotated DCS txt file"

execute "python $watson_VarDict_annotate --variants_infile $VARDICT_DCS_TXT --annotated_infile ${annovar_annotated_VarDict_DCS}.hg19_multianno.csv"\
       " --min_family_size $CALLMIN --sample-name $SAMPLE --gene_transcripts $TWIST_transcript_details --outfile $annotated_VarDict_DCS"

execute "rm $VARDICT_DCS_TXT"
execute "rm ${annovar_annotated_VarDict_DCS}.hg19_multianno.csv"
execute "rm ${annovar_annotated_VarDict_DCS}.refGene.invalid_input"
execute "rm ${annovar_annotated_VarDict_DCS}.invalid_input"

############### DCS Vardict Java CALLING - without Vardic local indel realignment ###############################################
################################################################################

execute "vardict-java -G $REF -f 0.0000001 -b $DCS_GATK_bam -z 0 -c 1 -S 2 -E 3 -g 4 -r 1 -q 0 -o 0 -k 0 -u -y -th $THREADS $SNV_bed |"\
       " python $watson_VarDict_to_text --infile /dev/stdin --sample-name $SAMPLE --path-to-reference-genome $REF"\
       " --variants_txt $VARDICT_DCS_TXT_UNALIGNED"

banner "Annotate with ANNOVAR"

execute "perl $ANNOVAR_annotate $VARDICT_DCS_TXT_UNALIGNED $ANNOVAR_humandb/ -buildver hg19"\
       " -out $annovar_annotated_VarDict_DCS_unaligned"\
       " -remove -protocol refGene,cosmic92_coding,cosmic92_noncoding,exac03,gnomad_genome,clinvar_20200316 -operation g,f,f,f,f,f -nastring . -csvout"

banner "Create annotated DCS txt file"

execute "python $watson_VarDict_annotate --variants_infile $VARDICT_DCS_TXT_UNALIGNED --annotated_infile ${annovar_annotated_VarDict_DCS_unaligned}.hg19_multianno.csv"\
       " --min_family_size $CALLMIN --sample-name $SAMPLE --gene_transcripts $TWIST_transcript_details --outfile $annotated_VarDict_DCS_unaligned"

execute "rm $VARDICT_DCS_TXT_UNALIGNED"
execute "rm ${annovar_annotated_VarDict_DCS_unaligned}.hg19_multianno.csv"
execute "rm ${annovar_annotated_VarDict_DCS_unaligned}.refGene.invalid_input"
execute "rm ${annovar_annotated_VarDict_DCS_unaligned}.invalid_input"
execute "rm $DCS_bam_unmapped" #have the sorted version so don't need this as well
execute "rm $SSCS_bam_unmapped" #have the sorted version so don't need this as well

###############################################################################
execute "echo -ne '\007'" #make a beep noise

banner "Completed."

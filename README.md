# TETRIS-seq
This repository contains the code for ‘TETRIS-seq’ as described/ used in the manuscript manuscript **"Evolutionary dynamics in the decade preceding acute myeloid leukaemia"**.  TETRIS-seq is a duplex error-corrected sequencing approach which incorporates 4 (‘Tetris’) key components: 1) AML/ clonal haematopoiesis SNV/indel detection, 2) genome-wide mCA detection, 3) AML-associated chromosomal rearrangement detection, 4) in silico noise correction for low frequency SNV calling.  

Code for the data analysis and Figure generation in the manuscript manuscript **"Evolutionary dynamics in the decade preceding acute myeloid leukaemia"** can be found here: https://github.com/the-blundell-lab/preAML_evolutionary_dynamics/blob/main/README.md

## TETRIS-seq computational workflow:
The TETRIS-seq computational workflow (contained in this repository) consists of several key steps:

**1) Processing of SNV/ indel panel fastq files:**
**2) Processing of mCA/ chromosomal rearrangement panel fastq files:**
**3) Processing of ‘mapped merged’ BAM files (produced in 1 and 2 above) for generating error-corrected BAM files suitable for chromosomal rearrangement calling or FLT3-ITD calling.**
**4) FLT3-ITD calling**
**5) Chromosomal rearrangement calling**
**6) _In silico_ noise correction model for low frequency SNV detection**


- Github file: _Watson_code_SNV_panel_v1.7.sh_
**2) **
   - Github file: _Watson_code_CNV_panel_v2.3.sh_
**3) **
   - Github file: _Watson_code_SSCS_calling_for_translocations_and_FLT3_1.1.py_
   - Github file: _Watson_code_DCS_calling_1.4_for_FLT3_calling.py_
**4) **
   - Github file: _Watson_code_Pindel_FLT3_caller.sh_
**6) **
     - Manta used for chromosomal rearrangement calling: https://github.com/Illumina/manta
**7) **
   - Github file: _In silico noise correction model/ Duplex_Error_Model_initial_variant_calling_v5.py_
   - Github file: _In silico noise correction model/ Duplex_Error_Model_post-model_variant processing.ipynb_
  
The TETRIS-seq workflow uses a number of software packages, including Picard, fgio, BWA, GATK, SAMtools, VarDictJava, Pindel, ANNOVAR, Manta, as well as custom Python scripts (included in repository). 


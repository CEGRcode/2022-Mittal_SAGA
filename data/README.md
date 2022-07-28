
This directory stores both data (FASTQ & BAM files) and generally used reference files.


### data/Rossi_2021_Supplementary_Data_1.txt
Supplementary Table 1 from Rossi et al, 2020.

### data/ChexMix_Peak_Filter_List_190612.bed
Blacklist regions used in Rossi et al, 2020.

### data/FASTQ/*.fastq.gz
Aligned samples (FASTQ format) from this publication. These files are the inputs for running quality checks using GenoPipe's EpitopeID module.

Populated by `/00_Download_and_preprocessing/1_download_files.sh`.
Contents used by
- `/00_Download_and_preprocessing/2_align_rna_samples.pbs`
- `/01_Run_GenoPipe/run_epitopeid.pbs`

### data/BAM/*.bam
Aligned samples (BAM format) from both the Yeast Epigenome Project (Rossi et al, 2020. __Nature__) and this publication. These files are used to generate Mittal gene class reference files and figures throughout the paper

Populated by `/00_Download_and_preprocessing/1_download_files.sh` and `/00_Download_and_preprocessing/2_align_rna_samples.pbs`.
Contents used throughout the repo:
- `/00_Download_and_preprocessing/3_normalized_samples.sh`
- `/01_Run_GenoPipe/run_deletionid.pbs`
- `/01_Run_GenoPipe/run_strainid.pbs`
- `/01_Run_GenoPipe/run_strainid.pbs`
- `/02_Create_Gene_Classes/1_tally_occupancies.sh`
- `/02_Create_Gene_Classes/2_build_Mittal_GeneClasses.sh`
- `/04_Bulk_Processing/job/bulk_cross.pbs`

### data/NormalizationFactors/*XXXX_ScalingFactors.out
Scaling factor for every BAM file in `data/BAM/*.bam`. Different scaling methodologies applied depending on the sample's enrichment target.

Populated by `/00_Download_and_preprocessing/3_normalize_samples.sh`.

### data/RefPT-YEP/*.bed
Reference points generated based on Rossi_2021_Supplementary_Data_1.txt. Script for generating these reference points can be found at `data/RefPT-YEP/build_refs.sh`.

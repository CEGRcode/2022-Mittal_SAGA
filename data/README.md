
This directory stores both data (FASTQ & BAM files) and generally used reference points.


### data/Rossi_2021_Supplementary_Data_1.txt
Supplementary Table 1 from Rossi et al, 2020.

### data/ChexMix_Peak_Filter_List_190612.bed
Blacklist regions used in Rossi et al, 2020.

### data/FASTQ/*.fastq.gz
Aligned samples (FASTQ format) from this publication. These files are the inputs for running quality checks using GenoPipe's EpitopeID module.

### data/BAM/*.bam
Aligned samples (BAM format) from both the Yeast Epigenome Project (Rossi et al, 2020. __Nature__) and this publication. These files are used to generate Mittal gene class reference files and figures throughout the paper

### data/RefPT-YEP/*.bed
Reference points generated based on Rossi_2021_Supplementary_Data_1.txt. Script for generating these reference points can be found at `data/RefPT-YEP/build_refs.sh`.

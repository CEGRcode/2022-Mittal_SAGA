# 2022-Mittal_SAGA


## data
This directory stores both data (FASTQ & BAM files) and generally used reference files.

## bin
This directory stores generally used scripts and executables.

## 00_Download_and_preprocessing
This directory includes the bash and PBS scripts for downloading and preprocessing the sequencing samples for analysis throughout the rest of this repo.

## 01_Run_GenoPipe
This directory includes the results from running each of the three GenoPipe modules (StrainID, DeletionID, EpitopeID) on all of our samples. Also included are the scripts for running GenoPipe and parsing the results into a tab-delimited format that can be viewed and easily read in Excel.

## 02_Create_Gene_Classes
This directory contains the scripts for making the Mittal et al, 2022 M-series and H-series reference points.

## 03_Figures
This directory contains the scripts for making the figure from the paper.

## 04_Bulk_Processing
This directory contains scripts to pileup every BAM file against every RefPT BED file.

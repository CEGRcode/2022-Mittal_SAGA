# An integrated SAGA and TFIID PIC assembly pathway selective for poised and induced promoters

## Chitvan Mittal<sup>1,2</sup>, Olivia Lang<sup>2</sup>, William K.M. Lai<sup>1,2</sup>, and B. Franklin Pugh<sup>1,2</sup>

1 Department of Biochemistry and Molecular Biology, Pennsylvania State University, University Park, Pennsylvania 16801, USA
2 Department of Molecular Biology and Genetics, Cornell University, Ithaca, New York 14850, USA

### Correspondence: <fp265@cornell.edu>


### PMID : [36302553](https://pubmed.ncbi.nlm.nih.gov/36302553/)
### GEO ID : [GSE212655](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212655)
* ChIP-exo : [GSE212654](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212654)
* RNA-seq : [GSE212653](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212653)

## Overview
Genome-wide, little is understood about how proteins organize at inducible promoters before and after induction and to what extent inducible and constitutive architectures depend on cofactors. We report that sequence-specific transcription factors and their tethered cofactors (e.g., SAGA [Spt–Ada–Gcn5–acetyltransferase], Mediator, TUP, NuA4, SWI/SNF, and RPD3-L) are generally bound to promoters prior to induction (“poised”), rather than recruited upon induction, whereas induction recruits the preinitiation complex (PIC) to DNA. Through depletion and/or deletion experiments, we show that SAGA does not function at constitutive promoters, although a SAGA-inde- pendent Gcn5 acetylates +1 nucleosomes there. When inducible promoters are poised, SAGA catalyzes +1 nucleo- some acetylation but not PIC assembly. When induced, SAGA catalyzes acetylation, deubiquitylation, and PIC assembly. Surprisingly, SAGA mediates induction by creating a PIC that allows TFIID (transcription factor II-D) to stably associate, rather than creating a completely TFIID-independent PIC, as generally thought. These findings suggest that inducible systems, where present, are integrated with constitutive systems.

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

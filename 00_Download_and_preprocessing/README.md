

This directory includes the bash and PBS scripts for downloading and preprocessing the sequencing samples for analysis throughout the rest of this repo.


The following scripts should be executed in numerical order. The user should update the working directory (`$WRK`) filepath variable within every script before executing.

### 1_download_files.sh
This bash script downloads many initial files used by scripts throughout the repository:
- the ScriptManager binary executable, saved to `2022-Mittal_SAGA/bin/`
- the Master NoTag control sample (Rossi, 2021), saved to `2022-Mittal_SAGA/data/BAM`
- the Mittal 2022 RNA sequencing samples (`pegr_ids_RNA.txt`), `*.fastq.gz` files saved to `2022-Mittal_SAGA/data/FASTQ` while `*.bam` files saved to `2022-Mittal_SAGA/data/BAM`
- the Mittal 2022 & Rossi 2021 ChIP-exo samples (`sample_ids.txt`), `*.fastq.gz` files saved to `2022-Mittal_SAGA/data/FASTQ` while `*.bam` files saved to `2022-Mittal_SAGA/data/BAM`


Update file path using your favorite text editor, then run:
```
bash 1_download_files.sh
```

### 2_align_rna_samples.pbs
This PBS submission script uses the RNA FASTQ files downloaded in `1_download_files.sh` to create BAM alignments using Hisat2. BAM alignments are sorted and indexed using samtools.

This submission script is compatible with systems running the PBS scheduler.

Update file path and allocation name using your favorite text editor, then run:
```
qsub 2_align_rna_samples.pbs
```

### 3_normalize_samples.sh
Calculate scaling factor for every sample (RNA and ChIP-exo) and save the files to `data/NormalizationFactors/XXXX_ScalingFactors.out`.
- ChIP-exo samples (Histone targets): use NFR occupancy (50bp around center) of 03-TFO and 04-UNB gene classes (defined in Rossi, 2021)
  - targets include "H2B", "H2BK123ub", "H3", "H3K14ac", and "H3K9ac"
- RNA samples: use ScriptManager's implementation of the NCIS (Liang & Keles 2012, *BMC Bioinformatics*))
- ChIP-exo samples (TF targets): use ScriptManager's implementation of the NCIS (Liang & Keles 2012, *BMC Bioinformatics*))

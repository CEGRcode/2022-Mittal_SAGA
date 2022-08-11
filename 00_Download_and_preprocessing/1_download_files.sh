#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=14gb
#PBS -l walltime=01:00:00
#PBS -A open
#PBS -o logs/download.log.out
#PBS -e logs/download.log.err

module load gcc
module load samtools
module load anaconda3
source activate mittal

# Update path and these variables for accessing datasets
WRK=/path/to/2022-Mittal_SAGA
SCRIPTMANAGER=$WRK/bin/ScriptManager-v0.13.jar
USER_EMAIL=
GALAXY_API_KEY=
PEGR_API_KEY=

cd $WRK
[ -d logs ] || mkdir logs

# Download ScriptManager
wget https://github.com/CEGRcode/scriptmanager/releases/download/v0.13/ScriptManager-v0.13.jar
mv ScriptManager-v0.13.jar $SCRIPTMANAGER

# Download Control Master NoTag sample
wget https://www.datacommons.psu.edu/download/eberly/pughlab/yeast-epigenome-project/masterNoTag_20180928.zip
unzip masterNoTag_20180928.zip
mv masterNoTag_20180928.bam data/BAM/
samtools index data/BAM/masterNoTag_20180928.bam

# Download all sequencing samples (Mittal & YEP)
python bin/download_ALL_PEGR_py3.py -u $USER_EMAIL -g $GALAXY_API_KEY -p $PEGR_API_KEY -n tmp-galaxy -f 00_Download_and_preprocessing/sample_ids.txt
mv tmp-galaxy/*.bam data/BAM/
mv tmp-galaxy/*.fastq.gz data/FASTQ/

# Download all RNA samples (Mittal)
python bin/download_ALL_PEGR_py3.py -u $USER_EMAIL -g $GALAXY_API_KEY -p $PEGR_API_KEY -n tmp-rna-galaxy -f 00_Download_and_preprocessing/pegr_ids_rna.txt
mv tmp-rna-galaxy/*.bam data/RNA-BAM/
mv tmp-rna-galaxy/*.fastq.gz data/RNA-FASTQ/

# Clean-up
rm masterNoTag_20180928.zip
rm -r tmp-galaxy
rm -r tmp-rna-galaxy

#!/bin/bash

# Get MD5 checksums for all raw FASTQ files

# Fill in placeholder constants with your directories
WRK=/path/to/2022-Mittal_SAGA/03_Bulk_Processing
FQDIR=$WRK/../data/FQ/
OUTDIR=$WRK/Checksums

cd $WRK
[ -d $OUTDIR ] || mkdir $OUTDIR

for PBS_ARRAYID in {1..554};
do
	# Determine FASTQ file for the current job array index
	FILE=`ls $FQDIR/*fastq.gz | head -n $PBS_ARRAYID |tail -1`
	BASE=`basename $FILE`

	md5sum $FILE > $OUTDIR/$BASE.mdsum
done

paste <(cat $OUTDIR/*_R1.fastq.gz.mdsum |sort -k2,2 ) <(cat $OUTDIR/*_R2.fastq.gz.mdsum |sort -k2,2 ) > raw_mdsums.txt
paste <(cat $OUTDIR/*forward.bw.mdsum |sort -k2,2 ) <(cat $OUTDIR/*reverse.bw.mdsum |sort -k2,2 ) > processed_mdsums.txt

#!/bin/bash
#PBS -l nodes=1:ppn=4
#PBS -l pmem=14gb
#PBS -l walltime=01:00:00
#PBS -A open
#PBS -o logs/make_occ_table.log.out
#PBS -e logs/make_occ_table.log.err

module load gcc
module load samtools
module load anaconda3
source activate mittal

# Compile occupancies from bulk_occupancy.pbs into a tab-delimited file

# Fill in placeholder constants with your directories
WRK=/path/to/2022_Mittal
OLIBRARY=$WRK/OCCUPANCY/#Name the reference to build the occupancy table for

cd $WRK

# Process BAM file for each BED file in directory
OREF=`basename $OLIBRARY`

echo "GeneName" > COLUMN

[ -f $OREF.tab  ] && rm $OREF.tab

for OFILE in `ls $OLIBRARY/NormCounts/*`;
do
	[ -f $OREF.tab ] || cut -f4 $OFILE >> COLUMN
	[ -f $OREF.tab ] || mv COLUMN $OREF.tab
	echo $OFILE > $OREF.tmp
	cut -f5 $OFILE >> $OREF.tmp
	paste $OREF.tab $OREF.tmp > $OREF.tmp2
	mv $OREF.tmp2 $OREF.tab
done

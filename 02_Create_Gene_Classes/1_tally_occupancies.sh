#!/bin/bash

# This script builds the inital reference files and calculates occupancies used to make the Mittal 2022 gene classes

# Normalization factors are calculated for each of the 6 samples (Mock Heat Shock and Heat Shock for each of Sua7, Spt7, and Taf2 targets)
# Occupancy calculation is made for each gene (centered around Sua7 or STM feature in promoter) and saved as a score in a BED file for each sample.
#   tmp/Sua7_sort-Sua7-12275-Offset-NormalizedCount.bed
#   tmp/Sua7_sort-Sua7-26344-Offset-NormalizedCount.bed
#   tmp/STM_sort-Spt7-11960-Offset-NormalizedCount.bed
#   tmp/STM_sort-Spt7-20115-Offset-NormalizedCount.bed
#   tmp/Sua7_sort-Taf2-11846-Offset-NormalizedCount.bed
#   tmp/Sua7_sort-Taf2-28736-Offset-NormalizedCount.bed

module load gcc
module load samtools
module load anaconda3
source activate mittal

WRK=/path/to/2022-Mittal_SAGA/02_Create_Gene_Classes
BIN=$WRK/../bin
DBAM=$WRK/../data/BAM/
FDIR=$WRK/NormalizationFactors/
CONTROL=$WRK/../data/BAM/masterNoTag_20180928.bam
BLACKLIST=$WRK/../data/ChexMix_Peak_Filter_List_190612.bed
REFPT=$WRK/../data/RefPT-YEP/
TEMP=$WRK/tmp

# Script shortcuts
SCRIPTMANAGER=$BIN/ScriptManager-v0.13.jar
SORT=$BIN/sort_BED_by_Score_v2.pl
SUM=$BIN/sum_Row_CDT.pl
UPDATES=$BIN/update_BED_score_with_TAB_score.pl

cd $WRK
[ -d $TEMP ] || mkdir $TEMP
[ -d $FDIR ] || mkdir $FDIR

# Download ScriptManager
if [ -f $SCRIPTMANAGER ]; then
	wget https://github.com/CEGRcode/scriptmanager/releases/download/v0.13/ScriptManager-v0.13.jar
	mv ScriptManager-v0.13.jar $SCRIPTMANAGER
fi

##--------Calculate Normalization Values--------
for BAM in "11846_Taf2_i5006_BY4741_-_XO_FilteredBAM" \
"12275_Sua7_i5006_BY4741_-_XO_FilteredBAM" \
"26344_Sua7_i5006_BY4741_-_XO_FilteredBAM" \
"11960_Spt7_i5006_BY4741_-_XO_FilteredBAM" \
"20115_Spt7_i5006_BY4741_-_XO_FilteredBAM" \
"28736_Taf2_i5006_BY4741_-_XO_FilteredBAM";
do
	BAMFILE=$DBAM/$BAM.bam
	[ -f $BAMFILE.bai ] || samtools index $BAMFILE
	java -jar $SCRIPTMANAGER read-analysis scaling-factor $BAMFILE --ncis -c $CONTROL -f $BLACKLIST -w 500 -o $FDIR/$BAM
done


##--------Pileup occupancy scores for Spt7, Sua7, and Taf2--------

count_occupancy () {
	REFFILE=$1
	BAMFILE=$2
	OCCUPANCY=$3

	BAM=`basename $BAMFILE ".bam"`
	REF=`basename $REFFILE ".bed"`
	BED=$REF\_200bp
	BASE=$TEMP/$BAM\_$BED\_read1

	cd $TEMP
	echo "EXPAND"
	# Expand BED 200bp from center
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 200 $REFFILE -o $BED.bed
	# Tag Pileup
	java -jar $SCRIPTMANAGER read-analysis tag-pileup $BED.bed $BAMFILE -M $BASE

	echo "MERGE"
	# Merge offset slices from anti and sense
	paste <(cut -f1,2,16-116 $BASE\_sense.cdt) <(cut -f47-147 $BASE\_anti.cdt) > CDT
	perl $SUM CDT TAB
	perl $UPDATES $REFFILE TAB BED
	perl $SORT BED desc $OCCUPANCY-TagCount.bed

	# Merge normalized slices from anti and sense
	FACTOR=`grep 'Scaling factor' $FDIR/$BAM\_ScalingFactors.out | awk -F" " '{print $3}'`
	java -jar $SCRIPTMANAGER read-analysis scale-matrix CDT -s $FACTOR -o NORMALIZED
	perl $SUM NORMALIZED NTAB
	perl $UPDATES $REFFILE NTAB NBED
	perl $SORT NBED desc $OCCUPANCY-NormalizedCount.bed

	# Clean-up
	#rm CDT NORMALIZED TAB BED
}

#===Sua7===
for SUA7 in "$DBAM/12275_Sua7_i5006_BY4741_-_XO_FilteredBAM.bam" \
 "$DBAM/26344_Sua7_i5006_BY4741_-_XO_FilteredBAM.bam";
do
	ID=`basename $SUA7 ".bam" | awk -F"_" '{print $1}'`
	count_occupancy $REFPT/Sua7.bed $SUA7 $TEMP/Sua7_sort-Sua7-$ID-Offset
done

#===Taf2===
for TAF2 in "$DBAM/11846_Taf2_i5006_BY4741_-_XO_FilteredBAM.bam" \
 "$DBAM/28736_Taf2_i5006_BY4741_-_XO_FilteredBAM.bam";
do
	ID=`basename $TAF2 ".bam" | awk -F"_" '{print $1}'`
	count_occupancy $REFPT/Sua7.bed $TAF2 $TEMP/Sua7_sort-Taf2-$ID-Offset
done

#===Spt7===
for SPT7 in "$DBAM/11960_Spt7_i5006_BY4741_-_XO_FilteredBAM.bam" \
 "$DBAM/20115_Spt7_i5006_BY4741_-_XO_FilteredBAM.bam";
do
	ID=`basename $SPT7 ".bam" | awk -F"_" '{print $1}'`
	count_occupancy $REFPT/STM.bed $SPT7 $TEMP/STM_sort-Spt7-$ID-Offset
	# count_occupancy $REFPT/UASv2.bed $SPT7 $TEMP/UASv2_sort-Spt7-$ID-Offset
	# count_occupancy $REFPT/UASv3.bed $SPT7 $TEMP/UASv3_sort-Spt7-$ID-Offset
	# count_occupancy $REFPT/UASv4.bed $SPT7 $TEMP/UASv4_sort-Spt7-$ID-Offset
done

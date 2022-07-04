#!/bin/bash

# This script builds the the Mittal 2022 gene classes

WRK=/path/to/2022-Mittal_2022/02_Create_Gene_Classes
BIN=$WRK/../bin
ROSSI=$WRK/../data/RefPT-YEP
MITTAL=$WRK/../data/RefPT-Mittal
TEMP=$WRK/tmp


#===Initial Inputs===
TSS=$ROSSI/TSS.bed
SUA7_MHS=$TEMP/Sua7_sort-Sua7-12275-Offset-NormalizedCount.bed
SUA7_HS=$TEMP/Sua7_sort-Sua7-26344-Offset-NormalizedCount.bed
SPT7_MHS=$TEMP/STM_sort-Spt7-11960-Offset-NormalizedCount.bed
SPT7_HS=$TEMP/STM_sort-Spt7-20115-Offset-NormalizedCount.bed
TAF2_MHS=$TEMP/Sua7_sort-Taf2-11846-Offset-NormalizedCount.bed
TAF2_HS=$TEMP/Sua7_sort-Taf2-28736-Offset-NormalizedCount.bed


#===Script Shortcuts===
SCRIPTMANAGER=$BIN/ScriptManager-v0.13.jar
FILTERS=$BIN/filter_BED_by_string_ColumnSelect.pl
FILTERV=$BIN/filter_BED_by_value_ColumnSelect.pl
FILTERL=$BIN/filter_BED_by_list_ColumnSelect.pl
CLOSEST=$BIN/determine_closest_RefPoint_output_Both.pl
UPDATEC=$BIN/update_BED_coord.pl
UPDATES=$BIN/update_BED_score_with_TAB_score.pl
RATIO=$BIN/calculate_BED_ScoreRatio.pl
SORT=$BIN/sort_BED_by_Score_v2.pl
SUM=$BIN/sum_Row_CDT.pl
DEDUPLICATE=$BIN/deduplicate_BED_coord_keep_highest_score.py
STM=$BIN/get_STM_from_Rossi_STable1.py

cd $WRK
[ -d $MITTAL ] || mkdir $MITTAL
[ -d $MASTER ] || mkdir $MASTER
[ -d $TEMP ] || mkdir $TEMP
cd $TEMP


#===Build M01 and H01 gene classes===
echo Build M01/H01...

# Filter to keep only 01_RP
perl $FILTERS $TSS 01_RP 4 keep TSS_01-RP.bed

# Get geneID of YEP class=01_RP
cut -f4 TSS_01-RP.bed > 01-RP.tab

# Filter Sua7 bedfile for YEP classes=01_RP
perl $FILTERL $SUA7_MHS 01-RP.tab 3 keep Sua7_RP_sort-Sua7-12275.bed
perl $FILTERL $SUA7_HS 01-RP.tab 3 keep Sua7_RP_sort-Sua7-26344.bed

# Update BED file with Spt7 occupancy score
perl $UPDATES Sua7_RP_sort-Sua7-12275.bed <(cut -f4,5 $SPT7_MHS) Sua7_RP_score-Spt7-11960.bed
perl $UPDATES Sua7_RP_sort-Sua7-26344.bed <(cut -f4,5 $SPT7_HS) Sua7_RP_score-Spt7-20115.bed

# Sort BED file by the occupancy score
perl $SORT Sua7_RP_score-Spt7-11960.bed desc $MITTAL/Sua7_RP_sort-Spt7-11960.bed
perl $SORT Sua7_RP_score-Spt7-20115.bed desc $MITTAL/Sua7_RP_sort-Spt7-20115.bed

# Expand TFIIB/Sua7 coordinates to 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MITTAL/Sua7_RP_sort-Spt7-11960.bed -o $MITTAL/Sua7_RP_sort-Spt7-11960_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $MITTAL/Sua7_RP_sort-Spt7-20115.bed -o $MITTAL/Sua7_RP_sort-Spt7-20115_1000bp.bed


#===Build M02/H02 gene classes===
echo Build M02/H02...

# Filter to remove 01_RP
perl $FILTERL $SUA7_MHS 01-RP.tab 3 remove Sua7_STM-TFO-UNB_sort-Sua7-12275_unfiltered.bed
perl $FILTERL $SUA7_HS 01-RP.tab 3 remove Sua7_STM-TFO-UNB_sort-Sua7-26344-unfiltered.bed

# Remove blacklisted genes (YNR050C)
grep -v 'YNR050C' Sua7_STM-TFO-UNB_sort-Sua7-12275_unfiltered.bed > Sua7_STM-TFO-UNB_sort-Sua7-12275.bed
grep -v 'YNR050C' Sua7_STM-TFO-UNB_sort-Sua7-26344-unfiltered.bed > Sua7_STM-TFO-UNB_sort-Sua7-26344.bed

# Take top 1000 genes
head -n 1000 Sua7_STM-TFO-UNB_sort-Sua7-12275.bed > Sua7_top1000_sort-Sua7-12275.bed
head -n 1000 Sua7_STM-TFO-UNB_sort-Sua7-26344.bed > Sua7_top1000_sort-Sua7-26344.bed

# Update with Taf2 scores
perl $UPDATES Sua7_top1000_sort-Sua7-12275.bed <(cut -f4,5 $TAF2_MHS) Sua7_top1000_score-Taf2-11846.bed
perl $UPDATES Sua7_top1000_sort-Sua7-26344.bed <(cut -f4,5 $TAF2_HS) Sua7_top1000_score-Taf2-28736.bed

# Calculate Sua7/Taf2 ratio
perl $RATIO Sua7_top1000_sort-Sua7-12275.bed Sua7_top1000_score-Taf2-11846.bed Sua7_top1000_score-Sua7Taf2-ratio-MHS.bed
perl $RATIO Sua7_top1000_sort-Sua7-26344.bed Sua7_top1000_score-Taf2-28736.bed Sua7_top1000_score-Sua7Taf2-ratio-HS.bed

# Sort BED file by the ratio score
perl $SORT Sua7_top1000_score-Sua7Taf2-ratio-MHS.bed desc Sua7_top1000_sort-Sua7Taf2-ratio-MHS.bed
perl $SORT Sua7_top1000_score-Sua7Taf2-ratio-HS.bed desc Sua7_top1000_sort-Sua7Taf2-ratio-HS.bed

# Pull M02/H02 classes
head -n 300 Sua7_top1000_sort-Sua7Taf2-ratio-MHS.bed > Sua7_M02_sort-Sua7Taf2-ratio.bed
head -n 300 Sua7_top1000_sort-Sua7Taf2-ratio-HS.bed > Sua7_H02_sort-Sua7Taf2-ratio.bed

# Determine class labels
head -n 150 Sua7_M02_sort-Sua7Taf2-ratio.bed > Sua7_M02a_sort-Sua7Taf2-ratio.bed
tail -n 150 Sua7_M02_sort-Sua7Taf2-ratio.bed > Sua7_M02b_sort-Sua7Taf2-ratio.bed
head -n 150 Sua7_H02_sort-Sua7Taf2-ratio.bed > Sua7_H02a_sort-Sua7Taf2-ratio.bed
tail -n 150 Sua7_H02_sort-Sua7Taf2-ratio.bed > Sua7_H02b_sort-Sua7Taf2-ratio.bed

# Update with Sua7 scores
perl $UPDATES Sua7_M02a_sort-Sua7Taf2-ratio.bed <(cut -f4,5 $SUA7_MHS) Sua7_M02a_score-Sua7-12275.bed
perl $UPDATES Sua7_M02b_sort-Sua7Taf2-ratio.bed <(cut -f4,5 $SUA7_MHS) Sua7_M02b_score-Sua7-12275.bed
perl $UPDATES Sua7_H02a_sort-Sua7Taf2-ratio.bed <(cut -f4,5 $SUA7_HS) Sua7_H02a_score-Sua7-26344.bed
perl $UPDATES Sua7_H02b_sort-Sua7Taf2-ratio.bed <(cut -f4,5 $SUA7_HS) Sua7_H02b_score-Sua7-26344.bed

M02A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M02a__150_SORT-Sua7occ
M02B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M02b__150_SORT-Sua7occ
H02A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H02a__150_SORT-Sua7occ
H02B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H02b__150_SORT-Sua7occ

# Sort by Sua7 occupancy
perl $SORT Sua7_M02a_score-Sua7-12275.bed desc $M02A.bed
perl $SORT Sua7_M02b_score-Sua7-12275.bed desc $M02B.bed
perl $SORT Sua7_H02a_score-Sua7-26344.bed desc $H02A.bed
perl $SORT Sua7_H02b_score-Sua7-26344.bed desc $H02B.bed

# Expand TFIIB/Sua7 coordinates to 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $M02A.bed -o $M02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $M02B.bed -o $M02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $H02A.bed -o $H02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $H02B.bed -o $H02B\_1000bp.bed

#cut -f4 tmp_makerefs/Sua7_Mittal-M02a_sort-Sua7Taf2-ratio.bed Chitvan_BED/M02a_TFIIB_1000bp.bed |sort |uniq -c |awk '{print $1}'  |sort |uniq -c


#===Build M03/H03 gene classes===
echo Build M03/H03...

# Get geneIDs of M02/H02 genes
cut -f4 Sua7_M02_sort-Sua7Taf2-ratio.bed > M02.tab
cut -f4 Sua7_H02_sort-Sua7Taf2-ratio.bed > H02.tab

# Remove M02/H02 from top 1000 Sua7 occupied genes
perl $FILTERL Sua7_STM-TFO-UNB_sort-Sua7-12275.bed M02.tab 3 remove Sua7_minusRP-minusM02_sort-Sua7-12275.bed
perl $FILTERL Sua7_STM-TFO-UNB_sort-Sua7-26344.bed H02.tab 3 remove Sua7_minusRP-minusH02_sort-Sua7-26344.bed

# Update with Spt7 scores
perl $UPDATES Sua7_minusRP-minusM02_sort-Sua7-12275.bed <(cut -f4,5 $SPT7_MHS) Sua7_minusRP-minusM02_score-Spt7-11960.bed
perl $UPDATES Sua7_minusRP-minusH02_sort-Sua7-26344.bed <(cut -f4,5 $SPT7_HS) Sua7_minusRP-minusH02_score-Spt7-20115.bed

# Sort BED file by the Spt7 occupancy score
perl $SORT Sua7_minusRP-minusM02_score-Spt7-11960.bed desc Sua7_minusRP-minusM02_sort-Spt7-11960.bed
perl $SORT Sua7_minusRP-minusH02_score-Spt7-20115.bed desc Sua7_minusRP-minusH02_sort-Spt7-20115.bed

# Shift upstream 250 bp
#perl $SHIFT Sua7_top1000-minusM02_sort-Spt7-11960.bed -250

# Pull M03/H03 classes
head -n 300 Sua7_minusRP-minusM02_sort-Spt7-11960.bed > Sua7_M03_sort-Spt7-11960.bed
head -n 300 Sua7_minusRP-minusH02_sort-Spt7-20115.bed > Sua7_H03_sort-Spt7-20115.bed

# Final class labels
head -n 150 Sua7_M03_sort-Spt7-11960.bed > Sua7_M03a_sort-Spt7-11960.bed
tail -n 150 Sua7_M03_sort-Spt7-11960.bed > Sua7_M03b_sort-Spt7-11960.bed
head -n 150 Sua7_H03_sort-Spt7-20115.bed > Sua7_H03a_sort-Spt7-20115.bed
tail -n 150 Sua7_H03_sort-Spt7-20115.bed > Sua7_H03b_sort-Spt7-20115.bed

# Update with Sua7 scores
perl $UPDATES Sua7_M03a_sort-Spt7-11960.bed <(cut -f4,5 $SUA7_MHS) Sua7_M03a_score-Sua7-12275.bed
perl $UPDATES Sua7_M03b_sort-Spt7-11960.bed <(cut -f4,5 $SUA7_MHS) Sua7_M03b_score-Sua7-12275.bed
perl $UPDATES Sua7_H03a_sort-Spt7-20115.bed <(cut -f4,5 $SUA7_HS) Sua7_H03a_score-Sua7-26344.bed
perl $UPDATES Sua7_H03b_sort-Spt7-20115.bed <(cut -f4,5 $SUA7_HS) Sua7_H03b_score-Sua7-26344.bed

M03A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M03a__150_SORT-Sua7occ
M03B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M03b__150_SORT-Sua7occ
H03A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H03a__150_SORT-Sua7occ
H03B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H03b__150_SORT-Sua7occ

# Sort by Sua7 occupancy
perl $SORT Sua7_M03a_score-Sua7-12275.bed desc $M03A.bed
perl $SORT Sua7_M03b_score-Sua7-12275.bed desc $M03B.bed
perl $SORT Sua7_H03a_score-Sua7-26344.bed desc $H03A.bed
perl $SORT Sua7_H03b_score-Sua7-26344.bed desc $H03B.bed

# Expand TFIIB/Sua7 coordinates to 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $M03A.bed -o $M03A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $M03B.bed -o $M03B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $H03A.bed -o $H03A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $H03B.bed -o $H03B\_1000bp.bed


#===Build M04/H04 gene classes===
echo Build M04/H04...

M04=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M04__150_SORT-Sua7occ
H04=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H04__150_SORT-Sua7occ

# Get geneIDs of M03/H03 genes
cut -f4 Sua7_M03_sort-Spt7-11960.bed > M03.tab
cut -f4 Sua7_H03_sort-Spt7-20115.bed > H03.tab

# Remove M03/H03 from minusRP-minusM/H02 set
perl $FILTERL Sua7_minusRP-minusM02_sort-Sua7-12275.bed M03.tab 3 remove Sua7_minusRP-minusM02-minusM03_sort-Sua7-12275.bed
perl $FILTERL Sua7_minusRP-minusH02_sort-Sua7-26344.bed H03.tab 3 remove Sua7_minusRP-minusH02-minusH03_sort-Sua7-26344.bed

# Update with Taf2 scores
perl $UPDATES Sua7_minusRP-minusM02-minusM03_sort-Sua7-12275.bed <(cut -f4,5 $TAF2_MHS) Sua7_minusRP-minusM02-minusM03_score-Taf2-11846.bed
perl $UPDATES Sua7_minusRP-minusH02-minusH03_sort-Sua7-26344.bed <(cut -f4,5 $TAF2_HS) Sua7_minusRP-minusH02-minusH03_score-Taf2-28736.bed

# Calculate Sua7/Taf2 ratio
perl $RATIO Sua7_minusRP-minusM02-minusM03_sort-Sua7-12275.bed Sua7_minusRP-minusM02-minusM03_score-Taf2-11846.bed Sua7_minusRP-minusM02-minusM03_score-Sua7Taf2-ratio-MHS.bed
perl $RATIO Sua7_minusRP-minusH02-minusH03_sort-Sua7-26344.bed Sua7_minusRP-minusH02-minusH03_score-Taf2-28736.bed Sua7_minusRP-minusH02-minusH03_score-Sua7Taf2-ratio-HS.bed

# Sort BED file by the ratio score
perl $SORT Sua7_minusRP-minusM02-minusM03_score-Sua7Taf2-ratio-MHS.bed desc Sua7_minusRP-minusM02-minusM03_sort-Sua7Taf2-ratio-MHS.bed
perl $SORT Sua7_minusRP-minusH02-minusH03_score-Sua7Taf2-ratio-HS.bed desc Sua7_minusRP-minusH02-minusH03_sort-Sua7Taf2-ratio-HS.bed

# Get 75 (150/2) offset of middle row index
MHS_CT=`wc -l Sua7_minusRP-minusM02-minusM03_sort-Sua7Taf2-ratio-MHS.bed | awk '{print $1/2+75}'`
HS_CT=`wc -l Sua7_minusRP-minusH02-minusH03_sort-Sua7Taf2-ratio-HS.bed | awk '{print $1/2+75}'`

# Pull M04/H04 classes from middle
head -n $MHS_CT Sua7_minusRP-minusM02-minusM03_sort-Sua7Taf2-ratio-MHS.bed | tail -n 150 > $M04.bed
head -n $HS_CT Sua7_minusRP-minusH02-minusH03_sort-Sua7Taf2-ratio-HS.bed | tail -n 150 > $H04.bed

# Expand TFIIB/Sua7 coordinates to 1000bp
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $M04.bed -o $M04\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $H04.bed -o $H04\_1000bp.bed

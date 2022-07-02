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
STM_MHS=$TEMP/STM_sort-Spt7-11960-Offset-NormalizedCount.bed
STM_HS=$TEMP/STM_sort-Spt7-20115-Offset-NormalizedCount.bed
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
perl $UPDATES Sua7_RP_sort-Sua7-12275.bed <(cut -f4,5 $STM_MHS) Sua7_RP_score-Spt7-11960.bed
perl $UPDATES Sua7_RP_sort-Sua7-26344.bed <(cut -f4,5 $STM_HS) Sua7_RP_score-Spt7-20115.bed

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

# Pull Taf2 occupancy scores
cut -f4,5 $TAF2_MHS > Taf2-11846_occupancy.tab
cut -f4,5 $TAF2_HS > Taf2-28736_occupancy.tab

# Update with Taf2 scores
perl $UPDATES Sua7_top1000_sort-Sua7-12275.bed Taf2-11846_occupancy.tab Sua7_top1000_score-Taf2-11846.bed
perl $UPDATES Sua7_top1000_sort-Sua7-26344.bed Taf2-28736_occupancy.tab Sua7_top1000_score-Taf2-28736.bed

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

# Pull Sua7 occupancy scores
cut -f4,5 $SUA7_MHS > Sua7-12275_occupancy.tab
cut -f4,5 $SUA7_HS > Sua7-26344_occupancy.tab

# Update with Sua7 scores
perl $UPDATES Sua7_M02a_sort-Sua7Taf2-ratio.bed Sua7-12275_occupancy.tab Sua7_M02a_score-Sua7-12275.bed
perl $UPDATES Sua7_M02b_sort-Sua7Taf2-ratio.bed Sua7-12275_occupancy.tab Sua7_M02b_score-Sua7-12275.bed
perl $UPDATES Sua7_H02a_sort-Sua7Taf2-ratio.bed Sua7-26344_occupancy.tab Sua7_H02a_score-Sua7-26344.bed
perl $UPDATES Sua7_H02b_sort-Sua7Taf2-ratio.bed Sua7-26344_occupancy.tab Sua7_H02b_score-Sua7-26344.bed


#M02A=Sua7_M02a_sort-Sua7Taf2-ratio
#M02B=Sua7_M02b_sort-Sua7Taf2-ratio
#H02A=Sua7_H02a_sort-Sua7Taf2-ratio
#H02B=Sua7_H02b_sort-Sua7Taf2-ratio

# === Write M02/H02 gene classes ===
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


#===Write M02/H02 with Nucleosome RefPT===
REFPT_NUC=$ROSSI/Nuc.bed

NM02A=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M02a__150_SORT-Sua7occ
NM02B=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M02b__150_SORT-Sua7occ
NH02A=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H02a__150_SORT-Sua7occ
NH02B=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H02b__150_SORT-Sua7occ

# Get Nuc RefPT
perl $UPDATEC $M02A.bed $REFPT_NUC $NM02A.bed
perl $UPDATEC $M02B.bed $REFPT_NUC $NM02B.bed
perl $UPDATEC $H02A.bed $REFPT_NUC $NH02A.bed
perl $UPDATEC $H02B.bed $REFPT_NUC $NH02B.bed

#  Expand BED 200bp from center
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM02A.bed -o $NM02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM02B.bed -o $NM02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH02A.bed -o $NH02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH02B.bed -o $NH02B\_1000bp.bed


#===Write M02/H02 with Upstream Activation Sequence (UAS) RefPT===
REFPT_UAS=$ROSSI/STM.bed

UM02A=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M02a__150_SORT-Sua7occ
UM02B=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M02b__150_SORT-Sua7occ
UH02A=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H02a__150_SORT-Sua7occ
UH02B=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H02b__150_SORT-Sua7occ

# Get UAS RefPT
perl $UPDATEC $M02A.bed $REFPT_UAS $UM02A.bed
perl $UPDATEC $M02B.bed $REFPT_UAS $UM02B.bed
perl $UPDATEC $H02A.bed $REFPT_UAS $UH02A.bed
perl $UPDATEC $H02B.bed $REFPT_UAS $UH02B.bed

#  Expand BED 200bp from center
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM02A.bed -o $UM02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM02B.bed -o $UM02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH02A.bed -o $UH02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH02B.bed -o $UH02B\_1000bp.bed

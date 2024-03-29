#!/bin/bash

# This script rebuilds the inital reference files derived from Rossi et al, 2020. Nature. (Supplementary Table 1)

# Dependencies:
#  - Perl
#  - wget
#  - sed
#  - awk
#  - head/tail


# It is expected that this script execute from /path/to/2022-Mittal_SAGA/data/RefPT-YEP/

SHIFT_BP=150

# File shortcuts
BIN=../../bin
YEPTABLE=../Rossi_2021_Supplementary_Data_1.txt
SUA7CX=../Sua7_CX.bed
TEMP=tmp
ELEVENK=$TEMP/ElevenK_features

# Script shortcuts
SCRIPTMANAGER=$BIN/ScriptManager-v0.13.jar
CLOSEST=$BIN/determine_closest_RefPoint_output_Both.pl
DEDUPLICATE=$BIN/deduplicate_BED_coord_keep_highest_score.py
FILTERL=$BIN/filter_BED_by_list_ColumnSelect.pl
FILTERS=$BIN/filter_BED_by_string_ColumnSelect.pl
FILTERV=$BIN/filter_BED_by_value_ColumnSelect.pl
SHIFT=$BIN/shift_BED_center_v2.pl
SORT=$BIN/sort_BED_by_Score_v2.pl
STM=$BIN/get_STM_from_Rossi_STable1.py
SUM=$BIN/sum_Row_CDT.pl
UPDATEC=$BIN/update_BED_coord.pl
UPDATES=$BIN/update_BED_score_with_TAB_score.pl

[ -d $TEMP ] || mkdir $TEMP


# Download ScriptManager
if [ ! -f $SCRIPTMANAGER ]; then
	wget https://github.com/CEGRcode/scriptmanager/releases/download/v0.13/ScriptManager-v0.13.jar
	mv ScriptManager-v0.13.jar $SCRIPTMANAGER
fi

##--------Build Reference files for 5376 gene feaures (YEP)--------

#===+1 Nucleosome===
# Pull +1Nuc coordinates with Sua7 occ from Rossi Supplementary Table 1 (5873 genes)
sed 1d $YEPTABLE \
  | awk '{FS="\t"}{OFS="\t"}{print $1,$21,$21,$7,$41,$2}' \
  | head -n 5378 \
  > Nuc.bed

#===Nucleosome Free Region (NFR)===
# Pull NFR coordinates with Sua7 occ from Rossi Supplementary Table 1 (5873 genes)
sed 1d $YEPTABLE \
  | awk '{FS="\t"}{OFS="\t"}{print $1,$22,$23,$7,$41,$2}' \
  | head -n 5378 \
  > NFR.bed
# Expand 50bp window
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 50 NFR.bed -o NFR_50bp.bed

# Subset NFR into 04-only, 03-only, and 04+03 ref files
tail -n 4257 NFR_50bp.bed > NFR_03-04_50bp.bed
head -n 1783 NFR_03-04_50bp.bed > NFR_03_50bp.bed
tail -n 2474 NFR_03-04_50bp.bed > NFR_04_50bp.bed

#===TSS===
# Get TSS coordinates from the Rossi 2021 Supplementary Table 1 for the 01_RP, 02_STM, 03_TFO, and 04_UNB gene classes (5873 genes)
sed 1d $YEPTABLE \
  | awk '{FS="\t"}{OFS="\t"}{print $1,$15,$15,$7,$4,$2}' \
  | head -n 5378 \
  > TSS.bed

#===TSS_shifted===
perl $SHIFT TSS.bed 250 TSS+250.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed TSS+250.bed -c 500 -o TSS_0-500.bed

#===Sua7===
# Associate TSS to a Sua7 ChexMix peak (downloaded from Rossi 2021 Github)
perl $CLOSEST TSS.bed $SUA7CX $TEMP/TSS_closestSua7.bed
# Column select for new BED file keeping Sua7 coords, geneIDs, and dist scores
awk '{FS="\t"}{OFS="\t"}{print $1,$8,$9,$4,$13,$6}' $TEMP/TSS_closestSua7.bed > $TEMP/Sua7_score-distTSS.bed
# Filter by distance (upstream 100 and downstream 61)
perl $FILTERV $TEMP/Sua7_score-distTSS.bed -100 4 keep $TEMP/Sua7_filter100upstream.bed
perl $FILTERV $TEMP/Sua7_filter100upstream.bed 61 4 remove $TEMP/Sua7_filter.bed
# Deduplicate Sua7 peaks
python $DEDUPLICATE -i $TEMP/Sua7_filter.bed -o $TEMP/Sua7_deduplicate.bed
# Get identifiers of genes with Sua7 peaks
cut -f4 $TEMP/Sua7_deduplicate.bed > $TEMP/geneswpeaks.tab
# Subtract out genes with Sua7 peak from TSS coordinate file
perl $FILTERL TSS.bed $TEMP/geneswpeaks.tab 3 remove $TEMP/TSS_geneswopeaks.bed
# Update BED scores with default "imputed" value
awk '{FS="\t"}{OFS="\t"}{print $1,$2,$3,$4,"imputed",$6}' $TEMP/TSS_geneswopeaks.bed > $TEMP/TSS_imputed.bed
# Merge genes with & missing a Sua7_CX peak
cat $TEMP/TSS_imputed.bed $TEMP/Sua7_deduplicate.bed > Sua7.bed

#===STMv1===
# Use custom script to make Spt7 ref center
python $STM -i $YEPTABLE -p $TEMP/Sua7_deduplicate.bed -o STM.bed -s $SHIFT_BP

#===UASv2===
echo "UASv2"
# Identify manually annotated UAS features
cat ../Rossi_2021_Supplementary_Data_1.txt \
 | awk '{FS="\t";OFS="_"}{if($35!="" && $35!=".") print $35,$2,$7}' \
 | awk '{FS="_"}{OFS="\t"}{print $3,$4,$5,$7,$1"_"$2,$6}' \
 | sed 1d > $TEMP/UAS_manual.bed
# Create a -150bp (upstream shift) version of the Sua7 annotations
perl $SHIFT Sua7.bed -150 $TEMP/Sua7_upstream150.bed
# Filter out genes with a UAS annotation from shifted Sua7 set
perl $FILTERL $TEMP/Sua7_upstream150.bed <(cut -f4 $TEMP/UAS_manual.bed) 3 remove $TEMP/Sua7_upstream150_noManualUAS.bed
# Merge BED of shifted Sua7 & manual UAS
cat $TEMP/UAS_manual.bed $TEMP/Sua7_upstream150_noManualUAS.bed > $TEMP/UAS_unsort.bed
# Update coordinates to sorted BED
perl $UPDATEC TSS.bed $TEMP/UAS_unsort.bed $TEMP/UAS_unlabeled.bed
# Update coordinates to sorted BED
perl $UPDATES $TEMP/UAS_unlabeled.bed <(cut -f4,5 $TEMP/UAS_unsort.bed) UASv2.bed
# Expand coords for occupancy calculations (200bp window)
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed UASv2.bed -c 200 -o UASv2_200bp.bed

#===UASv3===
echo "UASv3"
# Identify manually annotated UAS features
head -n 5379 ../Rossi_2021_Supplementary_Data_1.txt \
  | awk '{FS="\t";OFS="\t"}{print $1,$32,$32,$7,$30,$2}' \
  | sed 1d > UASv3.bed


#===UASv4===
echo "UASv4"
# Get list of 04_UNB genes
perl $FILTERS TSS.bed "04_UNB" 4 keep $TEMP/TSS_04-UNB.bed
# Update with Sua7 coord
perl $UPDATEC $TEMP/TSS_04-UNB.bed Sua7.bed $TEMP/Sua7_04-UNB.bed
# Shift upstream 150
perl $SHIFT $TEMP/Sua7_04-UNB.bed -150 $TEMP/Sua7_04-UNB_upstream150.bed
# Filter 04_UNB from UASv3
perl $FILTERL UASv3.bed <(cut -f4 $TEMP/TSS_04-UNB.bed) 3 remove $TEMP/UASv3_minus04-UNB.bed
# Append shifted Sua7 04_UNB RefPT to fill out filtered UASv3 set
cat $TEMP/UASv3_minus04-UNB.bed $TEMP/Sua7_04-UNB_upstream150.bed > UASv4.bed


##--------Build Reference files for 11K gene feaures (YEP)--------

#===+1Nuc===
sed 1d $YEPTABLE \
  | awk '{FS="\t"}{OFS="\t"}{if($4=="01_RP" || $4=="02_STM" || $4=="03_TFO" || $4=="04_UNB") print $1,$21,$21,$7,$4,$2; else print $1,$21,$21,$3,$4,$2}' \
  | awk '{FS="\t"}{OFS="\t"}{if($4!="80021" && $4!="160127" && $4!="240035" && $4!="240036") print }' \
  > ElevenK_features_+1Nuc.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed ElevenK_features_+1Nuc.bed -c 200 -o ElevenK_features_+1Nuc_200bp.bed

#===Sua7===
# Create Sua7 ref in STable 1 feature order
perl $UPDATEC TSS.bed Sua7.bed $ELEVENK\_Sua7_core5378.bed
# Create TSS ref for remaining STable 1 features
tail -n 5734 $YEPTABLE \
  | awk '{FS="\t"}{OFS="\t"}{print $1, $15, $15, $3, $4, $2}' \
   > $ELEVENK\_TSS_subtract5378_unfiltered.bed
# Filter out features that don't expand outside chromosome
awk '{FS="\t"}{OFS="\t"}{if($4!="80021" && $4!="160127" && $4!="240035" && $4!="240036") print }' $ELEVENK\_TSS_subtract5378_unfiltered.bed > $ELEVENK\_TSS_subtract5378.bed
# Merge 5378 ordered genes and rest of 11k
cat $ELEVENK\_Sua7_core5378.bed $ELEVENK\_TSS_subtract5378.bed > ElevenK_features_Sua7.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed ElevenK_features_Sua7.bed -c 200 -o ElevenK_features_Sua7_200bp.bed

#===STM===
# Create STM ref in STable 1 feature order
perl $UPDATEC TSS.bed STM.bed $ELEVENK\_STM_core5378.bed
# Use TSS ref and shift by
perl $SHIFT $ELEVENK\_TSS_subtract5378.bed -150 $ELEVENK\_TSS-shiftedby$SHIFT_BP\_subtract5378.bed
# Merge 5378 ordered genes and rest of 11k
cat $ELEVENK\_STM_core5378.bed $ELEVENK\_TSS-shiftedby$SHIFT_BP\_subtract5378.bed > ElevenK_features_STM.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed ElevenK_features_STM.bed -c 200 -o ElevenK_features_STM_200bp.bed


## Clean-up
rm -r $TEMP

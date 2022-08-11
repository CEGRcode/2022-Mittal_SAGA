#!/bin/bash

# This script rebuilds the inital reference files derived from Rossi et al, 2020. Nature. (Supplementary Table 1)

# Dependencies:
#  - Perl
#  - wget
#  - sed
#  - awk
#  - head/tail


# It is expected that this script execute from /path/to/2022-Mittal_SAGA/data/RefPT-Rossi/

# File shortcuts
BIN=../../bin
YEPTABLE=../Rossi_2021_Supplementary_Data_1.txt
SUA7CX=../Sua7_CX.bed
TEMP=tmp
ELEVENK=$TEMP/ElevenK_features
DYEP=../RefPT-YEP

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

##--------Build Reference files for 5378 gene feaures (YEP-lvl 1)--------

# Subset 5378 annotations by each of the four main Rossi gene classes using the TSS coordinate RefPT
for level_one in "01_RP" "02_STM" "03_TFO" "04_UNB";
do
	sed 1d $YEPTABLE \
	  | awk -v group="$level_one" '{FS="\t"}{OFS="\t"}{if($4==group) print $1,$15,$15,$7,$41,$2}' \
	  > $TEMP/tmp.bed
	# sort -nk1,1 | cut -f2-7 > $level_one\_ORF.bed
	LNCT=`wc -l $TEMP/tmp.bed | awk '{print $1}'`
	TYPE=`echo $level_one | sed 's/_/-/g'`
	OFFICIAL=FEAT-Pol-II_RefPT-TSS__SubFEAT-YPD-25-37C-coding_$TYPE\___$LNCT\_UNSORTED
	echo $OFFICIAL
	mv $TEMP/tmp.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed
done

# Create the +1Nuc & Sua7 & STMv1 & UASv2 versions of these Reference files
for BED in `ls FEAT*RefPT-TSS__*coding*_UNSORTED.bed`;
do
	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT-Sua7_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/Sua7.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed

	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT+1Nuc_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/Nuc.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed

	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT-STM_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/STM.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed

	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT-UASv2_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/UASv2.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed
done


exit

# Subset 5378 annotations - LEVEL TWO - RefPT-TSS
for level_two in "02a_STM_RSTM" "02b_STM_YPD" "02c_STM_HS" "03a_TF" "03b_COF" "04a_RSC_HS" "04b_UNB";
do
	sed 1d $YEPTABLE \
	  | awk -v group="$level_two" '{FS="\t"}{OFS="\t"}{if($4==group) print $1,$15,$15,$7,group,$2}' \
	  > $TEMP/tmp.bed
	# sort -nk1,1 | cut -f2-7 > $level_one\_ORF.bed
	LNCT=`wc -l tmp.bed | awk '{print $1}'`
	mv $TEMP/tmp.bed FEAT-Pol-II_RefPT-TSS__SubFEAT-YPD-25-37C-coding_$level_two\___$LNCT\_1000bp.bed
done
# RefPT+1Nuc & RefPT-Sua7 & RefPT-UAS
for BED in `ls FEAT*RefPT-TSS__*coding*_UNSORTED.bed`;
do
	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT-Sua7_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/Sua7.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed

	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT+1Nuc_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/Nuc.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed

	OFFICIAL=`basename $BED ".bed"| sed "s/_RefPT-TSS_/_RefPT-UAS_/g" `
	echo $OFFICIAL
	perl $UPDATEC $BED $DYEP/UAS.bed $OFFICIAL.bed
	java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 1000 $OFFICIAL.bed -o $OFFICIAL\_1000bp.bed
done

##--------Build Reference files for remaining of 11K gene feaures (YEP)--------

#14_LTR		delta	9999	SORT-Sua7occ	FEAT-LTR_RefPT-Start____delta__9999_SORT-Sua7occ.bed
#14_LTR		omega	9999	SORT-Sua7occ	FEAT-LTR_RefPT-Start____omega__9999_SORT-Sua7occ.bed
#14_LTR		sigma	9999	SORT-Sua7occ	FEAT-LTR_RefPT-Start____sigma__9999_SORT-Sua7occ.bed
#14_LTR		tau	9999	SORT-Sua7occ	FEAT-LTR_RefPT-Start____tau__9999_SORT-Sua7occ.bed
#14_LTR		ALL	9999	SORT-Sua7occ	FEAT-LTR_RefPT-Start____ALL__9999_SORT-Sua7occ.bed

## LTR LTR Rossi 2019 STable 1
for LTR in "LTR_delta" "LTR_omega" "LTR_sigma" "LTR_tau";
do
	sed 1d $YEPTABLE \
	  | awk -v group="$LTR" '{FS="\t"}{OFS="\t"}{if($4=="14_LTR" && $5==group) print $1,$11,$11,$7,$41,$2}' \
	  > $TEMP/tmp.bed
	# | sort -rnk5,5
	TYPE=`echo $LTR | sed 's/LTR_//g'`
	LNCT=`wc -l tmp.bed | awk '{print $1}'`
	mv $TEMP/tmp.bed FEAT-Pol-II_RefPT-Start__SubFEAT-YPD-25-37C-coding_$TYPE\___$LNCT\_1000bp.bed
done
# Full LTR list
awk '{FS="\t"}{OFS="\t"}{if($4=="14_LTR") print $1,$11,$11+1,$7,$41,$2}' Rossi_2021_Supplementary_Data_1.txt \
	| sort -rnk5,5 > temp.bed
LNCT=`wc -l temp.bed |awk '{print $1}'`
mv temp.bed FEAT-LTR_RefPT-Start____ALL__$LNCT\_SORT-Sua7occ.bed


#16_ACS	16-ACS		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___16-ACS___9999_SORT-Length.bed
#17_CEN	17-CEN		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___17-CEN___9999_SORT-Length.bed
#18_MAT	18-MAT		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___18-MAT___9999_SORT-Length.bed
#19_XEL	19-XEL		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___19-XEL___9999_SORT-Length.bed
#20_YEL	20-YEL		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___20-YEL___9999_SORT-Length.bed
#21_TEL	21-TEL		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___21-TEL___9999_SORT-Length.bed
#23_TMR	23-TMR		9999	SORT-Length	FEAT-NonTranscribed_RefPT-Start___23-TMR___9999_SORT-Length.bed

## Other Features from Rossi 2019 STable 1
for GROUP in "16_ACS" "17_CEN" "18_MAT" "19_XEL" "20_YEL" "21_TEL" "23_TMR";
do
	awk -v group="$GROUP" '{FS="\t"}{OFS="\t"}{if($4==group) print $1,$11,$11+1,$7,$13,$2}' Rossi_2021_Supplementary_Data_1.txt \
		| sort -rnk5,5 > temp.bed
	TYPE=`echo $GROUP | sed 's/_/-/g'`
	LNCT=`wc -l temp.bed |awk '{print $1}'`
	mv temp.bed FEAT-NonTranscribed_RefPT-Start___$TYPE\___$LNCT\_SORT-Length.bed
done

## Clean-up
rm -r $TEMP

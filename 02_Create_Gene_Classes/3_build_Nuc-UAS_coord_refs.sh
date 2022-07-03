# This script builds the Mittal gene group (M02-M04,H02-H04) Nucleosome and Upstream Activation Sequence (UAS) reference BED files

WRK=/path/to/2022-Mittal_2022/02_Create_Gene_Classes
BIN=$WRK/../bin
ROSSI=$WRK/../data/RefPT-YEP
MITTAL=$WRK/../data/RefPT-Mittal

#===Initial Inputs===
M02A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M02a__150_SORT-Sua7occ
M02B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M02b__150_SORT-Sua7occ
M03A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M03a__150_SORT-Sua7occ
M03B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M03b__150_SORT-Sua7occ
M04=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-25C_M04__150_SORT-Sua7occ
H02A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H02a__150_SORT-Sua7occ
H02B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H02b__150_SORT-Sua7occ
H03A=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H03a__150_SORT-Sua7occ
H03B=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H03b__150_SORT-Sua7occ
H04=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT-37C_H04__150_SORT-Sua7occ

#===Script Shortcuts===
SCRIPTMANAGER=$BIN/ScriptManager-v0.13.jar
UPDATEC=$BIN/update_BED_coord.pl

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

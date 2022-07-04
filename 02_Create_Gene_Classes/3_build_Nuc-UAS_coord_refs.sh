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
H02ANOTM02AB=$MITTAL/FEAT-Pol-II_RefPT-Sua7___SubFEAT_Induced__72_SORT-Sua7occ

#===Script Shortcuts===
SCRIPTMANAGER=$BIN/ScriptManager-v0.13.jar
UPDATEC=$BIN/update_BED_coord.pl

#===Write M02/H02 with Nucleosome RefPT===
REFPT_NUC=$ROSSI/Nuc.bed

NM02A=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M02a__150_SORT-Sua7occ
NM02B=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M02b__150_SORT-Sua7occ
NM03A=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M03a__150_SORT-Sua7occ
NM03B=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M03b__150_SORT-Sua7occ
NM04=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-25C_M04__150_SORT-Sua7occ
NH02A=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H02a__150_SORT-Sua7occ
NH02B=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H02b__150_SORT-Sua7occ
NH03A=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H03a__150_SORT-Sua7occ
NH03B=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H03b__150_SORT-Sua7occ
NH04=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT-37C_H04__150_SORT-Sua7occ
NH02ANOTM02AB=$MITTAL/FEAT-Pol-II_RefPT+1Nuc___SubFEAT_Induced__72_SORT-Sua7occ

# Get Nuc RefPT
perl $UPDATEC $M02A.bed $REFPT_NUC $NM02A.bed
perl $UPDATEC $M02B.bed $REFPT_NUC $NM02B.bed
perl $UPDATEC $M03A.bed $REFPT_NUC $NM03A.bed
perl $UPDATEC $M03B.bed $REFPT_NUC $NM03B.bed
perl $UPDATEC $M04.bed $REFPT_NUC $NM04.bed
perl $UPDATEC $H02A.bed $REFPT_NUC $NH02A.bed
perl $UPDATEC $H02B.bed $REFPT_NUC $NH02B.bed
perl $UPDATEC $H03A.bed $REFPT_NUC $NH03A.bed
perl $UPDATEC $H03B.bed $REFPT_NUC $NH03B.bed
perl $UPDATEC $H04.bed $REFPT_NUC $NH04.bed
perl $UPDATEC $H02ANOTM02AB.bed $REFPT_NUC $NH02ANOTM02AB.bed

#  Expand BED 200bp from center
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM02A.bed -o $NM02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM02B.bed -o $NM02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM03A.bed -o $NM03A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM03B.bed -o $NM03B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NM04.bed -o $NM04\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH02A.bed -o $NH02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH02B.bed -o $NH02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH03A.bed -o $NH03A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH03B.bed -o $NH03B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH04.bed -o $NH04\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $NH02ANOTM02AB.bed -o $NH02ANOTM02AB\_1000bp.bed


#===Write M02/H02 with Upstream Activation Sequence (UAS) RefPT===
REFPT_UAS=$ROSSI/STM.bed

UM02A=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M02a__150_SORT-Sua7occ
UM02B=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M02b__150_SORT-Sua7occ
UM03A=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M03a__150_SORT-Sua7occ
UM03B=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M03b__150_SORT-Sua7occ
UM04=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-25C_M04__150_SORT-Sua7occ
UH02A=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H02a__150_SORT-Sua7occ
UH02B=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H02b__150_SORT-Sua7occ
UH03A=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H03a__150_SORT-Sua7occ
UH03B=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H03b__150_SORT-Sua7occ
UH04=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT-37C_H04__150_SORT-Sua7occ
UH02ANOTM02AB=$MITTAL/FEAT-Pol-II_RefPT-UAS___SubFEAT_Induced__72_SORT-Sua7occ

# Get UAS RefPT
perl $UPDATEC $M02A.bed $REFPT_UAS $UM02A.bed
perl $UPDATEC $M02B.bed $REFPT_UAS $UM02B.bed
perl $UPDATEC $M03A.bed $REFPT_UAS $UM03A.bed
perl $UPDATEC $M03B.bed $REFPT_UAS $UM03B.bed
perl $UPDATEC $M04.bed $REFPT_UAS $UM04.bed
perl $UPDATEC $H02A.bed $REFPT_UAS $UH02A.bed
perl $UPDATEC $H02B.bed $REFPT_UAS $UH02B.bed
perl $UPDATEC $H03A.bed $REFPT_UAS $UH03A.bed
perl $UPDATEC $H03B.bed $REFPT_UAS $UH03B.bed
perl $UPDATEC $H04.bed $REFPT_UAS $UH04.bed
perl $UPDATEC $H02ANOTM02AB.bed $REFPT_UAS $UH02ANOTM02AB.bed

#  Expand BED 200bp from center
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM02A.bed -o $UM02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM02B.bed -o $UM02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM03A.bed -o $UM03A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM03B.bed -o $UM03B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UM04.bed -o $UM04\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH02A.bed -o $UH02A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH02B.bed -o $UH02B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH03A.bed -o $UH03A\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH03B.bed -o $UH03B\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH04.bed -o $UH04\_1000bp.bed
java -jar $SCRIPTMANAGER coordinate-manipulation expand-bed -c 100 $UH02ANOTM02AB.bed -o $UH02ANOTM02AB\_1000bp.bed

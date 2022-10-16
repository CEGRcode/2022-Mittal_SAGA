

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/2022-Mittal_SAGA/04_Figures/F1
cd $WRK

LIBRARY=$WRK/../../LIBRARY/

while read LINE; do
	DIR=`echo $LINE | awk '{print $1}'`
	SAMPLE_SET=`echo $LINE | awk '{print $1}'`
	BED=`echo $LINE | awk '{print $1}'`

	[ -d $DIR ] || mkdir $DIR

	for COMPOSITE in `ls  $LIBRARY/$BED/NormComposites/* |grep -f <(awk '{print "/NormComposites/"$1"_"}' $SAMPLE_SET)`;
		cp $COMPOSITE $DIR
done <subfigure_info.txt

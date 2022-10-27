

# FIRST CHANGE PATH TO EXECUTE
WRK=/path/to/2022-Mittal_SAGA/04_Figures
cd $WRK

LIBRARY=$WRK/../LIBRARY/

for FIG in "F1" "F2" "F3" "F4" "F5" "F6";
do
	cd $WRK/$FIG
	while read LINE; do
		DIR=`echo $LINE | awk '{print $1}'`
		SAMPLE_SET=`echo $LINE | awk '{print $2}'`
		BED=`echo $LINE | awk '{print $3}'`

		[ -d $DIR ] || mkdir $DIR

		for COMPOSITE in `ls  $LIBRARY/$BED/NormComposites/* |grep -f <(awk '{print "/NormComposites/"$1"_"}' input/$SAMPLE_SET)`;
		do
			cp $COMPOSITE $DIR
		done
	done <subfigure_info.txt
done

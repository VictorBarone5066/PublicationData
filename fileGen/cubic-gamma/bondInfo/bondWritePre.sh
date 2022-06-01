#!/bin/bash

if [ -f bondInfoInit.csv ]; then
	rm bondInfoInit.csv
fi

i=0
while [ "$i" -lt 12870 ]; do
	echo "Working on $i..."

	#Get nth binary representation
	j=`echo "$i + 1" | bc -l`
        thisBin=`head -"$j" AgBiI4_Seeds | tail -1`

	#Create POSCAR file
	./AgBiI4_Maker.py -s AgBiI4_Parent -n poscarTmp -g "$i" -b "$thisBin"  

	#Describe it
	./bi -i poscarTmp -o bondInfoInit.csv -b bondDataInfo

	((i++))
done

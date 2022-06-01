#!/bin/bash

if [ -f bondInfoIni.csv ]; then
	rm bondInfoIni.csv
fi
if [ -f bondInfoFin.csv ]; then
	rm bondInfoFin.csv
fi

i=0
while [ "$i" -lt 12870 ]; do
	echo "Working on $i..."

	#Describe it
	./bi -i all/$i/POSCAR -o bondInfoIni.csv -b bondDataInfo
	./bi -i all/$i/CONTCAR -o bondInfoFin.csv -b bondDataInfo

	((i++))
done

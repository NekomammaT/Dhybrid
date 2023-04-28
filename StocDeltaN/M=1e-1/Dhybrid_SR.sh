#!/bin/sh

#rm -f data/maxPdata.dat
while read line
do
    ./Dhybrid_SR $line
done < ./EoI_5th.dat

#rm -f data/maxPdata.dat
#./Dhybrid_SR 50 1 1.667746503534577577277766e-8
#./Dhybrid_SR 2000 1 2.383344659324451720350286e-7


<<COMMENTOUT
rm -f data/maxPdata.dat
for Pi2 in 50 100 200 500 1000 2000
do
    for DW in 1 2 5 10 20 50 100 200 500 1000
    do
	./Dhybrid_SR $Pi2 $DW
    done
done
COMMENTOUT



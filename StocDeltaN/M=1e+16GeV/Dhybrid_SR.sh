#!/bin/sh

#rm -f data/maxPdata.dat
while read line
do
    ./Dhybrid_SR $line
done < ./EoI_3rd.dat

#rm -f data/maxPdata.dat
#./Dhybrid_SR 2000 1000 1.124616056335486036962e-15



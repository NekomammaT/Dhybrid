#1/bin/sh

./Dhybrid_SR 2000 5
./Dhybrid_SR 2000 20


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



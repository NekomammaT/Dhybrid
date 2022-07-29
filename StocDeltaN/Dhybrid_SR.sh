#1/bin/sh

./Dhybrid_SR 50 1 3.0692e-7
./Dhybrid_SR 50 10 2.8742643286367e-6
./Dhybrid_SR 50 100 3.3986231385252e-6
./Dhybrid_SR 50 1000 4.1104691899188e-6


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



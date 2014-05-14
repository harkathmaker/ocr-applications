
#!/bin/bash

#declare -a db_size=(100 1000 10000 100000 1000000 10000000)
db_size=()
for (( i = 100; i < 1000000000; i *= 10 ))
do
	db_size+=($i)
done

for i in "${db_size[@]}"
do
	echo "$i"
done

# set variables
iterations=30
declare -a split=(1 2 5 10 100 1000 10000)
ss_name=ss.csv
ps_name=ps.csv
pds_name=pds.csv
pdes_name=pdes.csv

# iterate through array
#for i in "${db_size[@]}"
#do
	#./serial_stream                     -d "$i" -e "$ss_name"   -i "$iterations" -r
	#./parallel_stream                   -d "$i" -e "$ps_name"   -i "$iterations" -p -r
	#./parallel_db_breakdown_stream      -d "$i" -e "$pds_name"  -i "$iterations" -p -r
	#./parallel_db_edt_breakdown_stream  -d "$i" -e "$pdes_name" -i "$iterations" -p -r
#done


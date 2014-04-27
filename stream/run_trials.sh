#!/bin/bash

#make all

min=100
max=1000000
i=100

# serial stream runs
for name in serial_stream serial_edt_breakdown_stream
	do
		for ((db=$min; db<=$max; db=db*10))
			do 
				echo ./$name -d $db -i $i -e $name.csv -r
				./$name -d $db -i $i -e $name.csv -r
			done
	done

# parallel stream runs
for name in parallel_stream parallel_db_breakdown_stream parallel_db_edt_breakdown_stream
	do
		for ((db=$min; db<=$max; db=db*10))
			do
				for ((p=1; p<=$db; p=p*10)) 
					do
						echo ./$name -d $db -i $i -e $name.csv -p $p -r
						./$name -d $db -i $i -e $name.csv -p $p -r
					done
			done
	done

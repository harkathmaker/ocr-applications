#!/bin/bash

#make all

min=1
#2^24
max=16777216
i=100

# serial stream runs
for name in serial_stream serial_edt_breakdown_stream
	do
		for ((db=$min; db<=$max; db=db*2))
			do 
				echo ./$name -d $db -i $i -e $name.csv -r
				#./$name -d $db -i $i -e $name.csv -r
			done
	done

# parallel stream runs
for name in parallel_stream parallel_db_breakdown_stream parallel_db_edt_breakdown_stream
	do
		for ((db=$min; db<=$max; db=db*2))
			do
				for ((p=$min; p<=$db; p=p*2)) 
					do
						echo ./$name -d $db -i $i -e $name.csv -p $p -r
						#./$name -d $db -i $i -e $name.csv -p $p -r
					done
			done
	done

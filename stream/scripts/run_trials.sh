#!/bin/bash

make all
min=100
max=10000000
i=30

# native STREAM runs
for ((db=$min; db<=$max; db=db*2))
	do 
		echo ./$name -d $db -i $i -e $name.csv -r
		./$name -d $db -i $i -e $name.csv -r
	done

# OCR STREAM runs
'
name=parallel_stream
for p in 1 2 3 4
	do
		for ((db=$min; db<=$max; db=db*10))
			do
				echo ./$name -d $db -i $i -e $name.csv -p $p -r
				./$name -d $db -i $i -e $name.csv -p $p -r
			done

	done
'

#!/bin/bash
make all
min=100
max=10000000
itr=30

name=native_stream
threads=4
# native STREAM runs
for ((db=$min; db<=$max; db=db*2))
	do
		echo gcc $name.c -fopenmp -DSTREAM_ARRAY_SIZE=$db -DNUMTHREADS=$threads -DNTIMES=$itr -o $name
		gcc $name.c -fopenmp -DSTREAM_ARRAY_SIZE=$db -DNUMTHREADS=$threads -DNTIMES=$itr -o $name
		./$name -e $name.csv
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

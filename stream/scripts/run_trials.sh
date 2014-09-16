#!/bin/bash

# Global Controls
make all
min=100
max=10000000
itr=30
inc=2

minthreads=1
maxthreads=4

native=native_stream
nativeparalllel=native_parallel_stream
ocr=parallel_stream

# native STREAM Runs
for ((t =$minthreads; t <=$maxthreads; t++))
	do
		for ((db=$min; db<=$max; db=db*inc))
			do
				echo gcc $native.c -fopenmp -DSTREAM_ARRAY_SIZE=$db -DNUMTHREADS=$t -DNTIMES=$itr -o $nativeparalllel
				gcc $native.c -fopenmp -DSTREAM_ARRAY_SIZE=$db -DNUMTHREADS=$t -DNTIMES=$itr -o $nativeparalllel
				./$nativeparalllel -e $nativeparalllel.csv
			done
	done

# OCR STREAM Runs
for ((p=$minthreads; p <=$maxthreads; p++))
	do
		for ((db=$min; db<=$max; db=db*inc))
			do
				echo ./$ocr -d $db -i $itr -e $ocr.csv -p $p -r
				./$ocr -d $db -i $itr -e $ocr.csv -p $p -r
			done
	done
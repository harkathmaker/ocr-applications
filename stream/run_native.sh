#!/bin/bash

#make all

min=1
#2^24
max=16777216
tmin=1
tmax=16
i=100

# serial stream runs
for name in native_serial_stream
	do
		for ((db=$min; db<=$max; db=db*2))
			do 
				gcc native_stream.c -DSTREAM_ARRAY_SIZE=$db -DNTIMES=$i -o $name
				./$name nss.csv 
				echo yes
			done
	done

# parallel stream runs
for name in native_parallel_stream 
	do
		for ((db=$min; db<=$max; db=db*2))
			do
				for ((t=$tmin; t<=$tmax; t=t*2)) 
					do
						gcc native_stream.c -DSTREAM_ARRAY_SIZE=$db -DNUMTHREADS=$t -DNTIMES=$i -fopenmp -o $name
						./$name nps.csv
					done
			done
	done

#!/bin/bash

make all

echo "PARALLEL STREAM TESTS"
# TEST -h option
./parallel_stream -h

# TEST other options
./parallel_stream -d 100 -i 5 -p 2 -e parallel_stream.csv -r -v

echo "SERIAL STREAM TESTS"
# TEST -h option
./serial_stream -h

# TEST -p option is disabled
./serial_stream -p 2

# TEST other options
./serial_stream -d 100 -i 5 -e serial_stream.csv -r -v


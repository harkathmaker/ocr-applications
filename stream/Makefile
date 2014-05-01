CC=gcc
CFLAGS=-I${OCR_ROOT}/include -g
LDFLAGS=-L${OCR_ROOT}/lib -locr

ifdef mic
        CC := icc
        CFLAGS := $(CFLAGS) -mmic
        LDFLAG := $(LDFLAGS) -mmic -L/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/mic -limf -lsvml -lirng -lintlc
endif

all : native_serial_stream.o native_parallel_stream.o serial_stream.o serial_edt_breakdown_stream.o parallel_stream.o\
      parallel_db_breakdown_stream.o parallel_db_edt_breakdown_stream.o 

native_serial_stream.o : native_stream.c
	$(CC) native_stream.c -o native_serial_stream

native_parallel_stream.o : native_stream.c
	$(CC) native_stream.c -fopenmp -o native_parallel_stream

serial_stream.o : serial_stream.c
	$(CC) $(CFLAGS) $(LDFLAGS) serial_stream.c -o serial_stream

serial_edt_breakdown_stream.o : serial_edt_breakdown_stream.c
	$(CC) $(CFLAGS) $(LDFLAGS) serial_edt_breakdown_stream.c -o serial_edt_breakdown_stream

parallel_stream.o : parallel_stream.c
	$(CC) $(CFLAGS) $(LDFLAGS) parallel_stream.c -o parallel_stream

parallel_db_breakdown_stream.o : parallel_db_breakdown_stream.c
	$(CC) $(CFLAGS) $(LDFLAGS) parallel_db_breakdown_stream.c -o parallel_db_breakdown_stream

parallel_db_edt_breakdown_stream.o : parallel_db_edt_breakdown_stream.c
	$(CC) $(CFLAGS) $(LDFLAGS) parallel_db_edt_breakdown_stream.c -o parallel_db_edt_breakdown_stream

clean :
	rm -rf native_serial_stream native_parallel_stream serial_stream serial_edt_breakdown_stream parallel_stream\
		   parallel_db_breakdown_stream parallel_db_edt_breakdown_stream

PROG=hpcg
PTHREAD_DIR=/usr/lib/x86_64-linux-gnu/
CFLAGS=-O0 -g -Werror -I. -I${OCR_INSTALL}/../../inc -c
OCR_FLAGS=-L${OCR_INSTALL}/lib -L${PTHREAD_DIR} -I${OCR_INSTALL}/include -locr -lpthread -static

ifdef mic
        CC := icc
        CFLAGS := $(CFLAGS) -mmic
        LDFLAG := $(LDFLAGS) -mmic -L/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/mic -limf -lsvml -lirng -lintlc
endif

ifndef OCR_INSTALL
$(error OCR_INSTALL not set)
endif

ifndef OCR_CONFIG
OCR_CONFIG=${OCR_INSTALL}/config/default.cfg
$(warning OCR_CONFIG not set, defaulting to ${OCR_CONFIG})
endif

OCR_RUN_FLAGS=-ocr:cfg ${OCR_CONFIG}

all-test: compile run

compile:
	g++ $(CFLAGS) timer.cpp main.cpp matrix.cpp
	g++ -o $(PROG).exe timer.o main.o matrix.o $(OCR_FLAGS)

run:
#	octave ./cgmatrix.m
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg50.mat 50
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg100.mat 100
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg150.mat 150
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg200.mat 200
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg250.mat 250
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg300.mat 300
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg350.mat 350
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg400.mat 400
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg450.mat 450
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg500.mat 500
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg50.mat 50
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg100.mat 100
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg150.mat 150
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg200.mat 200
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg250.mat 250
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg300.mat 300
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg350.mat 350
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg400.mat 400
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg450.mat 450
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg500.mat 500

runB:
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg1000.mat 1000
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg1000.mat 1000
	
clean:
	-rm -Rf *.o *.bc $(PROG).exe

PROG=hpcg
PTHREAD_DIR=/usr/lib/x86_64-linux-gnu/
CFLAGS=-O0 -g -Werror -I. -I${OCR_INSTALL}/../../inc -c
OCR_FLAGS=-L${OCR_INSTALL}/lib -L${PTHREAD_DIR} -I${OCR_INSTALL}/include -locr -lpthread -static

CLANG=clang++
LLC=llc

ifndef OCR_INSTALL
$(error OCR_INSTALL not set)
endif

ifndef OCR_CONFIG
OCR_CONFIG=${OCR_INSTALL}/config/default.cfg
$(warning OCR_CONFIG not set, defaulting to ${OCR_CONFIG})
endif

OCR_RUN_FLAGS=-ocr:cfg ${OCR_CONFIG}

all-test: compile_pp run

compile:
	$(CLANG) -emit-llvm $(CFLAGS) timer.cpp main.cpp matrix.cpp
	$(CLANG) -g -o $(PROG).exe timer.bc main.bc matrix.bc $(OCR_FLAGS)

compile_pp:
	g++ $(CFLAGS) timer.cpp main.cpp matrix.cpp
	g++ -g -o $(PROG).exe timer.o main.o matrix.o $(OCR_FLAGS)

run:
#	octave ./cgmatrix.m
	./$(PROG).exe $(OCR_RUN_FLAGS) 0 cg500.mat 500
	./$(PROG).exe $(OCR_RUN_FLAGS) 1 cg500.mat 500

clean:
	-rm -Rf *.o *.bc $(PROG).exe

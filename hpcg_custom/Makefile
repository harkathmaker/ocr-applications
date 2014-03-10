PROG=hpcg
PTHREAD_DIR=/usr/lib/x86_64-linux-gnu/
CFLAGS=-O0 -g -Werror -I. -I${OCR_INSTALL}/../../inc -c
OCR_FLAGS=-L${OCR_INSTALL}/lib -L${PTHREAD_DIR} -I${OCR_INSTALL}/include -locr -lpthread -static

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
	g++ $(CFLAGS) main.cpp matrix.cpp
	g++ -g -o $(PROG).exe main.o matrix.o $(OCR_FLAGS)

run:
	./$(PROG).exe $(OCR_RUN_FLAGS)

clean:
	-rm -Rf *.o $(PROG).exe

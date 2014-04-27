CFLAGS=-I${OCR_ROOT}/include -I${OCR_ROOT}/../../src/inc -O0
LDFLAGS=-L${OCR_ROOT}/lib -locr -lm
OTHER_SOURCES=OCRSparseMatrix.cpp OCRDenseMatrix.cpp OCRBaseMatrix.cpp

%.o : %.c
	echo "$(CC) $(CFLAGS) $(LDFLAGS) $@"
	$(CC) $(CFLAGS) $(LDFLAGS) $@

%.out: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

%.o : %.cpp
	echo "$(CC) $(CFLAGS) $(LDFLAGS) $@"
	$(CC) $(CFLAGS) $(LDFLAGS) $@ 

%.out: %.cpp
	echo "our rule"
	$(CC) -lstdc++ $(CFLAGS) $(LDFLAGS) $< $(OTHER_SOURCES) -o $@

clean:
	rm *.out

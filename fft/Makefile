CFLAGS=-I${OCR_ROOT}/include -O3
LDFLAGS=-L${OCR_ROOT}/lib -locr -lm
CPC := g++
CC := gcc

ifdef mic
        CC := icc
        CFLAGS := $(CFLAGS) -mmic
        LDFLAG := $(LDFLAGS) -mmic -L/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/mic -limf -lsvml -lirng -lintlc
endif

%.o : %.c
	echo "$(CC) $(CFLAGS) $(LDFLAGS) $@"
	$(CC) $(CFLAGS) $(LDFLAGS) $@

%.out: %.c
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $<

%.o : %.cpp
	echo "$(CC) $(CFLAGS) $(LDFLAGS) $@"
	$(CPC) $(CFLAGS) $(LDFLAGS) $@

%.out: %.cpp
	echo "our rule"
	$(CPC) -lstdc++ $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm *.out

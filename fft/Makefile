CFLAGS=-I${OCR_ROOT}/include -O0
LDFLAGS=-L${OCR_ROOT}/lib -locr -lm

%.o : %.c
	echo "$(CC) $(CFLAGS) $(LDFLAGS) $@"
	$(CC) $(CFLAGS) $(LDFLAGS) $@

CFLAGS=-I${OCR_ROOT}/include -g
LDFLAGS=-L${OCR_ROOT}/lib -locr

%.o : %.c
	echo "$(CC) $(CFLAGS) $(LDFLAGS) $@"
	$(CC) $(CFLAGS) $(LDFLAGS) $@

CFLAGS=-I${OCR_ROOT}/include -O0
LDFLAGS=-L${OCR_ROOT}/lib -locr -lm

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
	$(CC) -lstdc++ $(CFLAGS) $(LDFLAGS) $< -o $@

clean:
	rm *.out

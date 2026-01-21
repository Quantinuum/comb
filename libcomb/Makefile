CFLAGS = -lpng -lc -lm # -lm4ri for dynamic linking, and LD_LIBRARY_PATH=/usr/local/lib at runtime
CC = gcc #g++?
M4RI = /usr/local/lib/libm4ri.a

demo: demo.c comb.h
	$(CC) demo.c -o demo $(M4RI) $(CFLAGS)
	./demo 1111001010110000

clean:
	rm demo

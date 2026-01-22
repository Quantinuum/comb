CFLAGS = -lpng -lc -lm # -lm4ri for dynamic linking, and LD_LIBRARY_PATH=/usr/local/lib at runtime
CC = gcc #g++?
M4RI = /usr/local/lib/libm4ri.a

demo: demo.c comb.c
	$(CC) demo.c comb.c -o demo $(M4RI) $(CFLAGS)
	./demo 1111001010110000

clean:
	rm -f demo
	rm -f *.o

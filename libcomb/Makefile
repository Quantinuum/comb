CFLAGS = -lpng -lc -lm # -lm4ri for dynamic linking, and LD_LIBRARY_PATH=/usr/local/lib at runtime
CC = gcc #g++?
M4RI = /usr/local/lib/libm4ri.a

c_files := $(wildcard *.c)
obj_files := $(c_files:.c=.o)

%.o : %.c
	$(CC) -fPIC -c $^ -o $@

libcomb.so: $(obj_files)
	gcc -O3 -shared -o libcomb.so $(CFLAGS) -rdynamic $(obj_files) $(M4RI)

clean:
	rm -f demo
	rm -f *.o

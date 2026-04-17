CFLAGS =  -lpng -lc -lm # -lm4ri for dynamic linking, and LD_LIBRARY_PATH=/usr/local/lib at runtime
CC = gcc #g++?
M4RI = /usr/local/lib/libm4ri.a

c_files := $(wildcard *.c)
obj_files := $(c_files:.c=.o)

%.o : %.c
	$(CC) -Wall -Werror -O3 -fPIC -c $^ -o $@

libcomb.so: $(obj_files)
	gcc -shared -o libcomb.so $(CFLAGS) -rdynamic $(obj_files) $(M4RI)

test_% : test/%.c $(obj_files) comb.h lwlog.h
	$(CC) -I. $< -fPIC -o $@ $(M4RI) $(CFLAGS) $(obj_files)

clean:
	rm -f demo
	rm -f *.o
	rm -f *.so
	rm -rf docs

docs: comb_external.h
	doxygen Doxyfile

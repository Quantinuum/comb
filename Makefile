# https on the safe side in case ssh not configured

M4RI_REPO = https://github.com/malb/m4ri.git
BRANCH = release-20240729
TMPDIR := $(shell mktemp -d)
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
M4RI_CONFIGURE_OPTIONS = --with-pic --prefix=$(ROOT_DIR)/m4ri --disable-png --disable-dependency-tracking

CFLAGS =  -fPIC -I./m4ri/include -lm # -lm4ri for dynamic linking, and LD_LIBRARY_PATH=/usr/local/lib at runtime
CC = gcc #g++?
M4RI = m4ri/lib/libm4ri.a

c_files := $(wildcard src/*.c)
obj_files := $(patsubst src/%.c,%.o,$(wildcard src/*.c))


libcomb.so: m4ri $(obj_files)
	$(CC) -shared -o libcomb.so $(CFLAGS) -rdynamic $(obj_files) $(M4RI)

m4ri:
	cd $(TMPDIR) && git clone $(M4RI_REPO) -b $(BRANCH) &&	cd m4ri/ && autoreconf --install && ./configure	$(M4RI_CONFIGURE_OPTIONS) && $(MAKE) && $(MAKE) check && $(MAKE) install

%.o : src/%.c m4ri
	$(CC) $(CFLAGS) -Wall -Werror -O3 -fPIC -c $< -o $@

m4ri-clean:
	rm -rf m4ri/

clean:
	rm -f *.o
	rm -f *.so
	rm -rf docs

docs:
	doxygen Doxyfile

# TODO fix these
test/% : test/%.c $(obj_files) comb.h lwlog.h
	$(CC) $(CFLAGS) -I. $< -fPIC -o $@ $(M4RI) $(CFLAGS) $(obj_files)

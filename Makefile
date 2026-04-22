# https on the safe side in case ssh not configured

M4RI_REPO = https://github.com/malb/m4ri.git
BRANCH = release-20240729
TMPDIR := $(shell mktemp -d)
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
M4RI_CONFIGURE_OPTIONS = --with-pic --prefix=$(ROOT_DIR)/m4ri --disable-png --disable-dependency-tracking

CFLAGS =  -fPIC -I./m4ri/include  -Wall -Werror -O3 -fPIC # -lm4ri for dynamic linking, and LD_LIBRARY_PATH=/usr/local/lib at runtime
CC = gcc #g++?
M4RI = m4ri/lib/libm4ri.a

c_files := $(wildcard src/*.c)
test_files :=$(wildcard test/*.c)
obj_files := $(patsubst src/%.c,%.o,$(wildcard src/*.c))
test_executables := $(patsubst test/%.c,test_%,$(wildcard test/*.c))


libcomb.so: m4ri $(obj_files)
	$(CC) -shared -o libcomb.so $(CFLAGS) -lm -rdynamic $(obj_files) $(M4RI)

m4ri:
	cd $(TMPDIR) && git clone $(M4RI_REPO) -b $(BRANCH) &&	cd m4ri/ && autoreconf --install && ./configure	$(M4RI_CONFIGURE_OPTIONS) && $(MAKE) && $(MAKE) check && $(MAKE) install

%.o : src/%.c m4ri
	$(CC) $(CFLAGS) -c $< -o $@

m4ri-clean:
	rm -rf m4ri/

clean:
	rm -f *.o
	rm -f *.so
	rm -rf docs
	rm -rf test-bin
	rm -f test_*

docs:
	doxygen Doxyfile

test-bin: $(test_executables)
	mkdir test-bin
	mv $(test_executables) test-bin

test_%: test/%.c libcomb.so
	$(CC) $(CFLAGS) -I./src $< -o $@ $(M4RI) $(obj_files) -lm

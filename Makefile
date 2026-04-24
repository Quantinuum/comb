# https on the safe side in case ssh not configured

M4RI_REPO = https://github.com/malb/m4ri.git
BRANCH = release-20240729
TMPDIR := $(shell mktemp -d)
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
M4RI_CONFIGURE_OPTIONS = --with-pic --prefix=$(ROOT_DIR)/m4ri --disable-png --disable-dependency-tracking

CFLAGS =  -fPIC -I./m4ri/include -Wall -Werror -O3 -fPIC
M4RI = m4ri/lib/libm4ri.a

c_files := $(wildcard src/c/*.c)
header_files := $(wildcard src/c/*.h)
test_files :=$(wildcard test/*.c)
obj_files := $(patsubst src/c/%.c,%.o,$(wildcard src/c/*.c))
test_executables := $(patsubst test/%.c,test_%,$(wildcard test/*.c))


libcomb.so: build
	ln -s build/libcomb.so libcomb.so

build: m4ri $(obj_files)
	mkdir -p build
	mv $(obj_files) build
	cd build && $(CC) -shared -o libcomb.so $(CFLAGS) -lm -rdynamic $(obj_files) ../$(M4RI)

m4ri:
	cd $(TMPDIR) && git clone $(M4RI_REPO) -b $(BRANCH) && cd m4ri/ && autoreconf --install && ./configure $(M4RI_CONFIGURE_OPTIONS) && $(MAKE) && $(MAKE) check && $(MAKE) install

%.o : src/c/%.c m4ri $(header_files)
	$(CC) $(CFLAGS) -c $< -o $@

m4ri-clean:
	rm -rf m4ri/

clean:
	rm -f *.o
	rm -f *.so
	rm -rf build
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

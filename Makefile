# https on the safe side in case ssh not configured

M4RI_REPO = https://github.com/malb/m4ri.git
BRANCH = release-20240729
TMPDIR := $(shell mktemp -d)
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
M4RI_CONFIGURE_OPTIONS = --with-pic --prefix=$(ROOT_DIR)/m4ri --disable-png --disable-dependency-tracking

override CFLAGS += -fPIC -I./m4ri/include -I./src/c -O3
M4RI = m4ri/lib/libm4ri.a

cpu_c_files := src/c/comb.c src/c/comb_cuda_flat.c src/c/comb_dense.c \
	src/c/comb_rank.c src/c/comb_reference.c src/c/convert.c src/c/dcp.c \
	src/c/gaussian_elimination.c src/c/init_free.c src/c/iteration.c \
	src/c/print.c src/c/utils.c
c_files := $(wildcard src/c/*.c)
header_files := $(wildcard src/c/*.h)
test_files :=$(wildcard test/*.c)
obj_files := $(patsubst src/c/%.c,build/%.o,$(cpu_c_files))
test_executables := $(patsubst test/%.c,test_%,$(wildcard test/*.c))

libcomb.so: $(obj_files)
	$(CC) -shared -o build/libcomb.so $(CFLAGS) -lm -rdynamic $(obj_files) $(M4RI)
	rm -f libcomb.so
	ln -s build/libcomb.so libcomb.so

m4ri:
	cd $(TMPDIR) && git clone $(M4RI_REPO) -b $(BRANCH) && cd m4ri/ && autoreconf --install && ./configure $(M4RI_CONFIGURE_OPTIONS) && $(MAKE) && $(MAKE) check && $(MAKE) install

build:
	mkdir -p build

build/%.o : src/c/%.c m4ri $(header_files) | build
	$(CC) $(CFLAGS) -c $< -o $@

m4ri-clean:
	rm -rf m4ri/

clean:
	rm -f *.so
	rm -rf build
	rm -rf test-bin
	rm -f test_*

docs: | build
	doxygen Doxyfile

cmake-configure:
	$(CMAKE) -S . -B build -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) $(CMAKE_OPTIONS)

cmake-build: cmake-configure
	$(CMAKE) --build build --parallel

cmake-test: cmake-build
	$(CTEST) --test-dir build --output-on-failure

cuda-configure:
	$(CMAKE) -S . -B build -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) -DCOMB_ENABLE_CUDA=ON $(CMAKE_OPTIONS)

mpi-configure:
	$(CMAKE) -S . -B build -DCMAKE_BUILD_TYPE=$(CMAKE_BUILD_TYPE) -DCOMB_ENABLE_CUDA=ON -DCOMB_ENABLE_MPI=ON $(CMAKE_OPTIONS)

test-bin: $(test_executables)
	mkdir test-bin
	mv $(test_executables) test-bin

test_test_comb_cuda_stub: test/test_comb_cuda_stub.c src/c/comb_cuda_stub.c \
		src/c/comb_cuda_flat.c
	$(CC) $(CFLAGS) $^ -o $@

test_test_comb_cuda_mpi_single_rank: test/test_comb_cuda_mpi_single_rank.c \
		src/mpi/comb_cuda_mpi.c src/c/comb_cuda_flat.c
	$(MPICC) $(CFLAGS) -I./src/mpi $^ -o $@

test_%: test/%.c libcomb.so
	$(CC) $(CFLAGS) -I./src $< -o $@ $(M4RI) $(obj_files) -lm

format:
	clang-format --style=Google --verbose -i $(c_files) $(header_files)

CMAKE ?= cmake
CTEST ?= ctest
MPICC ?= mpicc
CMAKE_BUILD_TYPE ?= Release
CMAKE_OPTIONS ?=

.PHONY: clean cmake-configure cmake-build cmake-test cuda-configure docs format \
	m4ri-clean mpi-configure

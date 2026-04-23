import cffi

ffi = cffi.FFI()

header_file=open("../c/comb_external.h","r")

# Strip lines beginning with # - preprocessor macros - from the header file as FFI can't deal with them
lines = [line for line in header_file.readlines() if not line[0] == '#']
ffi.cdef("".join(lines))

# Load the foreign library (ensure that 'make' is issued from the project root to build libcomb.so)
lib = ffi.dlopen("../../libcomb.so");



print(ffi)
print(lib)

print(lib.comb_alloc_bit_array(6))
print(lib.comb_verify_css_magic_state) 

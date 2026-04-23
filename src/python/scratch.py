# Python FFI to Comb using CFFI

import cffi
import subprocess
import sys
from pathlib import Path

library_path = "../../libcomb.so"
ffi = cffi.FFI()

# Strip lines beginning with # - preprocessor macros - from the header file as FFI can't deal with them
header_file = open("../c/comb_external.h","r")
ffi.cdef("".join([line for line in header_file.readlines() if not line[0] == '#']))

# If libcomb is not found, then try to build it

# subprocess invocation from https://trstringer.com/python-external-commands/ (under cc by 4.0 https://creativecommons.org/licenses/by/4.0/deed.en)
if not Path(library_path).is_file():
    cp = subprocess.run(
    "cd ../.. && make -j8",
    shell=True,
    text=True,
    # Send stderr to stdout and don't automatically throw an exception in order
    # to replicate checking the output of a shell command.
    stdout=subprocess.PIPE,
    stderr=subprocess.STDOUT,
    )
    print(cp.stdout)
    if cp.returncode != 0:
        print(f"Error! {cp.returncode}")
        # Maybe exit if you shouldn't continue on failure
        sys.exit(cp.returncode)

# Load the library
lib = ffi.dlopen(library_path);

# Test to see if it works

print(lib.comb_alloc_bit_array(6))
print(lib.comb_verify_css_magic_state) 

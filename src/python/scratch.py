# Python FFI to Comb using CFFI

import cffi

ffi = cffi.FFI()

header_file=open("../c/comb_external.h","r")

# Strip lines beginning with # - preprocessor macros - from the header file as FFI can't deal with them
lines = [line for line in header_file.readlines() if not line[0] == '#']
ffi.cdef("".join(lines))

import subprocess
import sys
from pathlib import Path

# If libcomb is not found, then try to build it

# from https://trstringer.com/python-external-commands/
if not Path("../../libcomb.so").is_file():
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
lib = ffi.dlopen("../../libcomb.so");

# Test to see if it works

print(lib.comb_alloc_bit_array(6))
print(lib.comb_verify_css_magic_state) 

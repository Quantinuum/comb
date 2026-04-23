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

# Example usage (inspired by goto)
n_all = 8
n = 7
n_x = 3
n_z = 3
k = 1
x_stabilisers = lib.comb_calloc_vector(n_x, n_all)
x_lookup_table = lib.comb_alloc_bit_array(2 ** n_x)
z_stabilisers = lib.comb_calloc_vector(n_z, n_all)
z_lookup_table = lib.comb_alloc_bit_array(2 ** n_z)
dcp_stabilisers = lib.comb_calloc_dcp(k, n_all)
fault_location_count = 10
fault_location_partitions = ffi.new("uint32_t[{}]".format(fault_location_count))
s = 1
bad_fault_array = ffi.new("uint32_t[{}]".format(s))
output_error_count = 100
output_errors = lib.comb_calloc_dcp(output_error_count, n_all)
x_measurement_count = 1
x_measurement_paulis = lib.comb_calloc_vector(x_measurement_count,n_all)
z_measurement_count = 1
z_measurement_paulis = lib.comb_calloc_vector(z_measurement_count,n_all)

for stabilisers in [x_stabilisers, z_stabilisers]:
    for i in range(n_all):
        for j in range(n_x):
            lib.comb_set_cvector_aref(stabilisers, j, i, 0)
    for i, j in [(0, [0, 2, 4, 6]), (1, [1, 2, 4 ,5]), (2, [3, 4, 5, 6])]:
        for index in j:
            lib.comb_set_cvector_aref(stabilisers, i, index, 1)

# Here, you would use the comb_set_* methods to set the other values of the foreign data

print(lib.comb_verify_css_magic_state(n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers,
                                      z_lookup_table, dcp_stabilisers, fault_location_count,
                                      fault_location_partitions, s, bad_fault_array, output_error_count,
                                      output_errors, x_measurement_count, x_measurement_paulis,
                                      z_measurement_count, z_measurement_paulis))

module LibComb

# Exposes comb_external.h. Refer to the documentation generated from the C sources there.
# Ensure that libcomb.so is built beforehand

# TODO use symbols explicitly and dlopen libcomb beforehand, not happy with ccall's (library, function)

export comb_alloc_bit_array
export comb_calloc_vector
export comb_calloc_dcp

export comb_set_cvector_aref
export comb_set_cdcp_F_aref
export comb_set_cdcp_h_x_aref
export comb_set_cdcp_h_z_aref
export comb_set_bit_array_aref

export comb_free_bit_array
export comb_cfree_vector
export comb_cfree_dcp

export comb_verify_css_magic_state
export comb_demo

library_path="./libcomb.so"

function comb_alloc_bit_array(n)
    return @ccall library_path.comb_alloc_bit_array(n::Cint)::Ptr{Cvoid}
end

function comb_calloc_vector(count, n)
    return @ccall library_path.comb_calloc_vector(count::Cint,n::Cint)::Ptr{Cvoid}
end

function comb_calloc_dcp(count, n)
    return @ccall library_path.comb_calloc_dcp(count::Cint,n::Cint)::Ptr{Cvoid}
end

function comb_set_cvector_aref(vector, n, i, value)
    return @ccall library_path.comb_set_cvector_aref(vector::Ptr{Cvoid},n::Cint,i::Cint,value::Cint)::Cvoid
end

function comb_set_cdcp_F_aref(dcp, n, i, j, value)
    return @ccall library_path.comb_set_cdcp_F_aref(dcp::Ptr{Cvoid},n::Cint,i::Cint,j::Cint,value::Cint)::Cvoid
end

function comb_set_cdcp_h_x_aref(dcp, n, i, value)
    return @ccall library_path.comb_set_cdcp_h_x_aref(dcp::Ptr{Cvoid},n::Cint,i::Cint,value::Cint)::Cvoid
end

function comb_set_cdcp_h_z_aref(dcp, n, i, value)
    return @ccall library_path.comb_set_cdcp_h_z_aref(dcp::Ptr{Cvoid},n::Cint,i::Cint,value::Cint)::Cvoid
end

function comb_set_bit_array_aref(array, i, value)
    return @ccall library_path.comb_set_bit_array_aref(array::Ptr{Cvoid},n::Cint,i::Cint,value::Cint)::Cvoid
end

function comb_free_bit_array(pointer)
    return @ccall library_path.comb_free_bit_array(pointer::Ptr{Cvoid})::Cvoid
end

function comb_cfree_vector(count, pointer)
    return @ccall library_path.comb_cfree_vector(count::Cint,pointer::Ptr{Cvoid})::Cvoid
end

function comb_cfree_dcp(count, pointer)
    return @ccall library_path.comb_cfree_dcp(count::Cint,pointer::Ptr{Cvoid})::Cvoid
end

function comb_verify_css_magic_state(n_all,n,n_x,n_z,k,x_stabilisers,x_lookup_table,z_stabilisers,z_lookup_table,dcp_stabilisers,fault_location_count,fault_location_partitions,s,bad_fault_array,output_error_count,output_errors,x_measurement_count,x_measurement_paulis,z_measurement_count,z_measurement_paulis,thread_count)
    return @ccall library_path.comb_verify_css_magic_state(n_all::UInt32,n::UInt32,n_x::UInt32,n_z::UInt32,k::UInt32,x_stabilisers::Ptr{Cvoid},x_lookup_table::Ptr{Cvoid},z_stabilisers::Ptr{Cvoid},z_lookup_table::Ptr{Cvoid},dcp_stabilisers::Ptr{Cvoid},fault_location_count::UInt32,fault_location_partitions::Ptr{Cvoid},s::UInt32,bad_fault_array::Ptr{Cvoid},output_error_count::UInt32,output_errors::UInt32,x_measurement_count::UInt32,x_measurement_paulis::Ptr{Cvoid},z_measurement_count::UInt32,z_measurement_paulis::Ptr{Cvoid},thread_count::UInt32)::Cint
end

# Friendlier interface to comb_verify_css_magic_state
function verify_css_magic_state(n_all,n,n_x,n_z,k,x_stabilisers,x_lookup_table,z_stabilisers,z_lookup_table,dcp_stabilisers,fault_location_count,fault_location_partitions,s,output_error_count,output_errors,x_measurement_count,x_measurement_paulis,z_measurement_count,z_measurement_paulis,thread_count)
    bad_fault_array = Vector{UInt32}(undef, s)
    code=comb_verify_css_magic_state(n_all,n,n_x,n_z,k,x_stabilisers,x_lookup_table,z_stabilisers,z_lookup_table,dcp_stabilisers,fault_location_count,fault_location_partitions,s,bad_fault_array,output_error_count,output_errors,x_measurement_count,x_measurement_paulis,z_measurement_count,z_measurement_paulis,thread_count)
    if code == 0
        return :all_faults_good
    end
    if code == 1
        return :bad_fault_found, bad_fault_array
    end
    if code == 2
        return :invalid_input
    end
    if code == 3
        return :problem_too_large
    end
    if code == 4
        return :internal_error
    end
    if code == 5
        return :other_error
    end
end

function comb_demo()
    n_all = 8
    n = 7
    n_x = 3
    n_z = 3
    k = 1
    x_stabilisers = comb_calloc_vector(n_x, n_all)
    x_lookup_table = comb_alloc_bit_array(2 ^ (n- n_x))
    z_stabilisers = comb_calloc_vector(n_z, n_all)
    z_lookup_table = comb_alloc_bit_array(2 ^ (n - n_z))
    dcp_stabilisers = comb_calloc_dcp(k, n_all)
    fault_location_count = 10
    fault_location_partitions = Vector{UInt32}(undef, fault_location_count);
    s = 1
    bad_fault_array = Vector{UInt32}(undef, s)
    output_error_count = 100
    output_errors = comb_calloc_dcp(output_error_count, n_all)
    x_measurement_count = 1
    x_measurement_paulis = comb_calloc_vector(x_measurement_count,n_all)
    z_measurement_count = 1
    z_measurement_paulis = comb_calloc_vector(z_measurement_count,n_all)
    thread_count = 1
    return verify_css_magic_state(n_all, n, n_x, n_z, k, x_stabilisers, x_lookup_table, z_stabilisers,
                                       z_lookup_table, dcp_stabilisers, fault_location_count,
                                       fault_location_partitions, s, output_error_count,
                                       output_errors, x_measurement_count, x_measurement_paulis,
                                       z_measurement_count, z_measurement_paulis, thread_count)
end
end

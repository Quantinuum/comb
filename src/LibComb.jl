module LibComb

# Write your package code here.

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

function comb_alloc_bit_array(n)
    return @ccall "./libcomb.so".comb_alloc_bit_array(n::Cint)::Ptr{Cvoid}
end

function comb_calloc_vector(count, n)
    return @ccall "./libcomb.so".comb_calloc_vector(count::Cint,n::Cint)::Ptr{Cvoid}
end

function comb_calloc_dcp(count, n)
    return @ccall "./libcomb.so".comb_calloc_vector(count::Cint,n::Cint)::Ptr{Cvoid}
end

function comb_set_cvector_aref(vector, n, i, value)
    return @ccall "./libcomb.so".comb_set_cvector_aref()::Cvoid
end

function comb_set_cdcp_F_aref(dcp, n, i, j, value)
    return @ccall "./libcomb.so".comb_set_cdcp_F_aref()::Cvoid
end

function comb_set_cdcp_h_x_aref(dcp, n, i, value)
    return @ccall "./libcomb.so".comb_set_cdcp_h_x_aref()::Cvoid
end

function comb_set_cdcp_h_z_aref(dcp, n, i, value)
    return @ccall "./libcomb.so".comb_set_cdcp_h_z_aref()::Cvoid
end

function comb_set_bit_array_aref(array, i, value)
    return @ccall "./libcomb.so".comb_set_bit_array_aref()::Cvoid
end

function comb_free_bit_array(pointer)
    return @ccall "./libcomb.so".comb_free_bit_array()::Cvoid
end

function comb_cfree_vector(count, pointer)
    return @ccall "./libcomb.so".comb_cfree_vector()::Cvoid
end

function comb_cfree_dcp(count, pointer)
    return @ccall "./libcomb.so".comb_cfree_dcp()::Cvoid
end

function comb_verify_css_magic_state(n_all,n,n_x,n_z,k,x_stabilisers,x_lookup_table,z_stabilisers,z_lookup_table,dcp_stabilisers,fault_location_count,fault_location_partitions,s,bad_fault_array,output_error_count,output_errors,x_measurement_count,x_measurement_paulis,z_measurement_count,z_measurement_paulis,thread_count)
    return @ccall "./libcomb.so".comb_verify_css_magic_state()::Cint
end

end

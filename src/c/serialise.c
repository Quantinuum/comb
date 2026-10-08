#include "comb.h"
#include <fcntl.h>
#include <unistd.h>

void store_uint32_t(int fd, uint32_t data){
  ssize_t count = write(fd, (void*)&data, 4);
  if(count!=4)exit(INTERNAL_ERROR);
}
uint32_t read_uint32_t(int fd){
  uint32_t buf;
  ssize_t count = read(fd, (void*)&buf, 4);
  if(count!=4)exit(INTERNAL_ERROR);
  return buf;
}


void store_bit_array(int fd, comb_bit_array_t table, uint32_t n){
  uint32_t* as_int = (uint32_t*)table;
  for(int i = 0; i < (1+ (n/4));i++){
    store_uint32_t(fd,as_int[i]);
  }
}
void read_bit_array(int fd, comb_bit_array_t table, uint32_t n){
  uint32_t* as_int = (uint32_t*)table;
  for(int i = 0; i < (1+ (n/4));i++){
    as_int[i]=read_uint32_t(fd);
  }
}

void store_vector(int fd, comb_vector_t* vector){
  store_uint32_t(fd,vector->n);
  store_bit_array(fd,vector->data,vector->n);
}
void read_vector(int fd, comb_vector_t *vector){
  vector->n = read_uint32_t(fd);
  read_bit_array(fd,vector->data,vector->n);
}
void store_matrix(int fd, comb_matrix_t* vector){
  store_uint32_t(fd,vector->m);
  store_uint32_t(fd,vector->n);
  store_bit_array(fd,vector->data,vector->m * vector->n);
}
void read_matrix(int fd, comb_matrix_t *vector){
  vector->m = read_uint32_t(fd);
  vector->n = read_uint32_t(fd);
  read_bit_array(fd,vector->data,vector->m * vector->n);
}
void store_dcp(int fd, comb_dcp_t *dcp){
  store_matrix(fd,&(dcp->F));
  store_vector(fd,&(dcp->h_x));
  store_vector(fd,&(dcp->h_z));
}
void read_dcp(int fd, comb_dcp_t *dcp){
  read_matrix(fd,&(dcp->F));
  read_vector(fd,&(dcp->h_x));
  read_vector(fd,&(dcp->h_z));
}

int serialise_run_to_file(int fd,    uint32_t n_all,  ///< number of qubits in code + all ancillas
    uint32_t n,      // number of qubits in code
    uint32_t n_x,    // number of independent x stabiliser generators in code
    uint32_t n_z,    // number of independent z stabiliser generators in code
    uint32_t k,  // number of logical qubits in code. equivalently, number of
                 // logical dcp operators that stabilise the code
    comb_pauli_x_t *_x_stabilisers,   // n-x of these
    comb_bit_array_t x_lookup_table,  // x lookup table. of size 2^(n-n_x)
    comb_pauli_z_t *_z_stabilisers, comb_bit_array_t z_lookup_table,
    comb_dcp_t *_dcp_stabilisers,   // k of these
    uint32_t fault_location_count,  // number of fault locations
    uint32_t
        *fault_location_partitions,  // at least fault_location_count of these.
                                     // number of faults at each location
    uint32_t s,  // number of faults to which circuit should be resilient to
    uint32_t output_error_count,   // number of output errors that result from
                                   // propagating all the faults. should be at
                                   // least sum(fault_location_partitions)
    comb_dcp_t *_output_errors,    // at least output_error_count of these
    uint32_t x_measurement_count,  // number of x measurement paulis
    comb_pauli_x_t
        *_x_measurement_paulis,    // postselections made in x basis. at
                                   // least x_measurement_count of these
    uint32_t z_measurement_count,  // number of z measurement paulis
                          comb_pauli_z_t *_z_measurement_paulis){
  store_uint32_t(fd, n_all);
  store_uint32_t(fd, n);
  store_uint32_t(fd, n_x);
  store_uint32_t(fd, n_z);
  store_uint32_t(fd, k);
  store_uint32_t(fd, fault_location_count);
  store_uint32_t(fd, s);
  store_uint32_t(fd, output_error_count);
  store_uint32_t(fd, x_measurement_count);
  store_uint32_t(fd, z_measurement_count);
  for(int i = 0; i < n_x; i++){
    store_vector(fd,_x_stabilisers+i);
  }
  store_bit_array(fd, x_lookup_table, 1 << (n - n_x));
  for(int i = 0; i < n_z; i++){
    store_vector(fd,_z_stabilisers+i);
  }
  store_bit_array(fd, z_lookup_table, 1 << (n - n_z));
  for(int i = 0; i < k; i++){
    store_dcp(fd,_dcp_stabilisers+i);
  }
  for(int i = 0; i < fault_location_count; i++){
    store_uint32_t(fd,fault_location_partitions[i]);
  }
  for(int i = 0; i < output_error_count; i++){
    store_dcp(fd,_output_errors+i);
  }
  for(int i = 0; i < x_measurement_count; i++){
    store_vector(fd,_x_measurement_paulis+i);
  }
  for(int i = 0; i < z_measurement_count; i++){
    store_vector(fd,_z_measurement_paulis+i);
  }
  return close(fd); // success
}

int deserialise_run_from_file(char *input_file, int thread_count){
  int fd = open(input_file, O_RDONLY);
  uint32_t n_all = read_uint32_t(fd);
  uint32_t n = read_uint32_t(fd);
  uint32_t n_x = read_uint32_t(fd);
  uint32_t n_z = read_uint32_t(fd);
  uint32_t k = read_uint32_t(fd);
  uint32_t fault_location_count = read_uint32_t(fd);
  uint32_t s = read_uint32_t(fd);
  uint32_t output_error_count = read_uint32_t(fd);
  uint32_t x_measurement_count = read_uint32_t(fd);
  uint32_t z_measurement_count = read_uint32_t(fd);
  uint32_t *bad_fault_array = (uint32_t*)calloc(s,sizeof(uint32_t));

  comb_pauli_x_t *_x_stabilisers = comb_calloc_vector(n_x,n);
  for(int i = 0; i < n_x; i++){
    read_vector(fd,_x_stabilisers+i);
  }
  comb_bit_array_t x_lookup_table = comb_alloc_bit_array(1 << (n - n_x));
  read_bit_array(fd, x_lookup_table, (1 << (n - n_x)));
  comb_pauli_z_t *_z_stabilisers = comb_calloc_vector(n_z,n);
    for(int i = 0; i < n_z; i++){
    read_vector(fd,_z_stabilisers+i);
  }
  comb_bit_array_t z_lookup_table = comb_alloc_bit_array(1 << (n - n_z));
  read_bit_array(fd, z_lookup_table, (1 << (n - n_z)));
  comb_dcp_t *_dcp_stabilisers = comb_calloc_dcp(k,n_all);
  for(int i = 0; i < k; i++){
    read_dcp(fd,_dcp_stabilisers+i);
  }
  uint32_t* fault_location_partitions=calloc(fault_location_count,sizeof(uint32_t));
  for(int i = 0; i < fault_location_count; i++){
    fault_location_partitions[i]=read_uint32_t(fd);
  }
  comb_dcp_t *_output_errors = comb_calloc_dcp(output_error_count,n_all);
  for(int i = 0; i < output_error_count; i++){
    read_dcp(fd,_output_errors+i);
  }
  comb_pauli_x_t *_x_measurement_paulis = comb_calloc_vector(x_measurement_count,n_all);
  for(int i = 0; i < x_measurement_count; i++){
    read_vector(fd,_x_measurement_paulis+i);
  }
  comb_pauli_z_t *_z_measurement_paulis = comb_calloc_vector(z_measurement_count,n_all);
  for(int i = 0; i < z_measurement_count; i++){
    read_vector(fd,_z_measurement_paulis+i);
  }
  if(close(fd))exit(INTERNAL_ERROR);
  // data read; now do the check
  comb_code_t code = comb_verify_css_magic_state(n_all,n,n_x,n_z,k,_x_stabilisers,x_lookup_table,_z_stabilisers,
                                                 z_lookup_table,_dcp_stabilisers,fault_location_count,
                                                 fault_location_partitions,s,bad_fault_array,output_error_count,
                                                 _output_errors,x_measurement_count,_x_measurement_paulis,
                                                 z_measurement_count,_z_measurement_paulis,thread_count);
  printf("code = %d\n",code);
  if(code == BAD_FAULT_FOUND)for(int i = 0; i < s; i++){
    printf("bad_fault_array[%d]=%d\n",i,bad_fault_array[i]);
  }
  return code; // success
}

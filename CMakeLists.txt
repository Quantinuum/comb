cmake_minimum_required(VERSION 3.24)

project(comb VERSION 0.1.0 LANGUAGES C)

include(CheckLanguage)
include(CTest)
include(GNUInstallDirs)

option(COMB_BUILD_CPU "Build the CPU verifier" ON)
option(COMB_ENABLE_CUDA "Build the NVIDIA CUDA verifier" OFF)
option(COMB_ENABLE_MPI "Build MPI orchestration targets" OFF)
option(COMB_BUILD_LEGACY_TESTS
       "Build the pre-CTest diagnostic programs in test/" OFF)

set(COMB_M4RI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/m4ri" CACHE PATH
    "M4RI installation prefix (the legacy Makefile installs here)")

set(COMB_HAVE_CPU OFF)

if(COMB_BUILD_CPU)
  find_path(
    COMB_M4RI_INCLUDE_DIR
    NAMES m4ri/m4ri.h
    HINTS "${COMB_M4RI_ROOT}/include")
  find_library(
    COMB_M4RI_LIBRARY
    NAMES m4ri
    HINTS "${COMB_M4RI_ROOT}/lib" "${COMB_M4RI_ROOT}/lib64")

  if(COMB_M4RI_INCLUDE_DIR AND COMB_M4RI_LIBRARY)
    add_library(M4RI::m4ri UNKNOWN IMPORTED)
    set_target_properties(
      M4RI::m4ri
      PROPERTIES IMPORTED_LOCATION "${COMB_M4RI_LIBRARY}"
                 INTERFACE_INCLUDE_DIRECTORIES "${COMB_M4RI_INCLUDE_DIR}")

    set(COMB_CPU_SOURCES
        src/c/comb.c
        src/c/comb_cuda_flat.c
        src/c/comb_dense.c
        src/c/comb_rank.c
        src/c/comb_reference.c
        src/c/convert.c
        src/c/dcp.c
        src/c/gaussian_elimination.c
        src/c/init_free.c
        src/c/iteration.c
        src/c/print.c
        src/c/utils.c)

    add_library(comb SHARED ${COMB_CPU_SOURCES})
    add_library(Comb::comb ALIAS comb)
    target_compile_features(comb PUBLIC c_std_11)
    target_include_directories(
      comb
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/c>
             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/comb>)
    target_link_libraries(comb PRIVATE M4RI::m4ri m)
    set_target_properties(comb PROPERTIES C_VISIBILITY_PRESET default)
    set(COMB_HAVE_CPU ON)

    install(TARGETS comb LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
  else()
    message(
      FATAL_ERROR
        "COMB_BUILD_CPU=ON, but M4RI was not found. Run `make m4ri`, install "
        "M4RI system-wide, set COMB_M4RI_ROOT, or explicitly configure the "
        "dependency-free contract tests with -DCOMB_BUILD_CPU=OFF.")
  endif()
endif()

# comb_rank is intentionally independent of M4RI. This makes its checked
# weighted rank/unrank contract testable on development hosts without the full
# scientific dependency stack.
if(BUILD_TESTING)
  add_executable(comb_rank_test src/c/comb_rank.c test/test_comb_rank.c)
  target_compile_features(comb_rank_test PRIVATE c_std_11)
  target_include_directories(comb_rank_test PRIVATE src/c)
  add_test(NAME comb.rank COMMAND comb_rank_test)

  add_executable(comb_dense_test src/c/comb_dense.c test/test_comb_dense.c)
  target_compile_features(comb_dense_test PRIVATE c_std_11)
  target_include_directories(comb_dense_test PRIVATE src/c)
  add_test(NAME comb.dense COMMAND comb_dense_test)

  add_executable(comb_cuda_flat_test src/c/comb_cuda_flat.c
                                    test/test_comb_cuda_flat.c)
  target_compile_features(comb_cuda_flat_test PRIVATE c_std_11)
  target_include_directories(comb_cuda_flat_test PRIVATE src/c)
  add_test(NAME comb.cuda_flat COMMAND comb_cuda_flat_test)

  add_executable(comb_cuda_stub_test src/c/comb_cuda_stub.c
                                    src/c/comb_cuda_flat.c
                                    test/test_comb_cuda_stub.c)
  target_compile_features(comb_cuda_stub_test PRIVATE c_std_11)
  target_include_directories(comb_cuda_stub_test PRIVATE src/c)
  add_test(NAME comb.cuda_stub COMMAND comb_cuda_stub_test)

  add_executable(
    comb_reference_test
    src/c/comb_reference.c
    src/c/comb_cuda_flat.c
    src/c/comb_dense.c
    src/c/comb_rank.c
    test/test_comb_reference.c)
  target_compile_features(comb_reference_test PRIVATE c_std_11)
  target_include_directories(comb_reference_test PRIVATE src/c)
  add_test(NAME comb.reference COMMAND comb_reference_test)

  # Exercise the production MPI scheduler without requiring CUDA or M4RI.  This
  # is deliberately a direct single-process CTest, not an mpiexec invocation:
  # multi-rank correctness remains a target-system validation gate.
  find_package(MPI QUIET COMPONENTS C)
  if(MPI_C_FOUND)
    add_executable(
      comb_cuda_mpi_single_rank_test
      src/mpi/comb_cuda_mpi.c
      src/c/comb_cuda_flat.c
      test/test_comb_cuda_mpi_single_rank.c)
    target_compile_features(comb_cuda_mpi_single_rank_test PRIVATE c_std_11)
    target_include_directories(comb_cuda_mpi_single_rank_test
                               PRIVATE src/c src/mpi)
    target_link_libraries(comb_cuda_mpi_single_rank_test PRIVATE MPI::MPI_C)
    if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
      target_compile_options(
        comb_cuda_mpi_single_rank_test
        PRIVATE -Wall -Wextra -Wpedantic -Werror -Wconversion
                -Wsign-conversion)
    endif()
    add_test(NAME comb.cuda_mpi_single_rank
             COMMAND comb_cuda_mpi_single_rank_test)
    set_tests_properties(
      comb.cuda_mpi_single_rank
      PROPERTIES LABELS "mpi;single-rank;mock"
                 ENVIRONMENT "CUDA_VISIBLE_DEVICES=comb-mock-device")
  endif()

  # This end-to-end equivalence gate is meaningful only when the actual M4RI
  # legacy verifier is present.  It compares the legacy decision/first witness
  # with the independent sequential FlatProblem oracle on generated valid CSS
  # instances; it is intentionally absent from dependency-free builds.
  if(COMB_HAVE_CPU)
    add_executable(comb_legacy_reference_test
                   test/test_comb_legacy_reference.c)
    target_compile_features(comb_legacy_reference_test PRIVATE c_std_11)
    target_include_directories(comb_legacy_reference_test PRIVATE src/c)
    target_link_libraries(comb_legacy_reference_test PRIVATE Comb::comb)
    add_test(NAME comb.legacy_reference COMMAND comb_legacy_reference_test)
    set_tests_properties(comb.legacy_reference
                         PROPERTIES LABELS "legacy;differential")
  endif()
endif()

if(COMB_BUILD_LEGACY_TESTS)
  if(NOT COMB_HAVE_CPU)
    message(FATAL_ERROR "COMB_BUILD_LEGACY_TESTS requires the CPU comb target")
  endif()
  foreach(test_source IN ITEMS
          test/dcp.c
          test/dcp_mul.c
          test/dcp_postselect_coset.c
          test/dcp_to_coset.c
          test/dcp_to_state.c)
    get_filename_component(test_name "${test_source}" NAME_WE)
    add_executable("comb_legacy_${test_name}" "${test_source}")
    target_include_directories("comb_legacy_${test_name}" PRIVATE src/c)
    target_link_libraries("comb_legacy_${test_name}" PRIVATE comb M4RI::m4ri m)
  endforeach()
endif()

if(COMB_ENABLE_CUDA)
  check_language(CUDA)
  if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "COMB_ENABLE_CUDA=ON, but no CUDA compiler was found")
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)

  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  set(COMB_CUDA_SOURCES src/c/comb_cuda.cu)
  if(NOT COMB_HAVE_CPU)
    list(APPEND COMB_CUDA_SOURCES src/c/comb_cuda_flat.c)
  endif()
  add_library(comb_cuda SHARED ${COMB_CUDA_SOURCES})
  add_library(Comb::cuda ALIAS comb_cuda)
  target_include_directories(
    comb_cuda
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/c>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/comb>)
  target_compile_options(comb_cuda PRIVATE
                         $<$<COMPILE_LANGUAGE:CUDA>:-lineinfo>)
  target_link_libraries(comb_cuda PRIVATE CUDA::cudart)
  if(COMB_HAVE_CPU)
    # The CPU target already owns the one copy of FlatProblem packing and
    # validation.  Reuse it instead of exporting duplicate flat symbols from
    # both shared libraries.
    target_link_libraries(comb_cuda PRIVATE Comb::comb)
  endif()
  set_target_properties(
    comb_cuda
    PROPERTIES CUDA_ARCHITECTURES "80-real;80-virtual"
               CUDA_SEPARABLE_COMPILATION ON)

  install(TARGETS comb_cuda LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
else()
  set(COMB_CUDA_STUB_SOURCES src/c/comb_cuda_stub.c)
  if(NOT COMB_HAVE_CPU)
    list(APPEND COMB_CUDA_STUB_SOURCES src/c/comb_cuda_flat.c)
  endif()
  add_library(comb_cuda_stub STATIC ${COMB_CUDA_STUB_SOURCES})
  add_library(Comb::cuda_stub ALIAS comb_cuda_stub)
  target_compile_features(comb_cuda_stub PUBLIC c_std_11)
  target_include_directories(
    comb_cuda_stub
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/c>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/comb>)
  if(COMB_HAVE_CPU)
    target_link_libraries(comb_cuda_stub PRIVATE Comb::comb)
  endif()
  install(TARGETS comb_cuda_stub ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()

if(COMB_ENABLE_MPI)
  if(NOT COMB_ENABLE_CUDA)
    message(FATAL_ERROR "COMB_ENABLE_MPI requires COMB_ENABLE_CUDA=ON")
  endif()
  find_package(MPI REQUIRED COMPONENTS C)

  add_library(comb_cuda_mpi STATIC src/mpi/comb_cuda_mpi.c)
  add_library(Comb::cuda_mpi ALIAS comb_cuda_mpi)
  target_compile_features(comb_cuda_mpi PRIVATE c_std_11)
  target_include_directories(
    comb_cuda_mpi
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/mpi>
           $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/comb>
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src/c)
  target_link_libraries(comb_cuda_mpi PUBLIC MPI::MPI_C Comb::cuda)
  if(COMB_HAVE_CPU)
    # comb owns FlatProblem packing in this configuration, and the MPI object
    # calls that internal broadcast interface directly.
    target_link_libraries(comb_cuda_mpi PUBLIC Comb::comb)
  endif()
  install(TARGETS comb_cuda_mpi ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
  install(FILES src/mpi/comb_cuda_mpi.h
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/comb)
endif()

install(FILES src/c/comb_external.h src/c/comb_cuda.h
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/comb)
if(COMB_HAVE_CPU)
  install(FILES src/c/comb_reference.h
          DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/comb)
endif()

find_package(Doxygen QUIET)
if(Doxygen_FOUND)
  add_custom_target(
    docs
    COMMAND "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile"
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Generating API documentation under build/docs"
    VERBATIM)
endif()

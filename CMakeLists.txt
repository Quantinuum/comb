cmake_minimum_required(VERSION 3.24)

project(comb VERSION 0.1.0 LANGUAGES C)

include(CheckLanguage)
include(CTest)
include(GNUInstallDirs)

option(COMB_BUILD_CPU "Build the CPU verifier" ON)
option(COMB_ENABLE_CUDA "Build the NVIDIA CUDA verifier" OFF)
option(COMB_ENABLE_MPI "Build MPI orchestration targets" OFF)
option(COMB_BUILD_CUDA_MPI_TESTS
       "Run native MPI/CUDA tests requiring one visible GPU per local rank" OFF)
set(COMB_CUDA_MPI_TEST_RANKS "1;2;3" CACHE STRING
    "MPI rank counts for opt-in native MPI/CUDA tests")
option(COMB_BUILD_DIAGNOSTIC_TESTS
       "Build the original diagnostic programs in test/" OFF)

if(COMB_BUILD_CUDA_MPI_TESTS AND
   (NOT COMB_ENABLE_CUDA OR NOT COMB_ENABLE_MPI OR NOT BUILD_TESTING))
  message(FATAL_ERROR
          "COMB_BUILD_CUDA_MPI_TESTS requires COMB_ENABLE_CUDA, COMB_ENABLE_MPI, "
          "and BUILD_TESTING")
endif()

set(COMB_M4RI_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/m4ri" CACHE PATH
    "M4RI installation prefix (the compatibility Makefile installs here)")

set(COMB_HAVE_CPU OFF)

# The serialized .run reader has no M4RI or CUDA dependency.  Keep it usable
# by both the GPU launcher and input-validation tools on development hosts.
add_library(comb_run_file STATIC src/mpi/comb_run_file.c)
target_compile_features(comb_run_file PUBLIC c_std_11)
target_include_directories(comb_run_file PUBLIC src/mpi src/c)

if(COMB_BUILD_CPU)
  find_package(OpenMP REQUIRED COMPONENTS C)

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
        src/c/chunk.c
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
        src/c/serialise.c
        src/c/utils.c)

    add_library(comb SHARED ${COMB_CPU_SOURCES})
    add_library(Comb::comb ALIAS comb)
    target_compile_features(comb PUBLIC c_std_11)
    target_include_directories(
      comb
      PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/c>
             $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/comb>)
    target_link_libraries(comb PRIVATE M4RI::m4ri OpenMP::OpenMP_C m)
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

  add_executable(comb_run_file_test src/c/comb_cuda_flat.c
                                  test/test_comb_run_file.c)
  target_link_libraries(comb_run_file_test PRIVATE comb_run_file)
  add_test(NAME comb.run_file COMMAND comb_run_file_test)

  add_executable(comb_cuda_cli_fixture_writer test/write_comb_cuda_cli_fixture.c)
  target_compile_features(comb_cuda_cli_fixture_writer PRIVATE c_std_11)

  # Execute the same scientific functions used by the CUDA kernel directly
  # on the host.  This checks their arithmetic without pretending to test the
  # CUDA compiler, runtime, scheduler, or GPU memory accesses.
  add_executable(
    comb_cuda_pipeline_test
    src/c/comb_cuda_preprocess.c
    src/c/comb_cuda_flat.c
    src/c/comb_dense.c
    src/c/comb_rank.c
    src/c/comb_reference.c
    test/test_comb_cuda_pipeline.c)
  target_compile_features(comb_cuda_pipeline_test PRIVATE c_std_11)
  target_include_directories(comb_cuda_pipeline_test PRIVATE src/c)
  add_test(NAME comb.cuda_pipeline_host COMMAND comb_cuda_pipeline_test)
  set_tests_properties(comb.cuda_pipeline_host
                       PROPERTIES LABELS "host;scientific;differential")

  add_executable(comb_cuda_math_test src/c/comb_dense.c
                                    test/test_comb_cuda_math.c)
  target_compile_features(comb_cuda_math_test PRIVATE c_std_11)
  target_include_directories(comb_cuda_math_test PRIVATE src/c)
  add_test(NAME comb.cuda_math_host COMMAND comb_cuda_math_test)
  set_tests_properties(comb.cuda_math_host
                       PROPERTIES LABELS "host;scientific;differential")

  # Exercise the production MPI scheduler without requiring CUDA or M4RI.
  # Single-rank and multi-rank host/mock checks run here; actual multi-GPU
  # execution remains a target-system validation gate.
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

    add_executable(
      comb_cuda_mpi_multi_rank_test
      src/mpi/comb_cuda_mpi.c
      src/c/comb_cuda_flat.c
      test/test_comb_cuda_mpi_multi_rank.c)
    target_compile_features(comb_cuda_mpi_multi_rank_test PRIVATE c_std_11)
    target_include_directories(comb_cuda_mpi_multi_rank_test PRIVATE src/c src/mpi)
    target_link_libraries(comb_cuda_mpi_multi_rank_test PRIVATE MPI::MPI_C)
    foreach(ranks IN ITEMS 2 3)
      add_test(NAME comb.cuda_mpi_${ranks}_ranks
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                       ${MPIEXEC_PREFLAGS}
                       $<TARGET_FILE:comb_cuda_mpi_multi_rank_test>
                       ${MPIEXEC_POSTFLAGS})
      set_tests_properties(comb.cuda_mpi_${ranks}_ranks
                           PROPERTIES LABELS "mpi;multi-rank;mock"
                                      PROCESSORS ${ranks} TIMEOUT 60)
    endforeach()

    # Test the production MPI transport with the real scientific candidate
    # functions executing on the host.  The test-only adapter supplies CUDA
    # API entry points; it is never linked into a production target.
    add_library(comb_cuda_host_test_backend STATIC
                test/comb_cuda_host_test_backend.c
                src/c/comb_cuda_preprocess.c src/c/comb_cuda_flat.c)
    target_compile_features(comb_cuda_host_test_backend PRIVATE c_std_11)
    target_include_directories(comb_cuda_host_test_backend PRIVATE src/c)
    target_link_libraries(comb_cuda_host_test_backend PRIVATE MPI::MPI_C)
    add_executable(comb_cuda_mpi_pipeline_test
                   src/mpi/comb_cuda_mpi.c src/c/comb_reference.c
                   src/c/comb_dense.c src/c/comb_rank.c
                   test/test_comb_cuda_mpi_pipeline.c)
    target_compile_features(comb_cuda_mpi_pipeline_test PRIVATE c_std_11)
    target_include_directories(comb_cuda_mpi_pipeline_test PRIVATE src/c src/mpi)
    target_link_libraries(comb_cuda_mpi_pipeline_test
                          PRIVATE comb_cuda_host_test_backend MPI::MPI_C)
    foreach(ranks IN ITEMS 1 2 3)
      add_test(NAME comb.cuda_mpi_pipeline_${ranks}_ranks
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                       ${MPIEXEC_PREFLAGS}
                       $<TARGET_FILE:comb_cuda_mpi_pipeline_test>
                       ${MPIEXEC_POSTFLAGS})
      set_tests_properties(comb.cuda_mpi_pipeline_${ranks}_ranks
                           PROPERTIES LABELS "mpi;host;scientific;differential"
                                      PROCESSORS ${ranks} TIMEOUT 60)
    endforeach()

    add_executable(comb_cuda_mpi_cli_test src/mpi/main.c src/mpi/comb_cuda_mpi.c)
    target_compile_features(comb_cuda_mpi_cli_test PRIVATE c_std_11)
    target_include_directories(comb_cuda_mpi_cli_test PRIVATE src/c src/mpi)
    target_link_libraries(comb_cuda_mpi_cli_test PRIVATE comb_run_file
                          comb_cuda_host_test_backend MPI::MPI_C)
    add_test(NAME comb.cuda_mpi_cli_help COMMAND comb_cuda_mpi_cli_test --help)
    set_tests_properties(comb.cuda_mpi_cli_help PROPERTIES
                          PASS_REGULAR_EXPRESSION "Usage: comb-cuda-mpi")
    add_test(NAME comb.cuda_mpi_cli_missing_file
             COMMAND ${CMAKE_COMMAND}
                     "-DPROGRAM=$<TARGET_FILE:comb_cuda_mpi_cli_test>"
                     "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
                     "-DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}"
                     "-DPREFLAGS=${MPIEXEC_PREFLAGS}"
                     "-DPOSTFLAGS=${MPIEXEC_POSTFLAGS}"
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_comb_cuda_cli.cmake)
    set_tests_properties(comb.cuda_mpi_cli_missing_file PROPERTIES
                          LABELS "mpi;cli" PROCESSORS 2 TIMEOUT 60)
    add_test(NAME comb.cuda_mpi_cli_science_host
             COMMAND ${CMAKE_COMMAND}
                     "-DPROGRAM=$<TARGET_FILE:comb_cuda_mpi_cli_test>"
                     "-DFIXTURE_WRITER=$<TARGET_FILE:comb_cuda_cli_fixture_writer>"
                     "-DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/testdata/cli_host"
                     "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
                     "-DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}" -DRANKS=2
                     "-DPREFLAGS=${MPIEXEC_PREFLAGS}"
                     "-DPOSTFLAGS=${MPIEXEC_POSTFLAGS}"
                     -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_comb_cuda_cli_science.cmake)
    set_tests_properties(comb.cuda_mpi_cli_science_host PROPERTIES
                          LABELS "mpi;cli;host;scientific" PROCESSORS 2 TIMEOUT 90)
  endif()

  # This end-to-end equivalence gate is meaningful only when the main M4RI CPU
  # verifier is present.  It compares the main CPU decision and witness with
  # the independent sequential FlatProblem reference on generated valid CSS
  # instances that deliberately fit in one main-CPU chunk.  The main verifier
  # does not define a minimum-rank witness across multiple chunks.
  if(COMB_HAVE_CPU)
    add_executable(comb_main_reference_test test/test_comb_main_reference.c)
    target_compile_features(comb_main_reference_test PRIVATE c_std_11)
    target_include_directories(comb_main_reference_test PRIVATE src/c)
    target_link_libraries(comb_main_reference_test PRIVATE Comb::comb)
    add_test(NAME comb.main_reference COMMAND comb_main_reference_test)
    set_tests_properties(comb.main_reference
                         PROPERTIES LABELS "main-cpu;differential")
  endif()
endif()

if(COMB_BUILD_DIAGNOSTIC_TESTS)
  if(NOT COMB_HAVE_CPU)
    message(FATAL_ERROR
            "COMB_BUILD_DIAGNOSTIC_TESTS requires the CPU comb target")
  endif()
  foreach(test_source IN ITEMS
          test/dcp.c
          test/dcp_mul.c
          test/dcp_postselect_coset.c
          test/dcp_to_coset.c
          test/dcp_to_state.c)
    get_filename_component(test_name "${test_source}" NAME_WE)
    add_executable("comb_cpu_${test_name}" "${test_source}")
    target_include_directories("comb_cpu_${test_name}" PRIVATE src/c)
    target_link_libraries("comb_cpu_${test_name}" PRIVATE comb M4RI::m4ri m)
  endforeach()
endif()

if(COMB_ENABLE_CUDA)
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "80-real;80-virtual" CACHE STRING
        "CUDA target architectures")
  endif()
  check_language(CUDA)
  if(NOT CMAKE_CUDA_COMPILER)
    message(FATAL_ERROR "COMB_ENABLE_CUDA=ON, but no CUDA compiler was found")
  endif()
  enable_language(CUDA)
  find_package(CUDAToolkit REQUIRED)

  set(CMAKE_CUDA_STANDARD 17)
  set(CMAKE_CUDA_STANDARD_REQUIRED ON)

  set(COMB_CUDA_SOURCES src/c/comb_cuda.cu src/c/comb_cuda_preprocess.c)
  if(NOT COMB_HAVE_CPU)
    list(APPEND COMB_CUDA_SOURCES src/c/comb_cuda_flat.c)
  endif()
  add_library(comb_cuda SHARED ${COMB_CUDA_SOURCES})
  target_compile_features(comb_cuda PRIVATE c_std_11)
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
  # All device functions are defined in this single CUDA translation unit or
  # its inline headers.  Avoid an unnecessary device-link step: retaining the
  # original fatbinary also preserves the configured virtual PTX image.
  set_target_properties(comb_cuda PROPERTIES CUDA_SEPARABLE_COMPILATION OFF)

  if(BUILD_TESTING)
    add_executable(
      comb_cuda_device_test
      src/c/comb_reference.c
      src/c/comb_dense.c
      src/c/comb_rank.c
      test/test_comb_cuda_pipeline.c)
    target_compile_features(comb_cuda_device_test PRIVATE c_std_11)
    target_compile_definitions(comb_cuda_device_test PRIVATE COMB_TEST_CUDA_DEVICE=1)
    target_include_directories(comb_cuda_device_test PRIVATE src/c)
    target_link_libraries(comb_cuda_device_test PRIVATE Comb::cuda)
    if(COMB_HAVE_CPU)
      target_link_libraries(comb_cuda_device_test PRIVATE Comb::comb)
    endif()
    add_test(NAME comb.cuda_device COMMAND comb_cuda_device_test)
    set_tests_properties(comb.cuda_device
                         PROPERTIES LABELS "cuda;gpu;scientific;differential"
                                    SKIP_RETURN_CODE 77 RUN_SERIAL TRUE)
  endif()

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
  add_executable(comb-cuda-mpi src/mpi/main.c)
  target_link_libraries(comb-cuda-mpi PRIVATE Comb::cuda_mpi comb_run_file)
  if(COMB_BUILD_CUDA_MPI_TESTS)
    # Opt in only on execution hosts with the requested GPU allocation.
    # Unlike login-node tests, these checks must fail if CUDA is unavailable.
    add_executable(comb_cuda_mpi_device_test test/test_comb_cuda_mpi_pipeline.c
                   src/c/comb_reference.c src/c/comb_dense.c src/c/comb_rank.c)
    target_compile_features(comb_cuda_mpi_device_test PRIVATE c_std_11)
    target_compile_definitions(comb_cuda_mpi_device_test
                               PRIVATE COMB_TEST_CUDA_DEVICE=1)
    target_include_directories(comb_cuda_mpi_device_test PRIVATE src/c)
    target_link_libraries(comb_cuda_mpi_device_test PRIVATE Comb::cuda_mpi)
    foreach(ranks IN LISTS COMB_CUDA_MPI_TEST_RANKS)
      if(NOT ranks MATCHES "^[1-9][0-9]*$")
        message(FATAL_ERROR "Invalid native GPU MPI test rank count: ${ranks}")
      endif()
      add_test(NAME comb.cuda_mpi_device_${ranks}_ranks
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${ranks}
                       ${MPIEXEC_PREFLAGS} $<TARGET_FILE:comb_cuda_mpi_device_test>
                       ${MPIEXEC_POSTFLAGS})
      set_tests_properties(comb.cuda_mpi_device_${ranks}_ranks PROPERTIES
                            LABELS "mpi;cuda;gpu;scientific;differential"
                            PROCESSORS ${ranks} RUN_SERIAL TRUE TIMEOUT 120)
      add_test(NAME comb.cuda_mpi_cli_device_${ranks}_ranks
               COMMAND ${CMAKE_COMMAND}
                       "-DPROGRAM=$<TARGET_FILE:comb-cuda-mpi>"
                       "-DFIXTURE_WRITER=$<TARGET_FILE:comb_cuda_cli_fixture_writer>"
                       "-DWORKDIR=${CMAKE_CURRENT_BINARY_DIR}/testdata/cli_gpu_${ranks}"
                       "-DMPIEXEC=${MPIEXEC_EXECUTABLE}"
                       "-DNUMPROC_FLAG=${MPIEXEC_NUMPROC_FLAG}" "-DRANKS=${ranks}"
                       "-DPREFLAGS=${MPIEXEC_PREFLAGS}"
                       "-DPOSTFLAGS=${MPIEXEC_POSTFLAGS}"
                       -P ${CMAKE_CURRENT_SOURCE_DIR}/test/check_comb_cuda_cli_science.cmake)
      set_tests_properties(comb.cuda_mpi_cli_device_${ranks}_ranks PROPERTIES
                            LABELS "mpi;cuda;gpu;cli;scientific"
                            PROCESSORS ${ranks} RUN_SERIAL TRUE TIMEOUT 120)
    endforeach()
  endif()
  install(TARGETS comb-cuda-mpi RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
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

# The launcher must broadcast root's input error, exit on every process, and
# preserve Comb's INVALID_INPUT=2 status instead of hanging in a collective.
execute_process(COMMAND "${MPIEXEC}" "${NUMPROC_FLAG}" 2 ${PREFLAGS}
                        "${PROGRAM}" ${POSTFLAGS}
                        "${CMAKE_CURRENT_BINARY_DIR}/does-not-exist-comb-input.run"
                RESULT_VARIABLE result OUTPUT_VARIABLE output
                ERROR_VARIABLE error TIMEOUT 30)
if(NOT result EQUAL 2)
  message(FATAL_ERROR "Expected INVALID_INPUT=2; got ${result}\n${output}\n${error}")
endif()
if(NOT error MATCHES "Cannot read")
  message(FATAL_ERROR "Missing root input diagnostic: ${output}\n${error}")
endif()

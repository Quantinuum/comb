file(MAKE_DIRECTORY "${WORKDIR}")
execute_process(COMMAND "${FIXTURE_WRITER}" "${WORKDIR}"
                RESULT_VARIABLE generated ERROR_VARIABLE generate_error)
if(NOT generated EQUAL 0)
  message(FATAL_ERROR "Cannot generate CLI fixtures: ${generate_error}")
endif()
set(good_checkpoint "${WORKDIR}/good.checkpoint")
set(bad_checkpoint "${WORKDIR}/bad.checkpoint")
file(REMOVE "${good_checkpoint}" "${bad_checkpoint}")

function(run_cli name input expected_code expected_checked)
  execute_process(COMMAND "${MPIEXEC}" "${NUMPROC_FLAG}" "${RANKS}"
                          ${PREFLAGS} "${PROGRAM}" ${POSTFLAGS} "${input}"
                          --epoch-size 2 --chunk-size 1 ${ARGN}
                  RESULT_VARIABLE result OUTPUT_VARIABLE output
                  ERROR_VARIABLE error TIMEOUT 25)
  if(NOT result EQUAL expected_code)
    message(FATAL_ERROR "${name}: expected exit ${expected_code}, got ${result}\n${output}\n${error}")
  endif()
  if(NOT output MATCHES "code = ${expected_code}")
    message(FATAL_ERROR "${name}: missing scientific result\n${output}\n${error}")
  endif()
  if(NOT output MATCHES "checked = ${expected_checked} / 4")
    message(FATAL_ERROR "${name}: wrong checked count\n${output}\n${error}")
  endif()
  if(expected_code EQUAL 1 AND
     (NOT output MATCHES "bad_rank = 3" OR
      NOT output MATCHES "bad_fault_array\\[0\\]=3"))
    message(FATAL_ERROR "${name}: wrong bad witness\n${output}\n${error}")
  endif()
endfunction()

run_cli(good "${WORKDIR}/good.run" 0 4 --checkpoint "${good_checkpoint}")
file(READ "${good_checkpoint}" good_state)
if(NOT good_state MATCHES "next_epoch 4")
  message(FATAL_ERROR "Wrong all-good checkpoint: ${good_state}")
endif()
run_cli(good_resume "${WORKDIR}/good.run" 0 0
        --checkpoint "${good_checkpoint}" --resume)
run_cli(bad "${WORKDIR}/bad.run" 1 4 --checkpoint "${bad_checkpoint}")
file(READ "${bad_checkpoint}" bad_state)
if(NOT bad_state MATCHES "next_epoch 2")
  message(FATAL_ERROR "Wrong bad-case checkpoint: ${bad_state}")
endif()
run_cli(bad_resume "${WORKDIR}/bad.run" 1 2
        --checkpoint "${bad_checkpoint}" --resume)
run_cli(bad_fast "${WORKDIR}/bad.run" 1 "[1-4]" --fast-any)

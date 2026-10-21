#!/usr/bin/env bash
set -euo pipefail

# Build and run the CPU OpenMP tests using an installed M4RI.
source_root=$(cd "$(dirname "$0")/.." && pwd)
m4ri_prefix=${M4RI_PREFIX:-"$source_root/m4ri"}
compiler=${CC:-cc}
build_dir=$(mktemp -d "${TMPDIR:-/tmp}/comb-cpu-openmp.XXXXXX")
flags=(-std=gnu11 -D_GNU_SOURCE -UNDEBUG -O2 -g -fPIC -I"$source_root/src/c"
       -I"$m4ri_prefix/include")
link_flags=(-L"$m4ri_prefix/lib" -Wl,-rpath,"$m4ri_prefix/lib" -lm4ri -lm)
sanitizer=0
if [[ $(uname -s) == Darwin ]]; then
  omp_prefix=${OPENMP_PREFIX:-/opt/homebrew/opt/libomp}
  flags+=(-Xclang -fopenmp -I"$omp_prefix/include")
  link_flags+=(-L"$omp_prefix/lib" -Wl,-rpath,"$omp_prefix/lib" -lomp)
  library=libcomb.dylib
else
  flags+=(-fopenmp)
  link_flags+=(-fopenmp)
  library=libcomb.so
fi
if [[ ${1:-} == --tsan ]]; then
  sanitizer=1
  archer=${ARCHER_LIBRARY:-}
  if [[ ! -f $archer ]]; then
    printf '%s\n' 'Set ARCHER_LIBRARY to the LLVM libarcher shared library.' \
      'OpenMP-aware ThreadSanitizer checks require its barrier annotations.' >&2
    exit 2
  fi
  export OMP_TOOL=enabled OMP_TOOL_LIBRARIES="$archer" ARCHER_OPTIONS=verbose=1
  export TSAN_OPTIONS=halt_on_error=1:ignore_noninstrumented_modules=1
  flags+=(-fsanitize=thread -fno-omit-frame-pointer)
  link_flags+=(-fsanitize=thread)
elif [[ $# != 0 ]]; then
  printf 'Usage: %s [--tsan]\n' "$0" >&2
  exit 2
fi

sources=()
for source in "$source_root"/src/c/*.c; do
  [[ $source == */main.c ]] || sources+=("$source")
done
"$compiler" "${flags[@]}" -shared "${sources[@]}" "${link_flags[@]}" \
  -o "$build_dir/$library"
"$compiler" "${flags[@]}" "$source_root/test/cpu_openmp.c" \
  -L"$build_dir" -Wl,-rpath,"$build_dir" -lcomb "${link_flags[@]}" \
  -o "$build_dir/cpu_openmp"

printf 'Build and test logs: %s\n' "$build_dir"
cd "$build_dir" # main serialises each invocation into the current directory.
for threads in ${COMB_TEST_THREADS:-1 2 4 8}; do
  for mode in good bad mixed format; do
    for ((repeat = 1; repeat <= ${COMB_TEST_REPEATS:-5}; ++repeat)); do
      log="$mode-$threads-$repeat.log"
      if ! env OMP_NUM_THREADS="$threads" OMP_DYNAMIC=FALSE \
          "$build_dir/cpu_openmp" "$mode" \
          "$threads" >"$log" 2>&1; then
        cat "$log" >&2
        exit 1
      fi
      if [[ $sanitizer == 1 ]] && ! grep -q \
          'Archer detected OpenMP application with TSan, supplying OpenMP synchronization semantics' "$log"; then
        printf 'Archer did not activate; inspect %s\n' "$build_dir/$log" >&2
        exit 1
      fi
      if [[ $mode == format ]]; then
        continue
      fi
      # Check the printed counts: the public API only returns the decision.
      awk '
        $2 == "chunks" { chunks = $1 }
        /^Verified chunk / { completed++; sum += $4 }
        /^Took / { count = $4; totals++ }
        /^openmp-test / {
          split($4, expected, "="); oracle = expected[2]; seen++
        }
        END {
          if (chunks != 11 || completed != 11 || totals != 1 || seen != 1 ||
              count != oracle || sum != oracle) {
            print "OpenMP count/chunk regression failed"; exit 1
          }
        }' "$log"
    done
    printf 'PASS mode=%s threads=%s repetitions=%s\n' \
      "$mode" "$threads" "${COMB_TEST_REPEATS:-5}"
  done
done

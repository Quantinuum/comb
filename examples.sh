#!/usr/bin/env bash
set -u
set -o pipefail
shopt -s nullglob

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
if (( $# > 1 )); then
  echo "Usage: $0 [PRIVATE_INPUT_DIRECTORY]" >&2
  exit 2
fi
examples_dir="${1:-${COMB_EXAMPLES_DIR:-$script_dir/examples}}"
good_inputs=("$examples_dir"/fault-tolerant/*.run "$examples_dir"/fault-tolerant/*.run.gz)
bad_inputs=("$examples_dir"/non-fault-tolerant/*.run "$examples_dir"/non-fault-tolerant/*.run.gz)
if (( ${#good_inputs[@]} + ${#bad_inputs[@]} == 0 )); then
  echo "No private verification inputs found in: $examples_dir" >&2
  echo "Supply a directory containing fault-tolerant/ and non-fault-tolerant/ .run or .run.gz files:" >&2
  echo "  COMB_EXAMPLES_DIR=/path/to/private/inputs $0" >&2
  echo "Verification circuit inputs are not distributed with this repository." >&2
  exit 2
fi

if ! make -C "$script_dir" comb; then
  exit 1
fi
temporary_dir="$(mktemp -d "${TMPDIR:-/tmp}/comb-examples.XXXXXX")" || exit 1
trap 'rm -rf -- "$temporary_dir"' EXIT
failures=0
checked=0

check_input() {
  local input="$1" expected="$2" code
  local temporary_run="$temporary_dir/input.run"
  if [[ "$input" == *.gz ]]; then
    if ! gzip -dc -- "$input" > "$temporary_run"; then
      echo "FAIL: cannot decompress $input" >&2
      failures=$((failures + 1))
      return
    fi
  elif ! cp -- "$input" "$temporary_run"; then
    echo "FAIL: cannot read $input" >&2
    failures=$((failures + 1))
    return
  fi
  if (cd -- "$temporary_dir" &&
      LD_LIBRARY_PATH="$script_dir/m4ri/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}" \
        "$script_dir/comb" "$temporary_run" 1); then
    code=0
  else
    code=$?
  fi
  checked=$((checked + 1))
  if (( code == expected )); then
    echo "PASS: $input (exit $code)"
  else
    echo "FAIL: $input (expected exit $expected, got $code)" >&2
    failures=$((failures + 1))
  fi
}

for input in "${good_inputs[@]}"; do check_input "$input" 0; done
for input in "${bad_inputs[@]}"; do check_input "$input" 1; done
echo "Checked $checked verification inputs; $failures failures."
(( failures == 0 ))

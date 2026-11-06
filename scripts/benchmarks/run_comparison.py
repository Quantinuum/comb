#!/usr/bin/env python3
"""Time complete MPI commands, with scientific-result gates and raw logs.

Requires only Python 3.8+. Example configuration (paths relative to its file):

{
  "output_directory": "results",
  "warmups": 1,
  "repetitions": 3,
  "timeout_seconds": 600,
  "profiles": [
    {"name": "cpu_2workers", "kind": "cpu_mpi", "expected_ranks": 3,
     "argv": ["/path/mpiexec", "-n", "3", "/path/comb", "{input}"],
     "environment": {"OMP_NUM_THREADS": "1"}},
    {"name": "gpu_2workers", "kind": "cuda_mpi",
     "argv": ["/path/mpiexec", "-n", "2", "/path/comb-cuda-mpi", "{input}"],
     "environment": {"OMP_NUM_THREADS": "1", "CUDA_VISIBLE_DEVICES": "1,2"}}
  ],
  "workloads": [
    {"name": "steane", "input": "steane.run", "expected_count": 944}
  ]
}

Optional profile fields: cwd, timeout_seconds. Workload timeout_seconds takes
precedence over profile and global values. keep_going defaults to false: a
failed warmup or measurement stops the comparison. metadata may hold arbitrary
JSON provenance. No command is evaluated through a shell. The sole argv
placeholder is {input}. Output timings include MPI startup, input parsing,
preprocessing, GPU initialization/transfers, verification, output I/O, cleanup,
and MPI shutdown. They are not GPU-kernel-only timings.

One warmup and three measured repetitions are the defaults. Profiles rotate
order between rounds; with two profiles the order swaps. All repetitions and
warmups must pass the scientific gates before a summary is marked qualified.
"""

import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import signal
import statistics
import subprocess
import sys
import threading
import time


CPU_RESULT = re.compile(
    r"Process\s+(\d+):\s+(\d+)\s+chunks checked,\s+Took\s+"
    r"([0-9]+(?:\.[0-9]+)?(?:[eE][-+]?[0-9]+)?)\s+seconds,\s+"
    r"(\d+)\s+fault combinations checked\."
)
GPU_RESULT = re.compile(r"^checked\s*=\s*(\d+)\s*/\s*(\d+)\s*$", re.MULTILINE)
GPU_CODE = re.compile(r"^code\s*=\s*(-?\d+)\s*$", re.MULTILINE)
SAFE_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


def positive_number(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(label + " must be a positive finite number")
    if value <= 0 or not math.isfinite(value):
        raise ValueError(label + " must be a positive finite number")
    return float(value)


def integer(value, label, minimum=0):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError("{} must be an integer >= {}".format(label, minimum))
    return value


def resolve_path(value, base):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (base / path).resolve()


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_config(config, base, args):
    warmups = config.get("warmups", 1) if args.warmups is None else args.warmups
    repetitions = (config.get("repetitions", 3) if args.repetitions is None
                   else args.repetitions)
    config["warmups"] = integer(warmups, "warmups")
    config["repetitions"] = integer(repetitions, "repetitions", 1)
    config["timeout_seconds"] = positive_number(
        config.get("timeout_seconds", 600), "timeout_seconds")
    config["termination_grace_seconds"] = positive_number(
        config.get("termination_grace_seconds", 5), "termination_grace_seconds")
    config["keep_going"] = config.get("keep_going", False)
    if not isinstance(config["keep_going"], bool):
        raise ValueError("keep_going must be boolean")
    output = args.output_directory or config.get("output_directory", "results")
    config["output_directory"] = str(resolve_path(output, base))
    config["cwd"] = str(resolve_path(config.get("cwd", "."), base))
    for category in ("profiles", "workloads"):
        entries = config.get(category)
        if not isinstance(entries, list) or not entries:
            raise ValueError(category + " must be a nonempty list")
        names = set()
        for entry in entries:
            name = entry.get("name", "")
            if not isinstance(name, str) or not SAFE_NAME.fullmatch(name):
                raise ValueError(category + " names must use letters, digits, . _ -")
            if name in names:
                raise ValueError("duplicate {} name: {}".format(category, name))
            names.add(name)
            if "timeout_seconds" in entry:
                entry["timeout_seconds"] = positive_number(
                    entry["timeout_seconds"], name + ".timeout_seconds")
            if category == "profiles":
                if entry.get("kind") not in ("cpu_mpi", "cuda_mpi"):
                    raise ValueError(name + ".kind must be cpu_mpi or cuda_mpi")
                argv = entry.get("argv")
                if (not isinstance(argv, list) or not argv or
                    not all(isinstance(token, str) and token for token in argv)):
                    raise ValueError(name + ".argv must be a nonempty string list")
                if not any("{input}" in token for token in argv):
                    raise ValueError(name + ".argv needs the {input} placeholder")
                environment = entry.setdefault("environment", {})
                if (not isinstance(environment, dict) or
                    not all(isinstance(k, str) and isinstance(v, str)
                            for k, v in environment.items())):
                    raise ValueError(name + ".environment must map strings to strings")
                entry["cwd"] = str(resolve_path(entry.get("cwd", config["cwd"]), base))
                if "expected_ranks" in entry:
                    integer(entry["expected_ranks"], name + ".expected_ranks", 1)
            else:
                entry["expected_count"] = integer(entry.get("expected_count"),
                                                   name + ".expected_count", 1)
                entry["input"] = str(resolve_path(entry["input"], base))
                entry["input_sha256"] = sha256(Path(entry["input"]))
                entry["input_bytes"] = Path(entry["input"]).stat().st_size
    return config


def inspect_result(kind, stdout, stderr, expected, expected_ranks=None):
    """Validate real CLI output; return count evidence or raise ValueError."""
    text = stdout + "\n" + stderr
    if kind == "cuda_mpi":
        codes = GPU_CODE.findall(text)
        counts = GPU_RESULT.findall(text)
        if codes != ["0"] or len(counts) != 1:
            raise ValueError("CUDA output needs exactly one code=0 and checked=N/N")
        checked, total = map(int, counts[0])
        if checked != expected or total != expected:
            raise ValueError("CUDA checked/total {}/{} != expected {}".format(
                checked, total, expected))
        return {"checked": checked, "total": total, "verification_code": 0}
    matches = CPU_RESULT.findall(text)
    if not matches:
        raise ValueError("native CPU output has no per-process count summaries")
    ranks = [int(item[0]) for item in matches]
    if len(ranks) != len(set(ranks)):
        raise ValueError("native CPU output repeats a process summary")
    if expected_ranks is not None and set(ranks) != set(range(expected_ranks)):
        raise ValueError("native CPU rank summaries {} != expected ranks 0..{}".format(
            sorted(ranks), expected_ranks - 1))
    checked = sum(int(item[3]) for item in matches)
    if checked != expected:
        raise ValueError("native CPU summed checked {} != expected {}".format(
            checked, expected))
    return {"checked": checked, "verification_code": 0,
            "processes": [{"rank": int(rank), "chunks": int(chunks),
                           "internal_seconds": float(seconds),
                           "checked": int(count)}
                          for rank, chunks, seconds, count in matches]}


def terminate_group(process, grace):
    """Signal only the process group created for this exact launcher."""
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=grace)
    except subprocess.TimeoutExpired:
        pass
    # A launcher can exit before its descendants. Kill remaining group members.
    try:
        os.killpg(process.pid, signal.SIGKILL)
    except ProcessLookupError:
        pass
    process.wait()


def run_one(profile, workload, phase, repetition, config, output):
    directory = output / workload["name"] / profile["name"]
    directory.mkdir(parents=True, exist_ok=True)
    stem = "{}-{:02d}".format(phase, repetition + 1)
    stdout_path = directory / (stem + ".stdout.log")
    stderr_path = directory / (stem + ".stderr.log")
    argv = [token.replace("{input}", workload["input"])
            for token in profile["argv"]]
    environment = os.environ.copy()
    environment.update(profile["environment"])
    timeout = workload.get("timeout_seconds", profile.get(
        "timeout_seconds", config["timeout_seconds"]))
    record = {"workload": workload["name"], "profile": profile["name"],
              "phase": phase, "repetition": repetition + 1, "argv": argv,
              "cwd": profile["cwd"], "environment_overrides": profile["environment"],
              "timeout_seconds": timeout, "stdout_log": str(stdout_path),
              "stderr_log": str(stderr_path), "timed_out": False,
              "interrupted": False, "valid": False}
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = None
        watchdog = None
        timeout_fired = threading.Event()
        start = time.perf_counter()
        try:
            process = subprocess.Popen(argv, cwd=profile["cwd"], env=environment,
                                       stdout=stdout, stderr=stderr,
                                       start_new_session=True, shell=False)

            def on_timeout():
                if process.poll() is None:
                    timeout_fired.set()
                    terminate_group(process, config["termination_grace_seconds"])

            # A separate watchdog allows waitpid to block normally. Using
            # Popen.wait(timeout=...) would poll, adding up to 50 ms of avoidable
            # latency to short benchmark measurements on POSIX systems.
            watchdog = threading.Timer(timeout, on_timeout)
            watchdog.daemon = True
            watchdog.start()
            process.wait()
        except KeyboardInterrupt:
            record["interrupted"] = True
            record["error"] = "comparison interrupted"
            if process is not None:
                terminate_group(process, config["termination_grace_seconds"])
        except OSError as error:
            record["error"] = "cannot launch command: " + str(error)
        finally:
            if watchdog is not None:
                watchdog.cancel()
                watchdog.join()
            record["wall_seconds"] = time.perf_counter() - start
            record["return_code"] = process.returncode if process is not None else None
            if timeout_fired.is_set():
                record["timed_out"] = True
                record["error"] = "command exceeded timeout"
    if (record["return_code"] == 0 and not record["timed_out"] and
            not record["interrupted"]):
        try:
            record["scientific_result"] = inspect_result(
                profile["kind"], stdout_path.read_text(errors="replace"),
                stderr_path.read_text(errors="replace"), workload["expected_count"],
                profile.get("expected_ranks"))
            record["valid"] = True
        except ValueError as error:
            record["error"] = str(error)
    elif "error" not in record:
        record["error"] = "command returned nonzero status"
    return record


def summarize(config, records):
    summaries = []
    for workload in config["workloads"]:
        for profile in config["profiles"]:
            group = [r for r in records if r["workload"] == workload["name"] and
                     r["profile"] == profile["name"]]
            measured = [r for r in group if r["phase"] == "measured" and r["valid"]]
            warmups = [r for r in group if r["phase"] == "warmup" and r["valid"]]
            qualified = (len(measured) == config["repetitions"] and
                         len(warmups) == config["warmups"] and
                         all(r["valid"] for r in group))
            entry = {"workload": workload["name"], "profile": profile["name"],
                     "qualified": qualified, "valid_measured_runs": len(measured),
                     "valid_warmup_runs": len(warmups)}
            if measured:
                seconds = [r["wall_seconds"] for r in measured]
                median = statistics.median(seconds)
                entry.update(median_wall_seconds=median,
                             min_wall_seconds=min(seconds),
                             max_wall_seconds=max(seconds),
                             median_tuples_per_second=workload["expected_count"] / median)
            summaries.append(entry)
    return summaries


def write_report(path, config, records, status):
    report = {"schema_version": 1,
              "generated_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "timing_scope": "complete mpiexec command including setup, I/O, cleanup",
              "status": status, "configuration": config,
              "runs": records, "summaries": summarize(config, records)}
    temporary = path.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("configuration", type=Path)
    parser.add_argument("--output-directory")
    parser.add_argument("--warmups", type=int)
    parser.add_argument("--repetitions", type=int)
    args = parser.parse_args()
    try:
        config_path = args.configuration.resolve()
        config = validate_config(json.loads(config_path.read_text()),
                                 config_path.parent, args)
        output = Path(config["output_directory"])
        output.mkdir(parents=True, exist_ok=True)
        report = output / "results.json"
        if report.exists():
            raise ValueError("output already contains results.json; choose a fresh directory")
    except (OSError, ValueError, TypeError, KeyError) as error:
        parser.error(str(error))
    records = []
    write_report(report, config, records, "running")
    for workload in config["workloads"]:
        for phase, count in (("warmup", config["warmups"]),
                             ("measured", config["repetitions"])):
            for repetition in range(count):
                profiles = config["profiles"]
                rotation = repetition % len(profiles)
                order = profiles[rotation:] + profiles[:rotation]
                for profile in order:
                    record = run_one(profile, workload, phase, repetition, config, output)
                    records.append(record)
                    write_report(report, config, records, "running")
                    print("{} {} {} {}: {:.3f}s {}".format(
                        workload["name"], profile["name"], phase, repetition + 1,
                        record["wall_seconds"], "PASS" if record["valid"] else
                        "FAIL " + record.get("error", "unknown error")), flush=True)
                    if record["interrupted"] or (not record["valid"] and
                                                   not config["keep_going"]):
                        write_report(report, config, records,
                                     "interrupted" if record["interrupted"] else "failed")
                        return 130 if record["interrupted"] else 1
    success = all(item["qualified"] for item in summarize(config, records))
    write_report(report, config, records, "complete" if success else "failed")
    print("Results: " + str(report), flush=True)
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())

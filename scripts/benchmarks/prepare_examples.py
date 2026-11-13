#!/usr/bin/env python3
"""Extract privately stored scientific .run inputs without changing contents.

This script records native-endian headers, exact weighted candidate counts,
and hashes. The directory names supply expected decisions, not validation.
Run check_input.c and both executables before using those labels as results.
Provide --examples or COMB_EXAMPLES_DIR; the compatibility default is the
ignored local examples/ directory. Circuit input data is not distributed here.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import sys


HEADER_NAMES = (
    "n_all", "n", "n_x", "n_z", "k", "fault_location_count", "s",
    "output_error_count", "x_measurement_count", "z_measurement_count",
)


class Reader:
    def __init__(self, data: bytes):
        self.data = data
        self.position = 0

    def u32(self) -> int:
        value, = struct.unpack_from("=I", self.data, self.position)
        self.position += 4
        return value

    def skip_bits(self, count: int) -> None:
        self.position += 4 * (1 + count // 4)
        if self.position > len(self.data):
            raise ValueError("truncated serialized bits")

    def vector(self, width: int) -> None:
        if self.u32() != width:
            raise ValueError("inconsistent vector dimension")
        self.skip_bits(width)

    def dcp(self, width: int) -> None:
        if self.u32() != width or self.u32() != width:
            raise ValueError("inconsistent DCP dimension")
        self.skip_bits(width * width)
        self.vector(width)
        self.vector(width)


def metadata(path: Path, data: bytes) -> dict:
    reader = Reader(data)
    h = {key: reader.u32() for key in HEADER_NAMES}
    for _ in range(h["n_x"]):
        reader.vector(h["n"])
    reader.skip_bits(1 << (h["n"] - h["n_x"]))
    for _ in range(h["n_z"]):
        reader.vector(h["n"])
    reader.skip_bits(1 << (h["n"] - h["n_z"]))
    for _ in range(h["k"]):
        reader.dcp(h["n_all"])
    partitions = [reader.u32() for _ in range(h["fault_location_count"])]
    if sum(partitions) != h["output_error_count"]:
        raise ValueError("partition sum does not equal output error count")
    for _ in range(h["output_error_count"]):
        reader.dcp(h["n_all"])
    for _ in range(h["x_measurement_count"] + h["z_measurement_count"]):
        reader.vector(h["n_all"])
    if reader.position != len(data):
        raise ValueError("trailing bytes")
    counts = [1] + [0] * h["s"]
    for partition in partitions:
        for selected in range(h["s"], 0, -1):
            counts[selected] += partition * counts[selected - 1]
    return dict(
        source=str(path),
        expected_code_from_directory=0 if path.parent.name == "fault-tolerant" else 1,
        uncompressed_sha256=hashlib.sha256(data).hexdigest(),
        compressed_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
        uncompressed_bytes=len(data),
        byteorder=sys.byteorder,
        header=h,
        candidate_count=counts[-1],
        location_subset_count=math.comb(h["fault_location_count"], h["s"]),
        native_mpi_chunk_delta=10000 // (15 ^ h["s"]),
        native_mpi_chunk_remainder=(
            math.comb(h["fault_location_count"], h["s"]) %
            (10000 // (15 ^ h["s"]))
        ),
        partitions=partitions,
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--examples", type=Path,
                        default=Path(os.environ.get(
                            "COMB_EXAMPLES_DIR",
                            Path(__file__).resolve().parents[2] / "examples")),
                        help="Private input directory (COMB_EXAMPLES_DIR or ignored examples/)")
    parser.add_argument("--names", nargs="+", help="Exact compressed filenames")
    args = parser.parse_args()
    sources = [source for source in sorted(args.examples.glob("*/*.run.gz"))
               if not args.names or source.name in args.names]
    if not sources:
        parser.error("No private .run.gz inputs found; supply --examples or COMB_EXAMPLES_DIR")
    if args.names and len(sources) != len(set(args.names)):
        parser.error("one or more requested example names were not found")
    args.output.mkdir(parents=True, exist_ok=True)
    entries = []
    for source in sources:
        data = gzip.decompress(source.read_bytes())
        entry = metadata(source, data)
        destination = args.output / source.with_suffix("").name
        destination.write_bytes(data)
        entry["output"] = str(destination)
        entries.append(entry)
        print(f"{destination.name}: n_all={entry['header']['n_all']}, "
              f"s={entry['header']['s']}, candidates={entry['candidate_count']:,}")
    (args.output / "inputs.json").write_text(json.dumps(entries, indent=2) + "\n")


if __name__ == "__main__":
    main()

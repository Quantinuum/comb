# Using comb from Python

We can invoke the comb library from Python using [cffi](https://cffi.readthedocs.io/en/stable/).

The [uv package manager](https://docs.astral.sh/uv/getting-started/installation/) is recommended for development of comb-py at this stage.

## Setting things up

Install the dependencies

```shell
uv sync --verbose
```

After we have installed the dependencies, we can run some demo scripts.

To run a script that demonstrates calling comb from Python using cffi, run the following from the `comb-py` directory

```shell
uv run src/comb_py/scratch.py
```

To run a script that calculates the pauli faults and the corresponding output errors, run

```shell
uv run src/comb_py/dcp.py
```
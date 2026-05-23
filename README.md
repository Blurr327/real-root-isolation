# Real Root Isolation

This project implements real root isolation for univariate integer
polynomials using exact arithmetic. The main algorithm is based on Descartes'
rule of signs, changes of variables, recursive subdivision, and FLINT
polynomial arithmetic.

The repository also contains a parallel version of the subdivision algorithm.
Parallelism is exposed in two places:

- OpenMP tasks for recursive subdivision;
- FLINT internal threads for dense polynomial operations such as Taylor shifts.

Experiments in this project showed that, for large dense polynomials, using
FLINT threads is usually more effective than relying only on OpenMP recursive
task parallelism.

## Repository Layout

```text
include/                 Public headers
src/                     Implementation files
test/                    Test and benchmark harnesses
CCA_Project_Report/      LaTeX report and generated PDF
testing/                 Extra scripts/data used during experimentation
vis/                     Visualization and plotting experiments
Makefile                 Build rules
```

Important source files:

- `src/rri_algo.c`: sequential subdivision algorithm.
- `src/parallel_rri_algo.c`: OpenMP/FLINT-threaded subdivision algorithm.
- `src/poly_utils.c`: polynomial transformations and utility functions.
- `src/fmpq_vec.c`: dynamic vector for rational interval endpoints.
- `test/test_rri_algo.c`: sequential root-isolation tests.
- `test/test_parallel_rri_algo.c`: parallel tests and timing harness.

## Dependencies

The project is written in C and depends on:

- GCC with OpenMP support;
- GMP;
- FLINT;
- Make;
- LaTeX tools, only if you want to rebuild the report.

On Ubuntu/Debian-like systems, the core dependencies are typically:

```sh
sudo apt install build-essential libgmp-dev libflint-dev
```

To rebuild the report:

```sh
sudo apt install texlive-latex-extra texlive-bibtex-extra biber latexmk
```

## Build

Build the main executable:

```sh
make
```

Build a specific test executable:

```sh
make test_rri_algo
make test_parallel_rri_algo
make test_poly_utils
```

Remove build artifacts:

```sh
make clean
```

By default the project is compiled with:

```text
-O3 -march=native -fopenmp
```

You can override optimization flags:

```sh
make clean
make OPTFLAGS="-O2"
```

## Running the Sequential Algorithm

Build and run the sequential test harness:

```sh
make test_rri_algo
./test_rri_algo --degree 20 --bits 16 --tests 5
```

This generates random square-free dense integer polynomials and checks that
the returned intervals isolate the real roots.

Useful options:

```text
--degree N      Degree of the random polynomial
--bits N        Coefficient bitsize
--tests N       Number of random tests
--random 0      Disable random generation and use --poly
--poly "..."    FLINT polynomial string
```

Example with an explicit polynomial:

```sh
./test_rri_algo --random 0 --poly "4  -8 0 106 -52"
```

FLINT polynomial strings start with the number of coefficients, followed by
coefficients in increasing degree order. Thus

```text
4  -8 0 106 -52
```

represents

```text
-8 + 0*x + 106*x^2 - 52*x^3.
```

## Running the Parallel Algorithm

Build and run the parallel benchmark harness:

```sh
make test_parallel_rri_algo
./test_parallel_rri_algo --degree 1000 --bits 32 --tests 3
```

The parallel version exposes two independent thread parameters:

```text
--threads N          Number of OpenMP subdivision threads
--flint-threads N    Number of FLINT internal arithmetic threads
```

Examples:

```sh
./test_parallel_rri_algo --degree 1000 --bits 32 --tests 3 --threads 4 --flint-threads 1
./test_parallel_rri_algo --degree 1000 --bits 32 --tests 3 --threads 2 --flint-threads 2
./test_parallel_rri_algo --degree 1000 --bits 32 --tests 3 --threads 1 --flint-threads 4
```

For larger dense inputs, the best observed configuration on the development
machine was often:

```sh
./test_parallel_rri_algo --degree 5000 --bits 32 --tests 1 --threads 1 --flint-threads 4
```

This is because the dominant cost moves into FLINT's dense polynomial
operations, especially Taylor shifts.

## Benchmarking with `perf`

Example:

```sh
perf stat ./test_parallel_rri_algo --degree 5000 --bits 32 --tests 1 --threads 1 --flint-threads 4
```

The report discusses CPU utilization, elapsed time, and thread-split behavior
for representative runs.

## Memory Checking with Valgrind

The Makefile provides Valgrind targets for test executables:

```sh
make vg_rri_algo
make vg_parallel_rri_algo
make vg_poly_utils
```

For example:

```sh
make vg_parallel_rri_algo
```

This runs:

```sh
valgrind --track-origins=yes --leak-check=full ./test_parallel_rri_algo
```

For shorter focused checks, build the test and pass arguments directly:

```sh
make test_parallel_rri_algo
valgrind --track-origins=yes --leak-check=full ./test_parallel_rri_algo --degree 10 --bits 16 --tests 1 --threads 1 --flint-threads 1
```

In the documented run, Valgrind reported `definitely lost: 0 bytes`. Some
`possibly lost` blocks were reported through GMP, FLINT, and OpenMP call
stacks, so the result is interpreted as no confirmed leak in the project code.

## Report

The project report is in:

```text
CCA_Project_Report/main.tex
CCA_Project_Report/main.pdf
```

To rebuild it:

```sh
cd CCA_Project_Report
latexmk -pdf main.tex
```

The report contains:

- notation and problem statement;
- input/output size analysis;
- Descartes' rule of signs;
- subdivision algorithm and correctness;
- complexity discussion;
- implementation notes;
- experimental timings and parallelization results.

## Notes and Limitations

- The implementation targets square-free input polynomials. The test harnesses
  generate square-free random polynomials.
- Root isolation is performed with exact rational endpoints.
- Coefficient growth after changes of variables is a major practical
  bottleneck.
- The `real_root_isolation` executable currently runs Taylor-shift
  experimentation code from `src/real_root_isolation.c`; the most useful
  entry points for root isolation are the test harnesses listed above.

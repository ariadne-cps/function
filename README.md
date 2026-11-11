# Ariadne Function

[![Format](https://github.com/ariadne-cps/function/actions/workflows/format.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/format.yml)
[![macOS Status](https://github.com/ariadne-cps/function/actions/workflows/macos.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/macos.yml)
[![Ubuntu Status](https://github.com/ariadne-cps/function/actions/workflows/ubuntu.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/ubuntu.yml)
[![Windows Status](https://github.com/ariadne-cps/function/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/function/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/function/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/function/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/function)

Ariadne Function is the standalone C++20 function layer of Ariadne. It provides symbolic, affine, validated and approximate functions, interfaces and function models, Taylor models, procedures, domains, and multifunctions.
Although independently maintained and usable, this repository is also part of the [Ariadne repository structure](https://github.com/ariadne-cps/ariadne/blob/main/doc/structure.md).

## Layout

Public headers live under `include/function`, implementations under `src`, C++ tests under `tests`, Python bindings and tests under `python`, and standalone examples under `tutorials`.

The umbrella header is `include/ariadne-function.hpp`.

## Dependencies

Function depends directly on:

- [ariadne-cps/geometry](https://github.com/ariadne-cps/geometry), included as a Git submodule.
- [ariadne-cps/build](https://github.com/ariadne-cps/build), used for shared build configuration.

Geometry supplies Algebra, Interval, Numeric and other lower-level Ariadne libraries transitively. Function does **not** depend on Solving, Set, Paving or BDD.

## Build

```bash
git clone --recurse-submodules https://github.com/ariadne-cps/function.git
cd function
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel --target tests
ctest --test-dir build --output-on-failure
```

For a release build and installation:

```bash
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel --target ariadne-function bindings
cmake --install build-release --prefix /path/to/install
```

## Python bindings

With a compatible Python >= 3.11 development environment, the standalone build produces the `pyariadne` extension by combining Function bindings with the transitive Geometry, Algebra, Interval and Numeric bindings.

The extension installs under the chosen prefix (`lib/pythonX.Y/site-packages` on Unix, `Lib/site-packages` on Windows), not into a global Python installation.

## Tutorials

- C++: `tutorial/tutorial_function.cpp`, built separately against an installed `ariadne-function`.
- Python: `python/tutorial.py`.

## Contribution guidelines

If you would like to contribute to Ariadne Function, please contact the developer:

- Luca Geretti <luca.geretti@univr.it>

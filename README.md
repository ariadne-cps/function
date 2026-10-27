# Ariadne Function

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/function/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/function/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/function/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/function/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/ariadne-cps/function/branch/main/graph/badge.svg)](https://codecov.io/gh/ariadne-cps/function)

Ariadne Function is the standalone C++20 function layer of Ariadne. It provides symbolic, affine, validated and approximate functions, interfaces and function models, Taylor models, procedures, domains, and multifunctions.

## Layout

Public headers live under `include/function`, implementations under `src`, C++ tests under `tests`, Python bindings and tests under `python`, and standalone examples under `tutorials`.

The umbrella header is `include/ariadne-function.hpp`.

## Dependencies

Function depends directly on:

- [ariadne-cps/geometry](https://github.com/ariadne-cps/geometry), included as a Git submodule.
- [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration), used for shared build configuration.

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

- C++: `tutorials/tutorial_function/tutorial_function.cpp`, built separately against an installed `ariadne-function`.
- Python: `python/tutorials/tutorial_function.py`.

## Extraction

Source files and tests were extracted from `ariadne-cps/kernel`, branch `split`. Function-based constrained-image sets remain in the **Set** module and were not imported. Replacing Function inside Kernel is a separate integration step.

## Contribution guidelines

If you would like to contribute to Ariadne Function, please contact the developer:

- Luca Geretti <luca.geretti@univr.it>

## License

GNU General Public License v3.0 (or later where individual source files specify it).

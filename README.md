# Ariadne Function

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Unix Status](https://github.com/ariadne-cps/function/actions/workflows/unix.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/unix.yml)
[![Windows Status](https://github.com/ariadne-cps/function/actions/workflows/win.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/win.yml)
[![Sanitizers Status](https://github.com/ariadne-cps/function/actions/workflows/sanitizers.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/sanitizers.yml)
[![Coverage Status](https://github.com/ariadne-cps/function/actions/workflows/coverage.yml/badge.svg)](https://github.com/ariadne-cps/function/actions/workflows/coverage.yml)

Standalone C++20 function layer of Ariadne: symbolic, affine, validated and approximate functions, function interfaces, Taylor models, procedures, domains and generic multifunctions.

## Dependencies

Function depends directly on [ariadne-cps/geometry](https://github.com/ariadne-cps/geometry) and uses [ariadne-cps/configuration](https://github.com/ariadne-cps/configuration) for the build. Geometry provides Algebra, Interval, Numeric and other lower libraries transitively. Function does **not** depend on Solving, Set, Paving or BDD.

## Build and test

```sh
git clone --recurse-submodules https://github.com/ariadne-cps/function.git
cd function
cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug
cmake --build build --parallel --target tests
ctest --test-dir build --output-on-failure
```

## Install

```sh
cmake -S . -B build-release -DCMAKE_BUILD_TYPE=Release
cmake --build build-release --parallel --target ariadne-function bindings
cmake --install build-release --prefix /path/to/install
```

The C++ library installs as `ariadne-function` along with the public headers and transitive dependency bundle. With Python >= 3.11, the installable `pyariadne` extension contains Function and lower-level Geometry/Algebra/Interval/Numeric bindings.

## Tutorials

- C++: `tutorials/tutorial_function/tutorial_function.cpp` (build separately against installed `ariadne-function`).
- Python: `python/tutorials/tutorial_function.py`.

## Layering

The source code and original tests were extracted from `ariadne-cps/kernel`, branch `split`. Function-based constrained-image set implementations remain in the **Set** module and were not imported. Replacing the Function module inside Kernel is deliberately deferred.

## License

GNU General Public License v3.0 (or later where individual sources specify it).

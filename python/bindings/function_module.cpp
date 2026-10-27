#include "pybind11.hpp"

void function_submodule(pybind11::module& module);

void function_module(pybind11::module& module) {
    function_submodule(module);
}

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "repair.hpp"
#include "runner.hpp"


namespace py = pybind11;


PYBIND11_MODULE(engine, m) {
    m.doc() = "c++ backend engine";
    m.def("run", &run, "run algo");
}
#include <pybind11/pybind11.h>
namespace py = pybind11;

char const* greet()
{
   return "hello, world";
}

PYBIND11_MODULE(gressure, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("greet", &greet, "A function which adds two numbers");
}

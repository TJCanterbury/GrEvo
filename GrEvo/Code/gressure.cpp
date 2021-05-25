#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include <pybind11/stl.h>
#include <pybind11/pybind11.h>
#include <iostream>

using namespace std;
namespace py = pybind11;

char const* greet()
{
   return "hello, world";
}

struct graphEdge {
    string start_ver, end_ver;
};

py::object print_vector( list<string> &v) {
    
    for ( string &s : v){
        if (s[0] == 'L')
            s[0] = 'R';
        else if (s[0] == 'R')
        {
            s[0] = 'L';
        }
    }

    py::object obj = py::cast(v);
    return v;
}

PYBIND11_MODULE(gressure, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("greet", &greet, "Says hello");
    m.def("print_vector", &print_vector, "prints a given python list from C++");
}

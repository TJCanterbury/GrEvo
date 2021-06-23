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

py::object Reflect( std::list<string> &v) {
    
    for ( string &s : v){
        if (s[0] == 'L')
            s[0] = 'R';
        else if (s[0] == 'R')
        {
            s[0] = 'L';
        }
    }
    py::list v2 = py::cast(v);
    return v2;
}

PYBIND11_MODULE(gressure, m)
{
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("greet", &greet, "Says hello");
    m.def("Reflect", &Reflect, "prints a given python list from C++");
}

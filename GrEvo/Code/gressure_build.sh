c++ -O3 -Wall -shared -std=c++11 -fPIC $(python3 -m pybind11 --includes) gressure.cpp -o gressure$(python3-config --extension-suffix)

python3 test.py
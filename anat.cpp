#include <iostream>
using namespace std;

void struct

double square(double x)
{
    return x*x;
}

void print_square(double x)
{
    cout << "The square of " << x << " is " << square(x) << "\n";
}

int main()
{
    cout << "Hello World!\n";
    print_square(1.234);
}
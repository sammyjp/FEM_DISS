#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"

int main(int argc, char* argv[])
{
    Vector* Gauss = new Vector (3);

    QuadratureLibrary().GetGaussQuadraturePoints(*Gauss);

    std::cout << *Gauss << std::endl;

    return 0;
}

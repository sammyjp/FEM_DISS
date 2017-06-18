#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"

int main(int argc, char* argv[])
{
    Matrix* GridPoints = new Matrix(2,5);
    (*GridPoints)(1,1) = 0;
    (*GridPoints)(1,2) = 1;
    (*GridPoints)(1,3) = 0;
    (*GridPoints)(1,4) = 0.5;
    (*GridPoints)(1,5) = 0.75;
    (*GridPoints)(2,1) = 0;
    (*GridPoints)(2,2) = 0;
    (*GridPoints)(2,3) = 0.75;
    (*GridPoints)(2,4) = 1;
    (*GridPoints)(2,5) = 0.5;


    std::cout << *GridPoints;

    Matrix* Connectivity = new Matrix(2,4);
    (*Connectivity)(1,1) = 1;
    (*Connectivity)(1,2) = 2;
    (*Connectivity)(1,3) = 3;
    (*Connectivity)(1,4) = 0;
    (*Connectivity)(2,1) = 2;
    (*Connectivity)(2,2) = 5;
    (*Connectivity)(2,3) = 4;
    (*Connectivity)(2,4) = 3;
    Mesh* M = new Mesh(*GridPoints, *Connectivity);

    std::cout << M->GetElement(0)->GetElementConnectivityArray() << std::endl;

    delete GridPoints;
    delete Connectivity;
    delete M;

    return 0;
}

#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"

int main(int argc, char* argv[])
{
    Matrix* Nodes = new Matrix(2,3);
    (*Nodes)(1,1) = 0;
    (*Nodes)(1,2) = 1;
    (*Nodes)(1,3) = 0;
    (*Nodes)(2,1) = 0;
    (*Nodes)(2,2) = 0;
    (*Nodes)(2,3) = 1;

    Matrix* GridPoints = new Matrix(2,3);
    (*GridPoints)(1,1) = 0;
    (*GridPoints)(1,2) = 1;
    (*GridPoints)(1,3) = 0;
    (*GridPoints)(2,1) = 0;
    (*GridPoints)(2,2) = 0;
    (*GridPoints)(2,3) = 1;

    Vector* Connectivity = new Vector(3);
    (*Connectivity)[0] = 1;
    (*Connectivity)[1] = 2;
    (*Connectivity)[2] = 3;

    Triangle* Tri = new Triangle(*Connectivity);

    Matrix* globalPoints = new Matrix(2,3);

    Tri->MapGlobalToLocal(*Nodes, *GridPoints, *globalPoints);

    std::cout << *globalPoints << std::endl;

    delete Nodes;
    delete GridPoints;
    delete Tri;

    return 0;
}

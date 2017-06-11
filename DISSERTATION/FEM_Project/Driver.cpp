#include <iostream>

#include "Mesh.hpp"
#include "FE_Solution.hpp"

int main(int argc, char* argv[])
{
    int dim = 1;
    int elem = 2;
    double h = 0.2;

    Mesh* H = new Mesh(elem, dim);

    std::cout << "Elements = " << H->GetNumElements() << std::endl;
    std::cout << "Nodes = " << H->GetNumNodes() << std::endl;
    std::cout << "Dimension = " << H->GetDimension() << std::endl;

    H->GenerateUniformMesh();

    std::cout << "X = " << std::endl << H->GetXGridPoints();

    FE_Solution* Sol = new FE_Solution(*H);

    Vector* x = new Vector(5);
    Vector* F = new Vector(5);
    (*x)[0] = 0;
    (*x)[1] = 0.25;
    (*x)[2] = 0.5;
    (*x)[3] = 0.75;
    (*x)[4] = 1;

    Sol->ComputeLinearBasisFunctionValues(1, *F,*x);
    std::cout << "F " << std::endl << *F;


    delete H;
    delete Sol;
    return 0;
}

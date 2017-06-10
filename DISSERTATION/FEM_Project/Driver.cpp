#include <iostream>

#include "Mesh.hpp"

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

    delete H;
    return 0;
}

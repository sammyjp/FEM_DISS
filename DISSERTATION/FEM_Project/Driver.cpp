#include <iostream>
#include <cmath>

#include "Matrix.hpp"
#include "Mesh.hpp"
#include "FE_Solution.hpp"

int main(int argc, char* argv[])
{
    int dim = 1;
    int elem = 2;

    Mesh* H = new Mesh(elem, dim);

    std::cout << "Elements = " << H->GetNumElements() << std::endl;
    std::cout << "Nodes = " << H->GetNumNodes() << std::endl;
    std::cout << "Dimension = " << H->GetDimension() << std::endl;


    Matrix* M = new Matrix(3,5);
    for (int i=1; i<=M->GetNumberOfRows(); i++)
    {
        for (int j=1; j<=M->GetNumberOfColumns(); j++)
        {
            (*M)(i,j) = pow(2*i+j,2.0);
        }
    }
    std::cout << "Matrix M =" << std::endl << *M;
    std::cout << "Row 2 =" << std::endl << M->GetRowAsVector(2);
    std::cout << "Column 4 =" << std::endl << M->GetColumnAsVector(4);
    delete M;


    delete H;
    return 0;
}

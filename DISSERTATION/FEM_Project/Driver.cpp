#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "SparseMatrix.hpp"

int main(int argc, char* argv[])
{

    Matrix* HELLO = new Matrix(3,4);
    (*HELLO)(1,1) = 1;
    (*HELLO)(1,2) = 0;
    (*HELLO)(1,3) = 3;
    (*HELLO)(1,4) = 1;
    (*HELLO)(2,1) = 0;
    (*HELLO)(2,2) = 0;
    (*HELLO)(2,3) = 2;
    (*HELLO)(2,4) = 0;
    (*HELLO)(3,1) = 1;
    (*HELLO)(3,2) = 5;
    (*HELLO)(3,3) = 0;
    (*HELLO)(3,4) = 1;

    Vector* vec = new Vector(4);
    (*vec)[0] = 1;
    (*vec)[1] = 2;
    (*vec)[2] = 3;
    (*vec)[3] = 4;

    SparseMatrix* HI = new SparseMatrix(*HELLO);

    for (int i=1; i<=HELLO->GetNumberOfRows(); i++)
    {
        for (int j=1; j<=HELLO->GetNumberOfColumns(); j++)
        {
            std::cout << HI->Read(i,j) << " ";
        }
        std::cout << std::endl;
    }

    return 0;

    Matrix* gridPoints = new Matrix (1,3);
    (*gridPoints)(1,1) = 0;
    (*gridPoints)(1,2) = 0.5;
    (*gridPoints)(1,3) = 2;

    Matrix* connectivity = new Matrix (2,2);
    (*connectivity)(1,1) = 1;
    (*connectivity)(1,2) = 2;
    (*connectivity)(2,1) = 2;
    (*connectivity)(2,2) = 3;

    Mesh* M = new Mesh (*gridPoints, *connectivity);

    for (int i=0; i<=1; i++)
    {
        Vector* quadraturePoints = new Vector (5);
        Vector* functionPoints = new Vector (5);

        Matrix* elementNodes = new Matrix (1,2);
        Matrix* jacobian = new Matrix(1,1);

        (*elementNodes)(1,1) = (*gridPoints)(1,i+1);
        (*elementNodes)(1,2) = (*gridPoints)(1,i+2);

        M->GetElement(i)->ComputeElementQuadraturePoints(*quadraturePoints);
        M->GetElement(i)->ComputeMappingJacobian(*elementNodes, *quadraturePoints, *jacobian);

        QuadratureLibrary().TransformGQPoints((*elementNodes)(1,1),(*elementNodes)(1,2),*quadraturePoints);
        for (int j=0; j<quadraturePoints->GetSize(); j++)
        {
            (*functionPoints)[j] = (*quadraturePoints)[j];
        }

        std::cout << M->GetElement(i)->PerformElementQuadrature(*quadraturePoints, *functionPoints, *jacobian) << std::endl;

        delete quadraturePoints;
        delete functionPoints;
        delete elementNodes;
        delete jacobian;
    }

    delete M;
    delete gridPoints;
    delete connectivity;

    return 0;
}

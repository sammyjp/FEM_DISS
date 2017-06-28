#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"

int main(int argc, char* argv[])
{
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

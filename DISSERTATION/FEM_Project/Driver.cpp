#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "SparseMatrix.hpp"

int main(int argc, char* argv[])
{
    int numElements = 1;

    Matrix* Grid = new Matrix (2,4);
    (*Grid)(1,1) = 0;
    (*Grid)(2,1) = 0;
    (*Grid)(1,2) = 1;
    (*Grid)(2,2) = 0;
    (*Grid)(1,3) = 0;
    (*Grid)(2,3) = 1;
    (*Grid)(1,4) = 1;
    (*Grid)(2,4) = 1;

    Vector* Connectivity1 = new Vector (4);
    (*Connectivity1)(1) = 1;
    (*Connectivity1)(2) = 2;
    (*Connectivity1)(3) = 4;
    (*Connectivity1)(4) = 3;

    Mesh* myMesh = new Mesh(*Grid, numElements);

    myMesh->InitialiseElement(1, *Connectivity1, 2);

    int n_q = 9;
    Vector* quadratureWeights = new Vector (n_q);
    Matrix* localQuadraturePoints = new Matrix (2,n_q);
    Matrix* globalQuadraturePoints = new Matrix (2,n_q);
    double I;

    Vector* functionPoints = new Vector (n_q);

    for (int i=1; i<=numElements; i++)
    {
        myMesh->GetElement(i)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
        for (int j=1; j<=functionPoints->GetSize(); j++)
        {
            (*functionPoints)(j) = pow((*globalQuadraturePoints)(1,j),2.0)+pow((*globalQuadraturePoints)(2,j),2.0);
        }

        I = myMesh->GetElement(i)->PerformElementQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *functionPoints);

        std::cout << "Element " << i << ":" << std::endl << std::endl;
        std::cout << "Gauss weights =\n" << *quadratureWeights;
        std::cout << "Local Gauss points =\n" << *localQuadraturePoints;
        std::cout << "Global Gauss points =\n" << *globalQuadraturePoints;
        std::cout << "Integral of f(x) = x on element = " << I << std::endl;
        std::cout << std::endl;
    }


    delete Grid;
    delete Connectivity1;
    delete myMesh;
    delete quadratureWeights;
    delete localQuadraturePoints;
    delete globalQuadraturePoints;
    delete functionPoints;

    return 0;
}

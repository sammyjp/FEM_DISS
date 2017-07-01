#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "SparseMatrix.hpp"

int main(int argc, char* argv[])
{
    int numElements = 3;

    Matrix* Grid = new Matrix (1,4);
    (*Grid)(1,1) = 0;
    (*Grid)(1,2) = 1;
    (*Grid)(1,3) = 2;
    (*Grid)(1,4) = 3;

    Vector* Connectivity1 = new Vector (2);
    Vector* Connectivity2 = new Vector (2);
    Vector* Connectivity3 = new Vector (2);
    (*Connectivity1)(1) = 1;
    (*Connectivity1)(2) = 2;
    (*Connectivity2)(1) = 2;
    (*Connectivity2)(2) = 3;
    (*Connectivity3)(1) = 3;
    (*Connectivity3)(2) = 4;

    Mesh* myMesh = new Mesh(*Grid, numElements);

    myMesh->InitialiseElement(1, *Connectivity1, 0);
    myMesh->InitialiseElement(2, *Connectivity2, 0);
    myMesh->InitialiseElement(3, *Connectivity3, 0);

    int n_q = 3;
    Vector* quadratureWeights = new Vector (n_q);
    Matrix* quadraturePoints = new Matrix (1,n_q);
    Matrix* jacobian = new Matrix (1,1);
    Vector* jacobianDeterminant = new Vector (n_q);
    double I;

    Vector* functionPoints = new Vector (n_q);

    for (int i=1; i<=numElements; i++)
    {
        I = 0;
        myMesh->GetElement(i)->GetQuadrature(n_q, *quadratureWeights, *quadraturePoints);
        for (int j=1; j<=functionPoints->GetSize(); j++)
        {
            (*functionPoints)(j) = pow((*quadraturePoints)(1,j),2.0);

            Vector* jacPoints = new Vector(quadraturePoints->GetColumnAsVector(j));
            myMesh->GetElement(i)->ComputeMappingJacobian(*jacPoints, *jacobian);
            delete jacPoints;
            (*jacobianDeterminant)(j) = jacobian->CalculateDeterminant();
        }
        for (int j=1; j<=quadratureWeights->GetSize(); j++)
        {
            I += (*functionPoints)(j)*(*quadratureWeights)(j)*(*jacobianDeterminant)(j);
        }

        std::cout << "Element " << i << ":" << std::endl << std::endl;
        std::cout << "Gauss weights =\n" << *quadratureWeights;
        std::cout << "Gauss points =\n" << *quadraturePoints;
        std::cout << "Jacobian determinant = " << jacobianDeterminant << std::endl;
        std::cout << "Integral of f(x) = x on element = " << I << std::endl;
        std::cout << std::endl;
    }


    delete Grid;
    delete Connectivity1;
    delete Connectivity2;
    delete Connectivity3;
    delete myMesh;
    delete quadratureWeights;
    delete quadraturePoints;
    delete functionPoints;

    return 0;
}

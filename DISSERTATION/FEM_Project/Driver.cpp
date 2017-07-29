#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

int main(int argc, char* argv[])
{
    int numElements = 4;

    Matrix* Grid = new Matrix (1,numElements+1);
    (*Grid)(1,1) = 0;
    (*Grid)(1,2) = 0.25;
    (*Grid)(1,3) = 0.5;
    (*Grid)(1,4) = 0.75;
    (*Grid)(1,5) = 1;

//    (*Grid)(1,1) = 0;
//    (*Grid)(1,2) = 0.125;
//    (*Grid)(1,3) = 0.25;
//    (*Grid)(1,4) = 0.375;
//    (*Grid)(1,5) = 0.5;
//    (*Grid)(1,6) = 0.625;
//    (*Grid)(1,7) = 0.75;
//    (*Grid)(1,8) = 0.875;
//    (*Grid)(1,9) = 1;

    Vector* Connectivity1 = new Vector (2);
    (*Connectivity1)(1) = 1;
    (*Connectivity1)(2) = 2;
    Vector* Connectivity2 = new Vector (2);
    (*Connectivity2)(1) = 2;
    (*Connectivity2)(2) = 3;
    Vector* Connectivity3 = new Vector (2);
    (*Connectivity3)(1) = 3;
    (*Connectivity3)(2) = 4;
    Vector* Connectivity4 = new Vector (2);
    (*Connectivity4)(1) = 4;
    (*Connectivity4)(2) = 5;
//    Vector* Connectivity5 = new Vector (2);
//    (*Connectivity5)(1) = 5;
//    (*Connectivity5)(2) = 6;
//    Vector* Connectivity6 = new Vector (2);
//    (*Connectivity6)(1) = 6;
//    (*Connectivity6)(2) = 7;
//    Vector* Connectivity7 = new Vector (2);
//    (*Connectivity7)(1) = 7;
//    (*Connectivity7)(2) = 8;
//    Vector* Connectivity8 = new Vector (2);
//    (*Connectivity8)(1) = 8;
//    (*Connectivity8)(2) = 9;

    Mesh* myMesh = new Mesh(*Grid, numElements);

    myMesh->InitialiseElement(1, *Connectivity1, 0);
    myMesh->InitialiseElement(2, *Connectivity2, 0);
    myMesh->InitialiseElement(3, *Connectivity3, 0);
    myMesh->InitialiseElement(4, *Connectivity4, 0);
//    myMesh->InitialiseElement(5, *Connectivity5, 0);
//    myMesh->InitialiseElement(6, *Connectivity6, 0);
//    myMesh->InitialiseElement(7, *Connectivity7, 0);
//    myMesh->InitialiseElement(8, *Connectivity8, 0);


    FE_Solution* FE = new FE_Solution(*myMesh, 1);

    SparseMatrix* A = new SparseMatrix(*FE, numElements);
    Vector* F = new Vector (FE->GetNumberOfDofs());

    std::cout << A->GetColumnIndexArray();
    return 0;

    int n_q = 5;

    int numElementDofs;

    for (int k=1; k<=myMesh->GetNumElements(); k++)
    {
        Vector* quadratureWeights = new Vector (n_q);
        Matrix* localQuadraturePoints = new Matrix (1,n_q);
        Matrix* globalQuadraturePoints = new Matrix (1,n_q);

        numElementDofs = FE->GetElementPolynomialSpace(k)->GetNumElementDofs();
        Matrix* A_loc = new Matrix (numElementDofs, numElementDofs);
        Vector* f_loc = new Vector (numElementDofs);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);

        for (int q=1; q<=n_q; q++)
        {
            Vector* basisValues = new Vector(FE->GetElementPolynomialSpace(1)->GetNumElementDofs());
            Matrix* basisGrad = new Matrix(myMesh->GetDimension(),FE->GetElementPolynomialSpace(1)->GetNumElementDofs());

            FE->GetElementPolynomialSpace(k)->ComputeBasis((*localQuadraturePoints)(1,q), *basisValues);
            FE->GetElementPolynomialSpace(k)->ComputeGradBasis((*localQuadraturePoints)(1,q), *basisGrad);

            for (int i=1; i<=A_loc->GetNumberOfRows(); i++)
            {
                for (int j=1; j<=A_loc->GetNumberOfColumns(); j++)
                {
                    (*A_loc)(i,j) += (*basisGrad)(1,i)*(*basisGrad)(1,j)*(*quadratureWeights)(q);
                }
            }
            for (int i=1; i<=f_loc->GetSize(); i++)
            {
                (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q);
            }

            delete basisValues;
            delete basisGrad;
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;

        for (int i=1; i<=A_loc->GetNumberOfRows(); i++)
        {
            for (int j=1; j<=A_loc->GetNumberOfColumns(); j++)
            {
                A->AddValue((*A_loc)(i,j), (FE->GetElementDofs(k))(i), (FE->GetElementDofs(k))(j));
            }
        }
        for (int i=1; i<=f_loc->GetSize(); i++)
        {
            (*F)((FE->GetElementDofs(k))(i)) += (*f_loc)(i);
        }

        delete A_loc;
        delete f_loc;
    }

    std::cout << *A;
    std::cout << *F;

    Vector* Sol = new Vector(FE->GetNumberOfDofs());

    A->CGSolveSystem(*F, *Sol, 1e-9, 1000);
    std::cout << *Sol;

    return 0;
}

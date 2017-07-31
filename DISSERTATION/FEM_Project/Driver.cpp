#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

void Example1D();
void Example2D();


int main(int argc, char* argv[])
{
    Example2D();
    return 0;
}

void Example1D()
{
    int numElements = 8;

    Matrix* Grid = new Matrix (1,numElements+1);
//    (*Grid)(1,1) = 0;
//    (*Grid)(1,2) = 0.25;
//    (*Grid)(1,3) = 0.5;
//    (*Grid)(1,4) = 0.75;
//    (*Grid)(1,5) = 1;

    (*Grid)(1,1) = 0;
    (*Grid)(1,2) = 0.125;
    (*Grid)(1,3) = 0.25;
    (*Grid)(1,4) = 0.375;
    (*Grid)(1,5) = 0.5;
    (*Grid)(1,6) = 0.625;
    (*Grid)(1,7) = 0.75;
    (*Grid)(1,8) = 0.875;
    (*Grid)(1,9) = 1;

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
    Vector* Connectivity5 = new Vector (2);
    (*Connectivity5)(1) = 5;
    (*Connectivity5)(2) = 6;
    Vector* Connectivity6 = new Vector (2);
    (*Connectivity6)(1) = 6;
    (*Connectivity6)(2) = 7;
    Vector* Connectivity7 = new Vector (2);
    (*Connectivity7)(1) = 7;
    (*Connectivity7)(2) = 8;
    Vector* Connectivity8 = new Vector (2);
    (*Connectivity8)(1) = 8;
    (*Connectivity8)(2) = 9;

    Mesh* myMesh = new Mesh(*Grid, numElements);

    myMesh->InitialiseElement(1, *Connectivity1, 0);
    myMesh->InitialiseElement(2, *Connectivity2, 0);
    myMesh->InitialiseElement(3, *Connectivity3, 0);
    myMesh->InitialiseElement(4, *Connectivity4, 0);
    myMesh->InitialiseElement(5, *Connectivity5, 0);
    myMesh->InitialiseElement(6, *Connectivity6, 0);
    myMesh->InitialiseElement(7, *Connectivity7, 0);
    myMesh->InitialiseElement(8, *Connectivity8, 0);


    FE_Solution* FE = new FE_Solution(*myMesh, 1);

    SparseMatrix* A = new SparseMatrix(*FE, numElements);
    Vector* F = new Vector (FE->GetNumberOfDofs());

    int n_q = 4;

    int numElementDofs;

    for (int k=1; k<=myMesh->GetNumElements(); k++)
    {
        Vector* quadratureWeights = new Vector (n_q);
        Matrix* localQuadraturePoints = new Matrix (1,n_q);
        Matrix* globalQuadraturePoints = new Matrix (1,n_q);

        numElementDofs = FE->GetNumElementDofs(k);
        Matrix* A_loc = new Matrix (numElementDofs, numElementDofs);
        Vector* f_loc = new Vector (numElementDofs);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);


        for (int q=1; q<=n_q; q++)
        {
            Vector* basisValues = new Vector(numElementDofs);
            Matrix* basisGrad = new Matrix(myMesh->GetDimension(),numElementDofs);
            Matrix* jacobian = new Matrix (myMesh->GetDimension(), myMesh->GetDimension());
            Vector* pointToEval = new Vector(localQuadraturePoints->GetColumnAsVector(q));

            myMesh->GetElement(k)->ComputeMappingJacobian(*pointToEval, *jacobian);
            delete pointToEval;

            FE->ComputeBasis(k, (*localQuadraturePoints)(1,q), *basisValues);
            FE->ComputeGradBasis(k, (*localQuadraturePoints)(1,q), *basisGrad);

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
            delete jacobian;
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;

        for (int i=1; i<=A_loc->GetNumberOfRows(); i++)
        {
            for (int j=1; j<=A_loc->GetNumberOfColumns(); j++)
            {
                A->AddValue(A->Read((FE->GetElementDofs(k))(i), (FE->GetElementDofs(k))(j)) + (*A_loc)(i,j), (FE->GetElementDofs(k))(i), (FE->GetElementDofs(k))(j));
            }
        }
        for (int i=1; i<=f_loc->GetSize(); i++)
        {
            (*F)((FE->GetElementDofs(k))(i)) += (*f_loc)(i);
        }

        delete A_loc;
        delete f_loc;
    }

    for (int i=1; i<=A->GetNumberOfRows(); i++)
    {
        for (int j=1; j<=A->GetNumberOfColumns(); j++)
        {
            if (i==1 || i==myMesh->GetNumNodes() || j==1 || j==myMesh->GetNumNodes())
            {
                if (i==j)
                {
                    A->AddValue(1, i, j);
                }
                else
                {
                    A->AddValue(0, i, j);
                }
            }
        }
    }
    for (int i=1; i<=F->GetSize(); i++)
    {
        if (i==1 || i==myMesh->GetNumNodes())
        {
            (*F)(i) = 0;
        }
    }
    std::cout << "A = \n" << *A;
    std::cout << "F = \n" << *F;

    A->CGSolveSystem(*F, FE->GetSolutionVector(), 1e-9, 1000);
    std::cout << "U = \n" << FE->GetSolutionVector();
}

void Example2D()
{
    int numElements = 4;

    Matrix* Grid = new Matrix (2,9);

    (*Grid)(1,1) = 0;
    (*Grid)(2,1) = 0;
    (*Grid)(1,2) = 0.5;
    (*Grid)(2,2) = 0;
    (*Grid)(1,3) = 1;
    (*Grid)(2,3) = 0;
    (*Grid)(1,4) = 0;
    (*Grid)(2,4) = 0.5;
    (*Grid)(1,5) = 0.5;
    (*Grid)(2,5) = 0.5;
    (*Grid)(1,6) = 1;
    (*Grid)(2,6) = 0.5;
    (*Grid)(1,7) = 0;
    (*Grid)(2,7) = 1;
    (*Grid)(1,8) = 0.5;
    (*Grid)(2,8) = 1;
    (*Grid)(1,9) = 1;
    (*Grid)(2,9) = 1;

    Vector* Connectivity1 = new Vector (4);
    (*Connectivity1)(1) = 1;
    (*Connectivity1)(2) = 2;
    (*Connectivity1)(3) = 5;
    (*Connectivity1)(4) = 4;
    Vector* Connectivity2 = new Vector (4);
    (*Connectivity2)(1) = 2;
    (*Connectivity2)(2) = 3;
    (*Connectivity2)(3) = 6;
    (*Connectivity2)(4) = 5;
    Vector* Connectivity3 = new Vector (4);
    (*Connectivity3)(1) = 4;
    (*Connectivity3)(2) = 5;
    (*Connectivity3)(3) = 8;
    (*Connectivity3)(4) = 7;
    Vector* Connectivity4 = new Vector (4);
    (*Connectivity4)(1) = 5;
    (*Connectivity4)(2) = 6;
    (*Connectivity4)(3) = 9;
    (*Connectivity4)(4) = 8;

    Mesh* myMesh = new Mesh(*Grid, numElements);

    myMesh->InitialiseElement(1, *Connectivity1, 2);
    myMesh->InitialiseElement(2, *Connectivity2, 2);
    myMesh->InitialiseElement(3, *Connectivity3, 2);
    myMesh->InitialiseElement(4, *Connectivity4, 2);


    FE_Solution* FE = new FE_Solution(*myMesh, 1);

    SparseMatrix* A = new SparseMatrix(*FE, numElements);
    Vector* F = new Vector (FE->GetNumberOfDofs());

    int n_q = 9;

    int numElementDofs;

    for (int k=1; k<=myMesh->GetNumElements(); k++)
    {
        Vector* quadratureWeights = new Vector (n_q);
        Matrix* localQuadraturePoints = new Matrix (myMesh->GetDimension(),n_q);
        Matrix* globalQuadraturePoints = new Matrix (myMesh->GetDimension(),n_q);

        numElementDofs = FE->GetNumElementDofs(k);
        Matrix* A_loc = new Matrix (numElementDofs, numElementDofs);
        Vector* f_loc = new Vector (numElementDofs);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);

        for (int q=1; q<=n_q; q++)
        {
            Vector* basisValues = new Vector(numElementDofs);
            Matrix* basisGrad = new Matrix(myMesh->GetDimension(),numElementDofs);
            Matrix* jacobian = new Matrix (myMesh->GetDimension(), myMesh->GetDimension());
            Vector* pointToEval = new Vector(localQuadraturePoints->GetColumnAsVector(q));

            myMesh->GetElement(k)->ComputeMappingJacobian(*pointToEval, *jacobian);
            delete pointToEval;

            Vector* localQuadPoints = new Vector (localQuadraturePoints->GetColumnAsVector(q));

            FE->ComputeBasis(k, *localQuadPoints, *basisValues);
            FE->ComputeGradBasis(k, *localQuadPoints, *basisGrad);

            delete localQuadPoints;

            for (int i=1; i<=A_loc->GetNumberOfRows(); i++)
            {
                for (int j=1; j<=A_loc->GetNumberOfColumns(); j++)
                {
                    (*A_loc)(i,j) += ((*basisGrad)(1,i)*(*basisGrad)(1,j) + (*basisGrad)(2,i)*(*basisGrad)(2,j))*(*quadratureWeights)(q);
                }
            }
            for (int i=1; i<=f_loc->GetSize(); i++)
            {
                (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q);
            }

            delete basisValues;
            delete basisGrad;
            delete jacobian;
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;

        for (int i=1; i<=A_loc->GetNumberOfRows(); i++)
        {
            for (int j=1; j<=A_loc->GetNumberOfColumns(); j++)
            {
                A->AddValue(A->Read((FE->GetElementDofs(k))(i), (FE->GetElementDofs(k))(j)) + (*A_loc)(i,j), (FE->GetElementDofs(k))(i), (FE->GetElementDofs(k))(j));
            }
        }
        for (int i=1; i<=f_loc->GetSize(); i++)
        {
            (*F)((FE->GetElementDofs(k))(i)) += (*f_loc)(i);
        }

        delete A_loc;
        delete f_loc;
    }

    for (int i=1; i<=A->GetNumberOfRows(); i++)
    {
        for (int j=1; j<=A->GetNumberOfColumns(); j++)
        {
            if (i!=5 || j!=5)
            {
                if (i==j)
                {
                    A->AddValue(1, i, j);
                }
                else
                {
                    A->AddValue(0, i, j);
                }
            }
        }
    }
    for (int i=1; i<=F->GetSize(); i++)
    {
        if (i!=5)
        {
            (*F)(i) = 0;
        }
    }
    std::cout << "A = \n" << *A;
    std::cout << "F = \n" << *F;

    A->CGSolveSystem(*F, FE->GetSolutionVector(), 1e-9, 1000);
    std::cout << "U = \n" << FE->GetSolutionVector();
}

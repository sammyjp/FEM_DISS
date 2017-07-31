#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

void Example1D();
void Example2D();


int main(int argc, char* argv[])
{
    Example1D();
    return 0;
}

void Example1D()
{
    int numElements = 4;
    int polynomialDegree = 1;

    Matrix* Grid = new Matrix (1,numElements+1);

    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = (j-1)*(1.0/numElements);
    }

    Matrix* Connectivity = new Matrix (numElements, 2);

    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        for (int j=1; j<=Connectivity->GetNumberOfColumns(); j++)
        {
            (*Connectivity)(i,j) = i + (j-1);
        }
    }

    Mesh* myMesh = new Mesh(*Grid, numElements);

    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (Connectivity->GetRowAsVector(i));
        myMesh->InitialiseElement(i, *Con_temp, 0);
        delete Con_temp;
    }

    FE_Solution* FE = new FE_Solution(*myMesh, polynomialDegree);

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
    int numElements = 16;

    Matrix* Grid = new Matrix (2,25);
    double value=0;
    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = ((j-1)%((int)(sqrt(numElements))+1))*(1.0/sqrt(numElements));
        if (((j-1)%((int)(sqrt(numElements))+1)) == 0 && j>1)
        {
            value += (1.0/sqrt(numElements));
        }
        (*Grid)(2,j) = value;
    }

    Matrix* Connectivity = new Matrix (numElements, 4);

    int k=0;
    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        k++;
        if ((i%((int)(sqrt(numElements))+1)) == 0)
        {
            k++;
        }
        for (int j=1; j<=Connectivity->GetNumberOfColumns(); j++)
        {
            if (j==1 || j==2)
            {
                (*Connectivity)(i,j) = k + (j-1);
            }
            if (j==3 || j==4)
            {
                (*Connectivity)(i,j) = k + (int)(sqrt(numElements)) + 4 - (j-1);
            }
        }
    }

    Mesh* myMesh = new Mesh(*Grid, numElements);

    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (Connectivity->GetRowAsVector(i));
        myMesh->InitialiseElement(i, *Con_temp, 2);
        delete Con_temp;
    }

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
            if (i<=5 || j<=1 || i>=21 || j<=21 || (i%5)==1 || (j%5)==1 || (i%5)==0 || (j%5)==0)
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
        if (i<=5 || i>=21 || (i%5)==1 || (i%5)==0)
        {
            (*F)(i) = 0;
        }
    }
    std::cout << "A = \n" << *A;
    std::cout << "F = \n" << *F;

    A->CGSolveSystem(*F, FE->GetSolutionVector(), 1e-9, 1000);
    std::cout << "U = \n" << FE->GetSolutionVector();
}

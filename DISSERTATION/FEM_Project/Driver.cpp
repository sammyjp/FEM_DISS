#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

void Example1D();
void ErrorAnalysis1D();
void Example2D();

const double Pi = 4*atan(1);

int main(int argc, char* argv[])
{
    Example1D();
    ErrorAnalysis1D();
    return 0;
}

void Example1D()
{
    /*
        f(x) = -2*pi^2*cos(2*pi*x)
    */

    double a = 0;
    double b = 1;
    int numElements = 15;
    int polynomialDegree = 1;

    Matrix* Grid = new Matrix (1,numElements+1);

    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = (j-1)*((b-a)/numElements) + a;
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

    delete Grid;

    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (Connectivity->GetRowAsVector(i));
        myMesh->InitialiseElement(i, *Con_temp, 0);
        delete Con_temp;
    }

    delete Connectivity;

    FE_Solution* FE = new FE_Solution(*myMesh, polynomialDegree);

    SparseMatrix* A = new SparseMatrix(*FE, numElements);
    Vector* F = new Vector (FE->GetNumberOfDofs());

    int n_q = 20;

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
                (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q) * (-2.0*pow(Pi,2.0)*cos(2.0*Pi*(*globalQuadraturePoints)(1,q)));
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

    A->CGSolveSystem(*F, FE->GetSolutionVector(), 1e-9, 1000);

    delete A;
    delete F;

    std::cout << std::setprecision(5) << std::scientific;
    std::cout << std::setw(8) << "U" << std::endl << std::setw(14) << FE->GetSolutionVector();

    // Creates an output file
    std::ofstream output_file;
    output_file.setf(std::ios::scientific,std::ios::floatfield);
    output_file.precision(6);

    // Open file (and perform a check)
    output_file.open("1DExample.dat");
    assert(output_file.is_open());
    for (int i=1; i<=myMesh->GetDimension(); i++)
    {
        std::cout << std::setw(8) << "x" << i << std::setw(13);
    }
    std::cout << "u" << std::endl;
    for (int j=1; j<=myMesh->GetNumNodes(); j++)
    {
        for (int i=1; i<=myMesh->GetDimension(); i++)
        {
            std::cout << std::setw(14) << (myMesh->GetAllGridPoints())(i,j) << std::setw(14);
            output_file << (myMesh->GetAllGridPoints())(i,j) << " ";
        }
        std::cout << (FE->GetSolutionVector())(j) << std::endl;
        output_file << (FE->GetSolutionVector())(j) << std::endl;
    }

    output_file.close();

    double Error_L2 = 0;

    for (int k=1; k<=myMesh->GetNumElements(); k++)
    {
        Vector* quadratureWeights = new Vector (n_q);
        Matrix* localQuadraturePoints = new Matrix (1,n_q);
        Matrix* globalQuadraturePoints = new Matrix (1,n_q);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
        for (int q=1; q<=n_q; q++)
        {
            Error_L2 += (pow((pow(sin(Pi*(*globalQuadraturePoints)(1,q)),2.0)) - (FE->ComputeUh(k, (*localQuadraturePoints)(1,q))),2.0) * (*quadratureWeights)(q));
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;
    }

    Error_L2 = sqrt(Error_L2);

    std::cout << std::endl;
    std::cout << std::setw(13) << "||u-uh|| = " << Error_L2 << std::endl;
    std::cout << std::endl;

    delete FE;
    delete myMesh;
}

void ErrorAnalysis1D()
{
    /*
        f(x) = -2*pi^2*cos(2*pi*x)
    */

    double a = 0;
    double b = 1;
    int numElements = 15;

    std::cout << "h = " << (b-a)/numElements << std::endl;

    Matrix* Grid = new Matrix (1,numElements+1);

    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = (j-1)*((b-a)/numElements) + a;
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

    delete Grid;

    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (Connectivity->GetRowAsVector(i));
        myMesh->InitialiseElement(i, *Con_temp, 0);
        delete Con_temp;
    }

    delete Connectivity;

    // Creates an output file
    std::ofstream error_output;
    error_output.setf(std::ios::scientific,std::ios::floatfield);
    error_output.precision(6);

    // Open file (and perform a check)
    error_output.open("1DError.dat");
    assert(error_output.is_open());

    std::cout << std::setprecision(5) << std::scientific;
    std::cout << std::setw(15) << "NumDofs^(1/d)" << std::setw(12) << "||u-uh||" << std::endl;

    for (int polynomialDegree=1; polynomialDegree<=10; polynomialDegree++)
    {
        FE_Solution* FE = new FE_Solution(*myMesh, polynomialDegree);

        SparseMatrix* A = new SparseMatrix(*FE, numElements);
        Vector* F = new Vector (FE->GetNumberOfDofs());

        int n_q = 20;

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
                    (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q) * (-2.0*pow(Pi,2.0)*cos(2.0*Pi*(*globalQuadraturePoints)(1,q)));
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

        A->CGSolveSystem(*F, FE->GetSolutionVector(), 1e-9, 1000);

        delete A;
        delete F;

        double Error_L2 = 0;

        for (int k=1; k<=myMesh->GetNumElements(); k++)
        {
            Vector* quadratureWeights = new Vector (n_q);
            Matrix* localQuadraturePoints = new Matrix (1,n_q);
            Matrix* globalQuadraturePoints = new Matrix (1,n_q);

            myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
            for (int q=1; q<=n_q; q++)
            {
                Error_L2 += (pow((pow(sin(Pi*(*globalQuadraturePoints)(1,q)),2.0)) - (FE->ComputeUh(k, (*localQuadraturePoints)(1,q))),2.0) * (*quadratureWeights)(q));
            }

            delete quadratureWeights;
            delete localQuadraturePoints;
            delete globalQuadraturePoints;
        }

        Error_L2 = sqrt(Error_L2);

        std::cout << std::setw(14) << pow(FE->GetNumberOfDofs(), 1.0/(myMesh->GetDimension())) << std::setw(15) << Error_L2 << std::endl;
        error_output << pow(FE->GetNumberOfDofs(), 1.0/(myMesh->GetDimension())) << " " << Error_L2 << std::endl;

        delete FE;
    }

    error_output.close();

    delete myMesh;
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

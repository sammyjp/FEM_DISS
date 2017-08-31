#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <cassert>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

void Example1D();
void ErrorAnalysis1D();
void Example2DQuads();
void Example2DTriangles();
void Example2DHybrid();

const double Pi = 4*atan(1);

int main(int argc, char* argv[])
{
    Example2DHybrid();
    return 0;
}

void Example1D()
{
    /*
        f(x) = -2*pi^2*cos(2*pi*x)
    */

    double a = 0;
    double b = 1;
    int numElements = 12;
    int polynomialDegree = 2;

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

    int n_q = 5;

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

    for (int polynomialDegree=1; polynomialDegree<=4; polynomialDegree++)
    {
        // Creates an output file
        std::ofstream error_output;
        error_output.setf(std::ios::scientific,std::ios::floatfield);
        error_output.precision(6);

        std::stringstream fileName;
        fileName << polynomialDegree;

        // Open file (and perform a check)
        error_output.open("1DError_p"+fileName.str()+".dat");
        assert(error_output.is_open());

        std::cout << std::setprecision(5) << std::scientific;
        std::cout << std::setw(15) << "NumDofs^(1/d)" << std::setw(12) << "||u-uh||" << std::endl;

        for (int numElements=4; numElements<=64; numElements=numElements*2)
        {
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
            delete myMesh;
        }
        std::cout << std::endl;
        error_output.close();
    }
}

void Example2DQuads()
{
    /*
        f(x) = (1 + 2 * Pi^2) * sin(Pi*x) * sin(Pi*y)
    */

    double h=1.0/8.0;
    double tempNumElements = pow(1.0/h,2.0);
    double tempNumGridPoints = pow(1.0/h + 1.0,2.0);
    int numRowPoints = (int)(sqrt(tempNumGridPoints));
    int numElements = (int)(tempNumElements);
    int numGridPoints = (int)(tempNumGridPoints);

    Matrix* Grid = new Matrix (2, numGridPoints);
    double value=0;
    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = ((j-1)%((int)(sqrt(numElements))+1))*(h);
        if (((j-1)%((int)(sqrt(numElements))+1)) == 0 && j>1)
        {
            value += h;
        }
        (*Grid)(2,j) = value;
    }

    Matrix* Connectivity = new Matrix (numElements, 4);

    int k=0;
    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        k++;
        if ((i%numRowPoints) == 0)
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

    int n_q = 64;

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
                    (*A_loc)(i,j) += (((*basisGrad)(1,i)*(*basisGrad)(1,j) + (*basisGrad)(2,i)*(*basisGrad)(2,j)) + (*basisValues)(i)*(*basisValues)(j))*(*quadratureWeights)(q);
                }
            }

            for (int i=1; i<=f_loc->GetSize(); i++)
            {
                (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q) * ((1+2*pow(Pi,2.0))*sin(Pi*(*globalQuadraturePoints)(1,q))*sin(Pi*(*globalQuadraturePoints)(2,q)));
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
            if (i<=(int)(sqrt(numGridPoints)) || j<=(int)(sqrt(numGridPoints)) || i>=numGridPoints-(int)(sqrt(numGridPoints))+1 || j>=numGridPoints-(int)(sqrt(numGridPoints))+1 || (i%(int)(sqrt(numGridPoints)))==1 || (j%(int)(sqrt(numGridPoints)))==1 || (i%(int)(sqrt(numGridPoints)))==0 || (j%(int)(sqrt(numGridPoints)))==0)
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
        if (i<=(int)(sqrt(numGridPoints)) || i>=numGridPoints-(int)(sqrt(numGridPoints))+1 || (i%(int)(sqrt(numGridPoints)))==1 || (i%(int)(sqrt(numGridPoints)))==0)
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
    output_file.open("2DQuadsExample.dat");
    assert(output_file.is_open());
    std::cout << std::setw(8);
    for (int i=1; i<=myMesh->GetDimension(); i++)
    {
        std::cout << "x" << i << std::setw(13);
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
        Matrix* localQuadraturePoints = new Matrix (2,n_q);
        Matrix* globalQuadraturePoints = new Matrix (2,n_q);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
        for (int q=1; q<=n_q; q++)
        {
            Vector* localQuadPoints = new Vector (localQuadraturePoints->GetColumnAsVector(q));
            Error_L2 += (pow((sin(Pi*(*globalQuadraturePoints)(1,q))*sin(Pi*(*globalQuadraturePoints)(2,q))) - (FE->ComputeUh(k, (*localQuadPoints))),2.0) * (*quadratureWeights)(q));
            delete localQuadPoints;
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;
    }

    Error_L2 = sqrt(fabs(Error_L2));

    std::cout << std::endl;
    std::cout << std::setw(13) << "||u-uh|| = " << Error_L2 << std::endl;
    std::cout << std::endl;

    delete FE;
    delete myMesh;
}

void Example2DTriangles()
{
    /*
        f(x) = (1 + 2 * Pi^2) * sin(Pi*x) * sin(Pi*y)
    */

    double h=1.0/8.0;
    double tempNumElements = pow(1.0/h,2.0);
    double tempNumGridPoints = pow(1.0/h + 1.0,2.0);
    int numRowPoints = (int)(sqrt(tempNumGridPoints));
    int numGridPoints = (int)(tempNumGridPoints);
    int numElements = (int)(tempNumElements)*2;

    Matrix* Grid = new Matrix (2, numGridPoints);
    double value=0;
    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = (j-1)%(numRowPoints)*(h);
        if ((j-1)%(numRowPoints) == 0 && j>1)
        {
            value += h;
        }
        (*Grid)(2,j) = value;
    }

    Matrix* Connectivity = new Matrix (numElements, 3);

    int k=1;
    int counter = 1;
    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        if ((i%2)==0)
        {
            k++;
        }
        if (((i)%((numRowPoints-1)*2)) == 1 && i!=1)
        {
            k=counter*numRowPoints + 1;
            counter++;
        }
        for (int j=1; j<=Connectivity->GetNumberOfColumns(); j++)
        {
            if ((i%2)!=0)
            {
                if (j==1 || j==2)
                {
                    (*Connectivity)(i,j) = k + (j-1);
                }
                if (j==3)
                {
                    (*Connectivity)(i,j) = k + numRowPoints + 2 - (j-1);
                }
            }
            else
            {
                if (j==1)
                {
                    (*Connectivity)(i,j) = k;
                }
                if (j==2 || j==3)
                {
                    (*Connectivity)(i,j) = k + numRowPoints + 1 - (j-1);
                }
            }
        }
    }

    Mesh* myMesh = new Mesh(*Grid, numElements);

    for (int i=1; i<=Connectivity->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (Connectivity->GetRowAsVector(i));
        myMesh->InitialiseElement(i, *Con_temp, 1);
        delete Con_temp;
    }

    FE_Solution* FE = new FE_Solution(*myMesh, 1);

    SparseMatrix* A = new SparseMatrix(*FE, numElements);

    Vector* F = new Vector (FE->GetNumberOfDofs());

    int n_q = 16;

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
                    (*A_loc)(i,j) += (((*basisGrad)(1,i)*(*basisGrad)(1,j) + (*basisGrad)(2,i)*(*basisGrad)(2,j)) + (*basisValues)(i)*(*basisValues)(j))*(*quadratureWeights)(q);
                }
            }

            for (int i=1; i<=f_loc->GetSize(); i++)
            {
                (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q) * ((1+2*pow(Pi,2.0))*sin(Pi*(*globalQuadraturePoints)(1,q))*sin(Pi*(*globalQuadraturePoints)(2,q)));
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
            if (i<=(int)(sqrt(numGridPoints)) || j<=(int)(sqrt(numGridPoints)) || i>=numGridPoints-(int)(sqrt(numGridPoints))+1 || j>=numGridPoints-(int)(sqrt(numGridPoints))+1 || (i%(int)(sqrt(numGridPoints)))==1 || (j%(int)(sqrt(numGridPoints)))==1 || (i%(int)(sqrt(numGridPoints)))==0 || (j%(int)(sqrt(numGridPoints)))==0)
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
        if (i<=(int)(sqrt(numGridPoints)) || i>=numGridPoints-(int)(sqrt(numGridPoints))+1 || (i%(int)(sqrt(numGridPoints)))==1 || (i%(int)(sqrt(numGridPoints)))==0)
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
    output_file.open("2DTrianglesExample.dat");
    assert(output_file.is_open());
    std::cout << std::setw(8);
    for (int i=1; i<=myMesh->GetDimension(); i++)
    {
        std::cout << "x" << i << std::setw(13);
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
        Matrix* localQuadraturePoints = new Matrix (2,n_q);
        Matrix* globalQuadraturePoints = new Matrix (2,n_q);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
        for (int q=1; q<=n_q; q++)
        {
            Vector* localQuadPoints = new Vector (localQuadraturePoints->GetColumnAsVector(q));
            Error_L2 += (pow((sin(Pi*(*globalQuadraturePoints)(1,q))*sin(Pi*(*globalQuadraturePoints)(2,q))) - (FE->ComputeUh(k, (*localQuadPoints))),2.0) * (*quadratureWeights)(q));
            delete localQuadPoints;
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;
    }

    Error_L2 = sqrt(fabs(Error_L2));

    std::cout << std::endl;
    std::cout << std::setw(13) << "||u-uh|| = " << Error_L2 << std::endl;
    std::cout << std::endl;

    delete FE;
    delete myMesh;
}

void Example2DHybrid()
{
    /*
        f(x) = (1 + 2 * Pi^2) * sin(Pi*x) * sin(Pi*y)
    */

    double h=1.0/8.0; // Must be greater than or equal to 0.25
    double tempNumElements = pow(1.0/h,2.0);
    double tempNumGridPoints = pow(1.0/h + 1.0,2.0);
    int numRowPoints = (int)(sqrt(tempNumGridPoints));
    int numGridPoints = (int)(tempNumGridPoints);
    int numTriangleElements;
    int numQuadElements;
    if ((int)(tempNumElements)%2==0)
    {
        numTriangleElements = (int)(tempNumElements);
        numQuadElements = (int)(tempNumElements)/2;
    }
    else
    {
        numTriangleElements = (numRowPoints-2)*(numRowPoints-1) + (numRowPoints-1)*2;
        numQuadElements = (numRowPoints-2)/2*(numRowPoints-1);
    }
    int numElements = numTriangleElements + numQuadElements;

    Matrix* Grid = new Matrix (2, numGridPoints);
    double value=0;
    for (int j=1; j<=Grid->GetNumberOfColumns(); j++)
    {
        (*Grid)(1,j) = (j-1)%(numRowPoints)*(h);
        if ((j-1)%(numRowPoints) == 0 && j>1)
        {
            value += h;
        }
        (*Grid)(2,j) = value;
    }

    Matrix* ConnectivityTriangles = new Matrix (numTriangleElements, 3);
    Matrix* ConnectivityQuads = new Matrix (numQuadElements, 4);

    if (numTriangleElements==numQuadElements*2)
    {
        int k=(numRowPoints-1)/2;
        for (int i=1; i<=ConnectivityQuads->GetNumberOfRows(); i++)
        {
            k++;
            if ((i%((numRowPoints-1)/2)) == 1 && i!=1)
            {
                k=k+(numRowPoints-1)/2+1;
            }
            for (int j=1; j<=ConnectivityQuads->GetNumberOfColumns(); j++)
            {
                if (j==1 || j==2)
                {
                    (*ConnectivityQuads)(i,j) = k + (j-1);
                }
                if (j==3 || j==4)
                {
                    (*ConnectivityQuads)(i,j) = k + (numRowPoints-1) + 4 - (j-1);
                }
            }
        }

        k=1;
        int counter = 1;
        for (int i=1; i<=ConnectivityTriangles->GetNumberOfRows(); i++)
        {
            if ((i%2)==0)
            {
                k++;
            }
            if (((i)%((numRowPoints-1))) == 1 && i!=1)
            {
                k=counter*numRowPoints + 1;
                counter++;
            }
            for (int j=1; j<=ConnectivityTriangles->GetNumberOfColumns(); j++)
            {
                if ((i%2)!=0)
                {
                    if (j==1 || j==2)
                    {
                        (*ConnectivityTriangles)(i,j) = k + (j-1);
                    }
                    if (j==3)
                    {
                        (*ConnectivityTriangles)(i,j) = k + numRowPoints + 2 - (j-1);
                    }
                }
                else
                {
                    if (j==1)
                    {
                        (*ConnectivityTriangles)(i,j) = k;
                    }
                    if (j==2 || j==3)
                    {
                        (*ConnectivityTriangles)(i,j) = k + numRowPoints + 1 - (j-1);
                    }
                }
            }
        }
    }
    else
    {
        int k=(numRowPoints-2)/2 + 1;
        for (int i=1; i<=ConnectivityQuads->GetNumberOfRows(); i++)
        {
            k++;
            if ((i%((numRowPoints-2)/2)) == 1 && i!=1)
            {
                k=k+(numRowPoints-2)/2+2;
            }
            for (int j=1; j<=ConnectivityQuads->GetNumberOfColumns(); j++)
            {
                if (j==1 || j==2)
                {
                    (*ConnectivityQuads)(i,j) = k + (j-1);
                }
                if (j==3 || j==4)
                {
                    (*ConnectivityQuads)(i,j) = k + (numRowPoints-1) + 4 - (j-1);
                }
            }
        }
        k=1;
        int counter = 1;
        for (int i=1; i<=ConnectivityTriangles->GetNumberOfRows(); i++)
        {
            if ((i%2)==0)
            {
                k++;
            }
            if (((i)%((numRowPoints))) == 1 && i!=1)
            {
                k=counter*numRowPoints + 1;
                counter++;
            }
            for (int j=1; j<=ConnectivityTriangles->GetNumberOfColumns(); j++)
            {
                if ((i%2)!=0)
                {
                    if (j==1 || j==2)
                    {
                        (*ConnectivityTriangles)(i,j) = k + (j-1);
                    }
                    if (j==3)
                    {
                        (*ConnectivityTriangles)(i,j) = k + numRowPoints + 2 - (j-1);
                    }
                }
                else
                {
                    if (j==1)
                    {
                        (*ConnectivityTriangles)(i,j) = k;
                    }
                    if (j==2 || j==3)
                    {
                        (*ConnectivityTriangles)(i,j) = k + numRowPoints + 1 - (j-1);
                    }
                }
            }
        }
    }

    Mesh* myMesh = new Mesh(*Grid, numElements);

    int elementCounter = 1;
    for (int i=1; i<=ConnectivityTriangles->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (ConnectivityTriangles->GetRowAsVector(i));
        myMesh->InitialiseElement(elementCounter, *Con_temp, 1);
        delete Con_temp;
        elementCounter++;
    }
    for (int i=1; i<=ConnectivityQuads->GetNumberOfRows(); i++)
    {
        Vector* Con_temp = new Vector (ConnectivityQuads->GetRowAsVector(i));
        myMesh->InitialiseElement(elementCounter, *Con_temp, 2);
        delete Con_temp;
        elementCounter++;
    }

    FE_Solution* FE = new FE_Solution(*myMesh, 1);

    SparseMatrix* A = new SparseMatrix(*FE, numElements);

    Vector* F = new Vector (FE->GetNumberOfDofs());

    int n_q = 16;

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
                    (*A_loc)(i,j) += (((*basisGrad)(1,i)*(*basisGrad)(1,j) + (*basisGrad)(2,i)*(*basisGrad)(2,j)) + (*basisValues)(i)*(*basisValues)(j))*(*quadratureWeights)(q);
                }
            }

            for (int i=1; i<=f_loc->GetSize(); i++)
            {
                (*f_loc)(i) += (*basisValues)(i)*(*quadratureWeights)(q) * ((1+2*pow(Pi,2.0))*sin(Pi*(*globalQuadraturePoints)(1,q))*sin(Pi*(*globalQuadraturePoints)(2,q)));
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
            if (i<=(int)(sqrt(numGridPoints)) || j<=(int)(sqrt(numGridPoints)) || i>=numGridPoints-(int)(sqrt(numGridPoints))+1 || j>=numGridPoints-(int)(sqrt(numGridPoints))+1 || (i%(int)(sqrt(numGridPoints)))==1 || (j%(int)(sqrt(numGridPoints)))==1 || (i%(int)(sqrt(numGridPoints)))==0 || (j%(int)(sqrt(numGridPoints)))==0)
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
        if (i<=(int)(sqrt(numGridPoints)) || i>=numGridPoints-(int)(sqrt(numGridPoints))+1 || (i%(int)(sqrt(numGridPoints)))==1 || (i%(int)(sqrt(numGridPoints)))==0)
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
    output_file.open("2DHybridExample.dat");
    assert(output_file.is_open());
    std::cout << std::setw(8);
    for (int i=1; i<=myMesh->GetDimension(); i++)
    {
        std::cout << "x" << i << std::setw(13);
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
        Matrix* localQuadraturePoints = new Matrix (2,n_q);
        Matrix* globalQuadraturePoints = new Matrix (2,n_q);

        myMesh->GetElement(k)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
        for (int q=1; q<=n_q; q++)
        {
            Vector* localQuadPoints = new Vector (localQuadraturePoints->GetColumnAsVector(q));
            Error_L2 += (pow((sin(Pi*(*globalQuadraturePoints)(1,q))*sin(Pi*(*globalQuadraturePoints)(2,q))) - (FE->ComputeUh(k, (*localQuadPoints))),2.0) * (*quadratureWeights)(q));
            delete localQuadPoints;
        }

        delete quadratureWeights;
        delete localQuadraturePoints;
        delete globalQuadraturePoints;
    }

    Error_L2 = sqrt(fabs(Error_L2));

    std::cout << std::endl;
    std::cout << std::setw(13) << "||u-uh|| = " << Error_L2 << std::endl;
    std::cout << std::endl;

    delete FE;
    delete myMesh;
}


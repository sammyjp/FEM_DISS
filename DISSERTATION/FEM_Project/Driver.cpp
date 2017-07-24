#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

int main(int argc, char* argv[])
{

    //TODO create polynomial (deal 2 fe) class
    //TODO in polynomial class have the basis shape functions with methods to evaluate on functions
    //TODO add DOFs into FE_Solution class

    /*
        for element 1:numElements
            get quadrature;

            for idofIndex 1:dofNumbersPerElement.GetLength
                for jdofIndex 1:dofNumbersPerElement.GetLength
                    calculate basis and/or gradbasis values;
                    localMatrix(i,j) = performQuadrature on element with values;
                end
                calculate values;
                rhs(i) = performQuadrature on element with values;
            end

            add to global matrix;
            reset;
        end

        solve system;
    */

    int numElements = 8;

    Matrix* Grid = new Matrix (1,9);
    (*Grid)(1,1) = 0;
    (*Grid)(1,2) = 0.125;
    (*Grid)(1,3) = 0.25;
    (*Grid)(1,4) = 0.375;
    (*Grid)(1,5) = 0.5;
    (*Grid)(1,5) = 0.625;
    (*Grid)(1,5) = 0.75;
    (*Grid)(1,5) = 0.875;
    (*Grid)(1,5) = 1;

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


    FE_Solution* FE = new FE_Solution(*myMesh);

    int n_q = 10;
    Vector* quadratureWeights = new Vector (n_q);
    Matrix* localQuadraturePoints = new Matrix (1,n_q);
    Matrix* globalQuadraturePoints = new Matrix (1,n_q);
    double I;
    double F;

    Vector* functionTemp = new Vector(n_q);
    Vector* functionPoints = new Vector (n_q);

    Matrix* globalMatrix = new Matrix(9,9);
    Vector* globalF = new Vector(9);

    for (int i=1; i<=numElements; i++)
    {
        Matrix* elementStiffness = new Matrix(2,2);
        Vector* elementF = new Vector(2);

        myMesh->GetElement(i)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);
        for (int k=1; k<=myMesh->GetElement(i)->GetElementConnectivityArray().GetSize(); k++)
        {
            for (int l=1; l<=myMesh->GetElement(i)->GetElementConnectivityArray().GetSize(); l++)
            {
                FE->ComputeLinearBasisFunctionDerivativeValues((myMesh->GetElement(i)->GetElementConnectivityArray())(k), *functionTemp, *globalQuadraturePoints);
                FE->ComputeLinearBasisFunctionDerivativeValues((myMesh->GetElement(i)->GetElementConnectivityArray())(l), *functionPoints, *globalQuadraturePoints);
                for (int j=0; j<functionPoints->GetSize();j++)
                {
                    (*functionPoints)[j] = (*functionPoints)[j]*(*functionTemp)[j];
                }
                I = myMesh->GetElement(i)->PerformElementQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *functionPoints);
                (*elementStiffness)(k,l) = I;
            }

            FE->ComputeLinearBasisFunctionValues((myMesh->GetElement(i)->GetElementConnectivityArray())(k), *functionTemp, *globalQuadraturePoints);
            for (int j=0; j<functionTemp->GetSize();j++)
            {
                (*functionTemp)[j] = (*functionTemp)[j]*(*globalQuadraturePoints)(1,j+1);
            }
            F = myMesh->GetElement(i)->PerformElementQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *functionTemp);
            (*elementF)(k) = F;
        }

        for (int g=1; g<=elementStiffness->GetNumberOfRows(); g++)
        {
            for (int h=1; h<=elementStiffness->GetNumberOfRows(); h++)
            {
                (*globalMatrix)((myMesh->GetElement(i)->GetElementConnectivityArray())(g),(myMesh->GetElement(i)->GetElementConnectivityArray())(h)) = (*elementStiffness)(g,h);
            }
            (*globalF)((myMesh->GetElement(i)->GetElementConnectivityArray())(g)) = (*elementF)(g);
        }


        delete elementStiffness;
        delete elementF;
    }
    std::cout << *globalMatrix;

    std::cout << *globalF;



    return 0;
}

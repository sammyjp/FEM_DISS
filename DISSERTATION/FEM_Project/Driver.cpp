#include <iostream>
#include <cmath>

#include "SP_FEM.hpp"
#include "FE_Solution.hpp"
#include "SparseMatrix.hpp"

int main(int argc, char* argv[])
{
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


    FE_Solution* FE = new FE_Solution(*myMesh, 1);

    int n_q = 4;
    Vector* quadratureWeights = new Vector (n_q);
    Matrix* localQuadraturePoints = new Matrix (1,n_q);
    Matrix* globalQuadraturePoints = new Matrix (1,n_q);
    double I;

    myMesh->GetElement(1)->GetQuadrature(n_q, *quadratureWeights, *localQuadraturePoints, *globalQuadraturePoints);

    for (int i=1; i<=n_q; i++)
    {
        Vector* basisValues = new Vector(FE->GetElementPolynomialSpace(1)->GetNumElementDofs());
        Matrix* basisGrad = new Matrix(1,FE->GetElementPolynomialSpace(1)->GetNumElementDofs());

        FE->GetElementPolynomialSpace(1)->ComputeBasis((*localQuadraturePoints)(1,i), *basisValues);
        FE->GetElementPolynomialSpace(1)->ComputeGradBasis((*localQuadraturePoints)(1,i), *basisGrad);

        std::cout << *basisValues;
        std::cout << *basisGrad << std::endl;

        delete basisValues;
        delete basisGrad;
    }


    return 0;
}

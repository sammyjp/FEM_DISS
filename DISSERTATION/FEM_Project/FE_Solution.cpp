#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "FE_Solution.hpp"
#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Vector.hpp"

FE_Solution::FE_Solution(Mesh& mesh, int polynomialDegree)
{
    mMeshReference = &mesh;
    mPolynomialDegree = polynomialDegree;

    InitialiseElementDofs();
    mNumDofs = (*dofStart)(dofStart->GetSize()) - 1;

    mSolutionVector = new Vector (mNumDofs);
}

FE_Solution::~FE_Solution()
{
    delete mSolutionVector;
    delete dofStart;

    for (int i=0; i<mMeshReference->GetNumElements(); i++)
    {
        delete mElementPolySpaceArray[i];
    }
    delete[] mElementPolySpaceArray;
}

void FE_Solution::InitialiseElementDofs()
{
    mElementPolySpaceArray = new PolynomialSpace* [mMeshReference->GetNumElements()];
    dofStart = new Vector(mMeshReference->GetNumElements() + 1);

    (*dofStart)[0] = mMeshReference->GetNumNodes() + 1;

    for (int i=0; i<mMeshReference->GetNumElements(); i++)
    {
        mElementPolySpaceArray[i] = new PolynomialSpace(mPolynomialDegree, mMeshReference->GetElement(i+1)->GetElementType());

        if (mMeshReference->GetElement(i+1)->GetElementType() == 0)
        {
            (*dofStart)[i+1] = (*dofStart)[i] + (mPolynomialDegree - 1);
        }
        else if (mMeshReference->GetElement(i+1)->GetElementType() == 1)
        {

        }
        else if (mMeshReference->GetElement(i+1)->GetElementType() == 2)
        {
            (*dofStart)[i+1] = (*dofStart)[i] + ((mPolynomialDegree-1)*4);
        }
    }
}

PolynomialSpace* FE_Solution::GetElementPolynomialSpace(int elementNumber) const
{
    assert(elementNumber >= 1 && elementNumber <= mMeshReference->GetNumElements());
    return mElementPolySpaceArray[elementNumber - 1];
}

Vector FE_Solution::GetElementDofs(int elementNumber)
{
    int connectivitySize = mMeshReference->GetElement(elementNumber)->GetElementConnectivityArray().GetSize();
    Vector elementDofs(connectivitySize + (*dofStart)(elementNumber+1) - (*dofStart)(elementNumber));

    for (int i=0; i<connectivitySize; i++)
    {
        elementDofs[i] = (mMeshReference->GetElement(elementNumber)->GetElementConnectivityArray())[i];
    }
    for (int i=connectivitySize; i<elementDofs.GetSize(); i++)
    {
        elementDofs[i] = (*dofStart)(elementNumber) + (i-connectivitySize);
    }

    return elementDofs;
}

int FE_Solution::GetNumberOfDofs()
{
    return mNumDofs;
}

int FE_Solution::GetNumElementDofs(int elementNumber)
{
    return GetElementPolynomialSpace(elementNumber)->GetNumElementDofs();
}

Vector& FE_Solution::GetSolutionVector()
{
    return *mSolutionVector;
}

void FE_Solution::ComputeBasis(int elementNumber, double localGridPoint, Vector& basisValues)
{
    GetElementPolynomialSpace(elementNumber)->ComputeBasis(localGridPoint, basisValues);
}

void FE_Solution::ComputeBasis(int elementNumber, Vector& localGridPoint, Vector& basisValues)
{
    GetElementPolynomialSpace(elementNumber)->ComputeBasis(localGridPoint, basisValues);
}

void FE_Solution::ComputeGradBasis(int elementNumber, double localGridPoint, Matrix& gradBasisValues)
{
    GetElementPolynomialSpace(elementNumber)->ComputeGradBasis(localGridPoint, gradBasisValues);

    Matrix* jacobian = new Matrix (mMeshReference->GetDimension(), mMeshReference->GetDimension());

    mMeshReference->GetElement(elementNumber)->ComputeMappingJacobian(localGridPoint, *jacobian);

    gradBasisValues = gradBasisValues*(1.0/(jacobian->CalculateDeterminant()));

    delete jacobian;
}

void FE_Solution::ComputeGradBasis(int elementNumber, Vector& localGridPoint, Matrix& gradBasisValues)
{
    GetElementPolynomialSpace(elementNumber)->ComputeGradBasis(localGridPoint, gradBasisValues);

    Matrix* jacobian = new Matrix (mMeshReference->GetDimension(), mMeshReference->GetDimension());

    mMeshReference->GetElement(elementNumber)->ComputeMappingJacobian(localGridPoint, *jacobian);

    gradBasisValues = gradBasisValues*(1.0/(jacobian->CalculateDeterminant()));

    delete jacobian;
}

double FE_Solution::ComputeUh(int elementNumber, double localGridPoint)
{
    double Uh = 0;

    Vector* basisValues = new Vector (GetNumElementDofs(elementNumber));
    ComputeBasis(elementNumber, localGridPoint, *basisValues);
    Vector* elementDofs = new Vector (GetElementDofs(elementNumber));

    for (int i=1; i<=GetNumElementDofs(elementNumber); i++)
    {
        Uh += ((*mSolutionVector)((*elementDofs)(i)))*((*basisValues)(i));
    }

    delete basisValues;
    delete elementDofs;

    return Uh;
}


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
}

FE_Solution::~FE_Solution()
{
    delete mMeshReference;
    delete solutionVector;
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

    for (int i=0; i<mMeshReference->GetNumElements(); i++)
    {
        mElementPolySpaceArray[i] = new PolynomialSpace(mPolynomialDegree, mMeshReference->GetElement(i+1)->GetElementType());
    }

    dofStart = new Vector(mMeshReference->GetNumElements() + 1);
    (*dofStart)[0] = mMeshReference->GetNumNodes() + 1;
    for(int i=1; i<dofStart->GetSize(); i++)
    {
        (*dofStart)[i] = (*dofStart)[i-1] + (mPolynomialDegree - 1);
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

Vector FE_Solution::GetSolutionVector()
{
    return *solutionVector;
}

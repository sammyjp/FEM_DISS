#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "FE_Solution.hpp"
#include "Matrix.hpp"
#include "Mesh.hpp"
#include "Vector.hpp"

FE_Solution::FE_Solution(Mesh& mesh, int polynomialDegree)
{
    mMesh = &mesh;
    mPolynomialDegree = polynomialDegree;

    InitialiseElementDofs();
}

FE_Solution::~FE_Solution()
{
    delete mMesh;
    delete dofStart;

    for (int i=0; i<mMesh->GetNumElements(); i++)
    {
        delete mElementPolySpaceArray[i];
    }
    delete[] mElementPolySpaceArray;
}

void FE_Solution::InitialiseElementDofs()
{
    mElementPolySpaceArray = new PolynomialSpace* [mMesh->GetNumElements()];

    for (int i=0; i<mMesh->GetNumElements(); i++)
    {
        mElementPolySpaceArray[i] = new PolynomialSpace(mPolynomialDegree, mMesh->GetElement(i+1)->GetElementType());
    }

    dofStart = new Vector(mMesh->GetNumElements() + 1);
    (*dofStart)[0] = mMesh->GetNumNodes() + 1;
    for(int i=1; i<dofStart->GetSize(); i++)
    {
        (*dofStart)[i] = (*dofStart)[i-1] + (mPolynomialDegree - 1);
    }
}

PolynomialSpace* FE_Solution::GetElementPolynomialSpace(int elementNumber) const
{
    assert(elementNumber >= 1 && elementNumber <= mMesh->GetNumElements());
    return mElementPolySpaceArray[elementNumber - 1];
}


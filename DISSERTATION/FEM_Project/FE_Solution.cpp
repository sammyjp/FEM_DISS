#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "FE_Solution.hpp"
#include "Mesh.hpp"
#include "Vector.hpp"


// Specialised Constructor
FE_Solution::FE_Solution(Mesh& mesh)
{
    mMesh = new Mesh(mesh);
}

// Destructor
FE_Solution::~FE_Solution()
{
    delete mMesh;
}

void FE_Solution::ComputeLinearBasisFunctionValues(int i, Vector& functionValues, Vector& x)
{
    assert(functionValues.GetSize() == x.GetSize());
    assert(i>=0&&i<=mMesh->GetNumElements());

    int N = mMesh->GetNumElements();

    Vector* xGrid = new Vector(mMesh->GetXGridPoints());

    if (i==0)
    {
        for (int j=0; j<x.GetSize(); j++)
        {
            if (((*xGrid)[i] <= x[j]) && (x[j] <= (*xGrid)[i+1]))
            {
                functionValues[j] = 1.0 - x[j]/((*xGrid)[i+1]-(*xGrid)[i]) + i;
            }
            else
            {
                functionValues[j] = 0;
            }
        }
    }
    else if (i==N)
    {
        for (int j=0; j<x.GetSize(); j++)
        {
            if (((*xGrid)[i-1] <= x[j]) && (x[j] <= (*xGrid)[i]))
            {
                functionValues[j] = 1.0 + x[j]/((*xGrid)[i]-(*xGrid)[i-1]) - i;
            }
            else
            {
                functionValues[j] = 0;
            }
        }
    }
    else
    {
        for (int j=0; j<x.GetSize(); j++)
        {
            if (((*xGrid)[i-1] <= x[j]) && (x[j] <= (*xGrid)[i]))
            {
                functionValues[j] = 1.0 + x[j]/((*xGrid)[i]-(*xGrid)[i-1]) - i;
            }
            else if (((*xGrid)[i] < x[j]) && (x[j] <= (*xGrid)[i+1]))
            {
                functionValues[j] = 1.0 - x[j]/((*xGrid)[i+1]-(*xGrid)[i]) + i;
            }
            else
            {
                functionValues[j] = 0;
            }
        }
    }

    delete xGrid;
}

void FE_Solution::ComputeLinearBasisFunctionDerivativeValues(int i, Vector& derivativeValues, Vector& x)
{
    assert(derivativeValues.GetSize() == x.GetSize());
    assert(i>=0&&i<=mMesh->GetNumElements());

    int N = mMesh->GetNumElements();

    Vector* xGrid = new Vector(mMesh->GetXGridPoints());

    if (i==0)
    {
        for (int j=0; j<x.GetSize(); j++)
        {
            if (((*xGrid)[i] < x[j]) && (x[j] < (*xGrid)[i+1]))
            {
                derivativeValues[j] = -1.0/((*xGrid)[i+1]-(*xGrid)[i]);
            }
            else if ((x[j] == (*xGrid)[i]) || (x[j] == (*xGrid)[i+1]))
            {
                derivativeValues[j] = NAN;
            }
            else
            {
                derivativeValues[j] = 0;
            }
        }
    }
    else if (i==N)
    {
        for (int j=0; j<x.GetSize(); j++)
        {
            if (((*xGrid)[i-1] < x[j]) && (x[j] < (*xGrid)[i]))
            {
                derivativeValues[j] = 1.0/((*xGrid)[i]-(*xGrid)[i-1]);
            }
            else if ((x[j] == (*xGrid)[i-1]) || (x[j] == (*xGrid)[i]))
            {
                derivativeValues[j] = NAN;
            }
            else
            {
                derivativeValues[j] = 0;
            }
        }
    }
    else
    {
        for (int j=0; j<x.GetSize(); j++)
        {
            if (((*xGrid)[i-1] < x[j]) && (x[j] < (*xGrid)[i]))
            {
                derivativeValues[j] = 1.0/((*xGrid)[i]-(*xGrid)[i-1]);
            }
            else if (((*xGrid)[i] < x[j]) && (x[j] < (*xGrid)[i+1]))
            {
                derivativeValues[j] = -1.0/((*xGrid)[i+1]-(*xGrid)[i]);
            }
            else if ((x[j] == (*xGrid)[i-1]) || (x[j] == (*xGrid)[i]) || (x[j] == (*xGrid)[i+1]))
            {
                derivativeValues[j] = NAN;
            }
            else
            {
                derivativeValues[j] = 0;
            }
        }
    }
}


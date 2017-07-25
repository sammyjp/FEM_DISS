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
}

FE_Solution::~FE_Solution()
{
    delete mMesh;
}

void FE_Solution::ComputeBasisFunctionValues(int dofNumber, Vector& functionValues, Vector& localQuadraturePoints)
{
    assert(functionValues.GetSize() == localQuadraturePoints.GetSize());


}

void FE_Solution::ComputeLinearBasisFunctionValues(int i, Vector& functionValues, Matrix& x)
{
    assert(functionValues.GetSize() == x.GetNumberOfColumns());
    i--;
    assert(i>=0&&i<=mMesh->GetNumElements());

    int N = mMesh->GetNumElements();

    Vector* xGrid = new Vector(mMesh->GetAllGridPoints().GetRowAsVector(1));

    if (i==0)
    {
        for (int j=0; j<x.GetNumberOfColumns(); j++)
        {
            if (((*xGrid)[i] <= x(1,j+1)) && (x(1,j+1) <= (*xGrid)[i+1]))
            {
                functionValues[j] = 1.0 - x(1,j+1)/((*xGrid)[i+1]-(*xGrid)[i]) + i;
            }
            else
            {
                functionValues[j] = 0;
            }
        }
    }
    else if (i==N)
    {
        for (int j=0; j<x.GetNumberOfColumns(); j++)
        {
            if (((*xGrid)[i-1] <= x(1,j+1)) && (x(1,j+1) <= (*xGrid)[i]))
            {
                functionValues[j] = 1.0 + x(1,j+1)/((*xGrid)[i]-(*xGrid)[i-1]) - i;
            }
            else
            {
                functionValues[j] = 0;
            }
        }
    }
    else
    {
        for (int j=0; j<x.GetNumberOfColumns(); j++)
        {
            if (((*xGrid)[i-1] <= x(1,j+1)) && (x(1,j+1) <= (*xGrid)[i]))
            {
                functionValues[j] = 1.0 + x(1,j+1)/((*xGrid)[i]-(*xGrid)[i-1]) - i;
            }
            else if (((*xGrid)[i] < x(1,j+1)) && (x(1,j+1) <= (*xGrid)[i+1]))
            {
                functionValues[j] = 1.0 - x(1,j+1)/((*xGrid)[i+1]-(*xGrid)[i]) + i;
            }
            else
            {
                functionValues[j] = 0;
            }
        }
    }

    delete xGrid;
}

void FE_Solution::ComputeLinearBasisFunctionDerivativeValues(int i, Vector& derivativeValues, Matrix& x)
{
    assert(derivativeValues.GetSize() == x.GetNumberOfColumns());
    i--;
    assert(i>=0&&i<=mMesh->GetNumElements());

    int N = mMesh->GetNumElements();

    Vector* xGrid = new Vector(mMesh->GetAllGridPoints().GetRowAsVector(1));

    if (i==0)
    {
        for (int j=0; j<x.GetNumberOfColumns(); j++)
        {
            if (((*xGrid)[i] < x(1,j+1)) && (x(1,j+1) < (*xGrid)[i+1]))
            {
                derivativeValues[j] = -1.0/((*xGrid)[i+1]-(*xGrid)[i]);
            }
            else if ((x(1,j+1) == (*xGrid)[i]) || (x(1,j+1) == (*xGrid)[i+1]))
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
        for (int j=0; j<x.GetNumberOfColumns(); j++)
        {
            if (((*xGrid)[i-1] < x(1,j+1)) && (x(1,j+1) < (*xGrid)[i]))
            {
                derivativeValues[j] = 1.0/((*xGrid)[i]-(*xGrid)[i-1]);
            }
            else if ((x(1,j+1) == (*xGrid)[i-1]) || (x(1,j+1) == (*xGrid)[i]))
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
        for (int j=0; j<x.GetNumberOfColumns(); j++)
        {
            if (((*xGrid)[i-1] < x(1,j+1)) && (x(1,j+1) < (*xGrid)[i]))
            {
                derivativeValues[j] = 1.0/((*xGrid)[i]-(*xGrid)[i-1]);
            }
            else if (((*xGrid)[i] < x(1,j+1)) && (x(1,j+1) < (*xGrid)[i+1]))
            {
                derivativeValues[j] = -1.0/((*xGrid)[i+1]-(*xGrid)[i]);
            }
            else if ((x(1,j+1) == (*xGrid)[i-1]) || (x(1,j+1) == (*xGrid)[i]) || (x(1,j+1) == (*xGrid)[i+1]))
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


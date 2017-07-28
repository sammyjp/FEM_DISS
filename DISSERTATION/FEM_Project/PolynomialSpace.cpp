#include <cmath>
#include <cassert>

#include "PolynomialSpace.hpp"

PolynomialSpace::PolynomialSpace(int polynomialDegree, int elementType)
{
    mPolynomialDegree = polynomialDegree;
    mElementType = elementType;
    mNumElementDofs = polynomialDegree + 1;
}

double PolynomialSpace::EvaluateNthLegendrePolynomial(int n, double pointToEvaluate)
{
    assert (n >= 0);

    if (n==0)
    {
        return 1;
    }
    else if (n==1)
    {
        return pointToEvaluate;
    }
    else
    {
        double* pn_prev = new double;
        (*pn_prev) = 1;
        double* pn = new double;
        (*pn) = pointToEvaluate;
        double legendrePoint;

        for (int i=1; i<n; i++)
        {
            legendrePoint = ((2.0*i)+1.0)/(i+1.0)*pointToEvaluate*(*pn) - (i/(i+1.0))*(*pn_prev);
            (*pn_prev) = (*pn);
            (*pn) = legendrePoint;
        }

        delete pn_prev;
        delete pn;

        return legendrePoint;
    }
}

double PolynomialSpace::EvaluateNthLobattoPolynomial(int n, double pointToEvaluate)
{
    assert(n >= 0);

    if (n==0)
    {
        return (1.0 - pointToEvaluate)/2.0;
    }
    else if (n==1)
    {
        return (1.0 + pointToEvaluate)/2.0;
    }
    else
    {
        double I=0;
        int n_q = ceil((n+1.0)/2.0);
        Matrix* gaussPoints = new Matrix(1,n_q);
        Vector* gaussWeights = new Vector(n_q);
        Vector* functionPoints = new Vector(n_q);

        QuadratureLibrary().GetQuadrature(0, n_q, *gaussWeights, *gaussPoints);

        Vector* gaussPointsVec = new Vector(gaussPoints->GetRowAsVector(1));
        delete gaussPoints;

        for (int i=0; i<gaussPointsVec->GetSize(); i++)
        {
            (*gaussPointsVec)[i] = -1.0 + (((pointToEvaluate + 1.0)/2.0)*(1+(*gaussPointsVec)[i]));
        }

        EvaluateNthLegendrePolynomial(n-1, *gaussPointsVec, *functionPoints);
        delete gaussPointsVec;

        for (int i=1; i<=n_q; i++)
        {
            I += ((pointToEvaluate+1.0)/2.0) * (*gaussWeights)(i) * (*functionPoints)(i);
        }

        delete gaussWeights;
        delete functionPoints;

        I = sqrt((2.0*n - 1.0)/2.0)*I;

        return I;
    }
}
void PolynomialSpace::EvaluateNthLegendrePolynomial(int n, Vector& pointsToEvaluate, Vector& legendrePoints)
{
    assert(pointsToEvaluate.GetSize() == legendrePoints.GetSize());
    assert(n >= 0);

    Vector* pn_prev = new Vector(pointsToEvaluate.GetSize());
    (*pn_prev) = 1;
    Vector* pn = new Vector(pointsToEvaluate);
    Vector* pn_next = new Vector(pointsToEvaluate.GetSize());

    if(n==0)
    {
        (*pn) = (*pn_prev);
    }

    if (n>1)
    {
        for (int i=1; i<n; i++)
        {
            for (int j=0; j<pn->GetSize(); j++)
            {
                (*pn_next)[j] = (pointsToEvaluate[j])*((*pn)[j])*((2.0*i)+1.0)/(i+1.0) - ((*pn_prev)[j])*(i/(i+1.0));
            }
            (*pn_prev) = (*pn);
            (*pn) = (*pn_next);
        }
    }

    legendrePoints = (*pn);

    delete pn_prev;
    delete pn;
    delete pn_next;
}

void PolynomialSpace::EvaluateNthLobattoPolynomial(int n, Vector& pointsToEvaluate, Vector& lobattoPoints)
{
    assert(pointsToEvaluate.GetSize() == lobattoPoints.GetSize());
    assert(n >= 0);

    if (n==0)
    {
        for (int i=0; i<pointsToEvaluate.GetSize(); i++)
        {
            lobattoPoints[i] = (1.0 - pointsToEvaluate[i])/2.0;
        }
    }
    else if (n==1)
    {
        for (int i=0; i<pointsToEvaluate.GetSize(); i++)
        {
            lobattoPoints[i] = (1.0 + pointsToEvaluate[i])/2.0;
        }
    }
    else
    {
        Vector* legendrePoints_prev = new Vector(pointsToEvaluate.GetSize());
        Vector* legendrePoints = new Vector(pointsToEvaluate.GetSize());


        EvaluateNthLegendrePolynomial(mNumElementDofs - 1, pointsToEvaluate, *legendrePoints_prev);
        EvaluateNthLegendrePolynomial(mNumElementDofs, pointsToEvaluate, *legendrePoints);

        for (int i=0; i<pointsToEvaluate.GetSize(); i++)
        {
            lobattoPoints[i] = (mNumElementDofs-1) * ((*legendrePoints_prev)[i] - (pointsToEvaluate[i]*(*legendrePoints)[i]));
        }

        delete legendrePoints_prev;
        delete legendrePoints;
    }

}

int PolynomialSpace::GetNumElementDofs()
{
    return mNumElementDofs;
}

void PolynomialSpace::ComputeBasis(double localGridPoint, Vector& basisValues)
{
    assert(basisValues.GetSize() == mNumElementDofs);

    switch(mElementType)
    {
    case 0:
        {
            for (int i=0; i<mNumElementDofs; i++)
            {
                basisValues[i] = EvaluateNthLobattoPolynomial(i, localGridPoint);
            }
        } break;
    }
}

void PolynomialSpace::ComputeGradBasis(double localGridPoint, Matrix& gradBasisValues)
{
    assert(gradBasisValues.GetNumberOfColumns() == mNumElementDofs);

    switch(mElementType)
    {
    case 0:
        {
            assert(gradBasisValues.GetNumberOfRows() == 1);

            gradBasisValues(1,1) = -0.5;
            gradBasisValues(1,2) = 0.5;

            for (int i=3; i<=mNumElementDofs; i++)
            {
                gradBasisValues(1,i) =  sqrt((2.0*(i-1) - 1.0)/2.0)*EvaluateNthLegendrePolynomial(i-2, localGridPoint);
            }
        } break;
    }
}

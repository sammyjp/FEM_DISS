#include <cmath>
#include <cassert>

#include "Quadrature.hpp"
#include "Vector.hpp"

Quadrature::Quadrature(double startPoint, double endPoint, int numNodes)
{
    mStartPoint = startPoint;
    mEndPoint = endPoint;
    mNumNodes = numNodes;

    mGridPoints = new Vector(mNumNodes);
}

Quadrature::Quadrature(Quadrature& otherQuadrature)
{
    mStartPoint = otherQuadrature.mStartPoint;
    mEndPoint = otherQuadrature.mEndPoint;
    mNumNodes = otherQuadrature.mNumNodes;
}

Quadrature::~Quadrature()
{
    delete mGridPoints;
}

Vector Quadrature::GetUniformGridPoints()
{
    mUniformStepSize = (mEndPoint-mStartPoint)/(mNumNodes-1);

    for (int i=0; i<mGridPoints->GetSize(); i++)
    {
        (*mGridPoints)[i] = mStartPoint + (i*mUniformStepSize);
    }

    return (*mGridPoints);
}

//////////////////////////////////////////////////
/*               GAUSS QUADRATURE               */
//////////////////////////////////////////////////

void Quadrature::TransformGQPoints()
{
    for (int i=0; i<mGridPoints->GetSize(); i++)
    {
        (*mGridPoints)[i] = mStartPoint + 0.5*(mEndPoint-mStartPoint)*((*mGridPoints)[i]+1);
    }
}

void Quadrature::ComputeGQWeights(Vector& GQWeights)
{

}

void Quadrature::EvaluateNthLegendrePolynomial(int n, Vector& pointsToEvaluate, Vector& legendrePoints)
{
    assert(n>=1);
    Vector* pn_prev = new Vector(pointsToEvaluate.GetSize());
    (*pn_prev) = 1;
    Vector* pn = new Vector(pointsToEvaluate);
    Vector* pn_next = new Vector(pointsToEvaluate.GetSize());

    if (i>=2)
    {
        for (int i=2; i<=n; i++)
        {
            for (int j=0; j<pn_next->GetSize(); j++)
            {
                (*pn_next)[j] = (pointsToEvaluate[j])*((*pn)[j])*((2.0*n)+1.0)/(n+1.0) - ((*pn_prev)[j])*(n/(n+1.0));
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

void Quadrature::NewtonForLegendre(double tolerance, Vector& initialGuess)
{

}

Vector Quadrature::GetGQPoints()
{


    TransformGQPoints();

    return (*mGridPoints);
}

double Quadrature::GaussQuadrature(Vector& functionPoints)
{
    assert(functionPoints.GetSize() == mGridPoints->GetSize());
    double I = 0;
    Vector* GQWeights = new Vector(mGridPoints->GetSize());

    ComputeGQWeights((*GQWeights));

    for (int i=0; i<functionPoints.GetSize(); i++)
    {
        I += (*GQWeights)[i]*functionPoints[i];
    }

    delete GQWeights;

    return I;
}

//////////////////////////////////////////////////
/*                 Trapezoid Rule               */
//////////////////////////////////////////////////

double Quadrature::TrapezoidRule(Vector& functionPoints)
{
    assert(functionPoints.GetSize() == mGridPoints->GetSize());
    double I = 0;

    for (int i=0; i<functionPoints.GetSize()-1; i++)
    {
        I += ((*mGridPoints)[i+1] - (*mGridPoints)[i]) * (functionPoints[i+1] + functionPoints[i]);
    }

    return I;
}

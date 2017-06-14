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

    return I;
}

//////////////////////////////////////////////////
/*                 Trapezoid Rule               */
//////////////////////////////////////////////////

double Quadrature::TrapezoidRule(Vector& functionPoints)
{
    assert(functionPoints.GetSize() == mGridPoints->GetSize());
    double I = 0;

    for (int i=0; i<functionPoints.GetSize(); i++)
    {
        if (i==0 || i==(functionPoints.GetSize()-1))
        {
            I += (mUniformStepSize/2.0)*functionPoints[i];
        }
        else
        {
            I += mUniformStepSize*functionPoints[i];
        }
    }

    return I;
}

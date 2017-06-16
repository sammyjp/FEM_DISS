#include <cmath>
#include <cassert>

#include "Quadrature.hpp"
#include "Vector.hpp"

const double Pi = 4*atan(1);

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
    assert(GQWeights.GetSize() == mGridPoints->GetSize());

    Vector* localGridPoints = new Vector(mGridPoints->GetSize());

    for (int i=0; i<localGridPoints->GetSize(); i++)
    {
        (*localGridPoints)[i] = (2.0*((*mGridPoints)[i]-mStartPoint))/(mEndPoint-mStartPoint) -1.0;
    }

    Vector* legendreDeriv = new Vector(GQWeights.GetSize());

    EvaluateNthLegendrePolynomialFirstDerivative((*localGridPoints), (*legendreDeriv));

    for (int i=0; i<GQWeights.GetSize(); i++)
    {
        GQWeights[i] = 2.0/((1.0-pow((*localGridPoints)[i],2.0))*(pow((*legendreDeriv)[i],2.0)));
    }

    delete localGridPoints;
    delete legendreDeriv;
}

void Quadrature::EvaluateNthLegendrePolynomial(Vector& pointsToEvaluate, Vector& legendrePoints)
{
    Vector* pn_prev = new Vector(pointsToEvaluate.GetSize());
    (*pn_prev) = 1;
    Vector* pn = new Vector(pointsToEvaluate);
    Vector* pn_next = new Vector(pointsToEvaluate.GetSize());

    if (mNumNodes>1)
    {
        for (int i=1; i<mNumNodes; i++)
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

void Quadrature::EvaluateNthLegendrePolynomialFirstDerivative(Vector& pointsToEvaluate, Vector& legendrePoints)
{
    Vector* pn_prev = new Vector(pointsToEvaluate.GetSize());
    (*pn_prev) = 0;
    Vector* pn = new Vector(pointsToEvaluate.GetSize());
    (*pn) = 1;
    Vector* pn_next = new Vector(pointsToEvaluate.GetSize());

    if (mNumNodes>1)
    {
        for (int i=1; i<mNumNodes; i++)
        {
            for (int j=0; j<pn_next->GetSize(); j++)
            {
                (*pn_next)[j] = (pointsToEvaluate[j])*((*pn)[j])*(((2.0*i)+1.0)/i) - ((*pn_prev)[j])*((i+1.0)/i);
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

void Quadrature::NewtonForLegendre(double tolerance, int maxIterations, Vector& initialGuess, Vector& legendreRoots)
{
    assert(initialGuess.GetSize() == legendreRoots.GetSize());
    assert(initialGuess.GetSize() == mNumNodes);

    int iterationCounter = 0;
    double error;
    Vector* error_vec = new Vector (legendreRoots.GetSize());
    Vector* legendreEval = new Vector (legendreRoots.GetSize());
    Vector* legendreDeriv = new Vector (legendreRoots.GetSize());

    do
    {
        EvaluateNthLegendrePolynomial(initialGuess, *legendreEval);
        EvaluateNthLegendrePolynomialFirstDerivative(initialGuess, *legendreDeriv);

        for (int i=0; i<legendreRoots.GetSize(); i++)
        {
            legendreRoots[i] = initialGuess[i] - (*legendreEval)[i]/(*legendreDeriv)[i];
        }

        (*error_vec) = legendreRoots - initialGuess;
        initialGuess = legendreRoots;

        error = error_vec->CalculateNorm(2);
        iterationCounter++;
    } while (error >= tolerance && iterationCounter <= maxIterations);

    delete error_vec;
    delete legendreEval;
    delete legendreDeriv;
}

Vector Quadrature::GetGQPoints()
{
    double tolerance = 1e-12;
    int maxIterations = 1000;
    Vector* initialGuess = new Vector (mNumNodes);

    for (int i=1; i<=mNumNodes; i++)
    {
        (*initialGuess)[i-1] = -cos((((2.0*i)-1.0)/(2.0*mNumNodes))*Pi);
    }

    NewtonForLegendre(tolerance, maxIterations, (*initialGuess), (*mGridPoints));
    delete initialGuess;

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
        I += (*GQWeights)[i]*(0.5*(mEndPoint-mStartPoint)*functionPoints[i]);
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

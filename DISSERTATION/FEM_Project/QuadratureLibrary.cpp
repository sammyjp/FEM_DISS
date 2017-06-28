#include <cmath>
#include <cassert>

#include "QuadratureLibrary.hpp"
#include "Vector.hpp"

const double Pi = 4*atan(1);

QuadratureLibrary::QuadratureLibrary()
{
}

//////////////////////////////////////////////////
/*               GAUSS QUADRATURE               */
//////////////////////////////////////////////////

void QuadratureLibrary::TransformGQPoints(double xStart, double xEnd, Vector& gaussPoints)
{
    for (int i=0; i<gaussPoints.GetSize(); i++)
    {
        gaussPoints[i] = xStart + 0.5*(xEnd-xStart)*(gaussPoints[i] + 1);
    }
}

void QuadratureLibrary::ComputeGQWeights(Vector& GQWeights, Vector& gaussPoints)
{
    assert(GQWeights.GetSize() == gaussPoints.GetSize());

    Vector* legendreDeriv = new Vector(GQWeights.GetSize());

    EvaluateNthLegendrePolynomialFirstDerivative(gaussPoints, (*legendreDeriv));

    for (int i=0; i<GQWeights.GetSize(); i++)
    {
        GQWeights[i] = 2.0/((1.0-pow(gaussPoints[i],2.0))*(pow((*legendreDeriv)[i],2.0)));
    }

    delete legendreDeriv;
}

void QuadratureLibrary::EvaluateNthLegendrePolynomial(Vector& pointsToEvaluate, Vector& legendrePoints)
{
    Vector* pn_prev = new Vector(pointsToEvaluate.GetSize());
    (*pn_prev) = 1;
    Vector* pn = new Vector(pointsToEvaluate);
    Vector* pn_next = new Vector(pointsToEvaluate.GetSize());

    if (pointsToEvaluate.GetSize()>1)
    {
        for (int i=1; i<pointsToEvaluate.GetSize(); i++)
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

void QuadratureLibrary::EvaluateNthLegendrePolynomialFirstDerivative(Vector& pointsToEvaluate, Vector& legendrePoints)
{
    Vector* pn_prev = new Vector(pointsToEvaluate.GetSize());
    (*pn_prev) = 0;
    Vector* pn = new Vector(pointsToEvaluate.GetSize());
    (*pn) = 1;
    Vector* pn_next = new Vector(pointsToEvaluate.GetSize());

    if (pointsToEvaluate.GetSize()>1)
    {
        for (int i=1; i<pointsToEvaluate.GetSize(); i++)
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

void QuadratureLibrary::NewtonForLegendre(double tolerance, int maxIterations, Vector& initialGuess, Vector& legendreRoots)
{
    assert(initialGuess.GetSize() == legendreRoots.GetSize());

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

void QuadratureLibrary::GetGaussQuadraturePoints(Vector& gaussPoints)
{
    double tolerance = 1e-12;
    int maxIterations = 1000;
    Vector* initialGuess = new Vector (gaussPoints.GetSize());

    for (int i=1; i<=gaussPoints.GetSize(); i++)
    {
        (*initialGuess)[i-1] = -cos((((2.0*i)-1.0)/(2.0*gaussPoints.GetSize()))*Pi);
    }

    NewtonForLegendre(tolerance, maxIterations, (*initialGuess), gaussPoints);
    delete initialGuess;
}

double QuadratureLibrary::GaussQuadrature(int elementType, Vector& gaussPoints, Vector& functionPoints, double mappingJacobianDeterminant)
{
    Vector* GQWeights = new Vector(gaussPoints.GetSize());
    ComputeGQWeights((*GQWeights), gaussPoints);

    double I = 0;

    switch(elementType)
    {
        case 0:
        {
            assert(functionPoints.GetSize() == gaussPoints.GetSize());

            for (int i=0; i<functionPoints.GetSize(); i++)
            {
                I += (*GQWeights)[i]*mappingJacobianDeterminant*functionPoints[i];
            }
        } break;
    }

    delete GQWeights;
    return I;
}

//////////////////////////////////////////////////
/*                 Trapezoid Rule               */
//////////////////////////////////////////////////

double QuadratureLibrary::TrapezoidRule1D(Vector& gridPoints, Vector& functionPoints)
{
    assert(gridPoints.GetSize() == functionPoints.GetSize());

    double I = 0;

    for (int i=0; i<functionPoints.GetSize()-1; i++)
    {
        I += (gridPoints[i+1] - gridPoints[i]) * (functionPoints[i+1] + functionPoints[i]);
    }

    return I;
}

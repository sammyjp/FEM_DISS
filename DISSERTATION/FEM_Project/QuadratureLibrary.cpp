#include <cmath>
#include <cassert>

#include "QuadratureLibrary.hpp"

const double Pi = 4*atan(1);

QuadratureLibrary::QuadratureLibrary()
{
}

//////////////////////////////////////////////////
/*               GAUSS QUADRATURE               */
//////////////////////////////////////////////////

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

void QuadratureLibrary::ComputeGQPoints(Vector& gaussPoints)
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

void QuadratureLibrary::MapReferenceQuadrilateralToTriangle(Matrix& quadrilateralPoints, Matrix& trianglePoints)
{
    assert(trianglePoints.GetNumberOfRows() == quadrilateralPoints.GetNumberOfRows());
    assert(trianglePoints.GetNumberOfColumns() == quadrilateralPoints.GetNumberOfColumns());
    assert(trianglePoints.GetNumberOfRows() == 2);

    for (int j=1; j<=trianglePoints.GetNumberOfColumns(); j++)
    {
        trianglePoints(1,j) = 0.5*(-1 + quadrilateralPoints(1,j) - quadrilateralPoints(2,j) - quadrilateralPoints(1,j)*quadrilateralPoints(2,j));
        trianglePoints(2,j) = quadrilateralPoints(2,j);
    }
}


void QuadratureLibrary::GetQuadrature(const int elementType, const int n_q, Vector& weights, Matrix& gaussPoints)
{
    int n;
    assert(gaussPoints.GetNumberOfColumns() == n_q);
    assert(weights.GetSize() == n_q);

    switch(elementType)
    {
    case 0:
        {
            assert(gaussPoints.GetNumberOfRows() == 1);

            n = n_q;

            Vector* gaussVec = new Vector(n);
            Vector* weightVec = new Vector(n);
            ComputeGQPoints(*gaussVec);
            ComputeGQWeights(*weightVec, *gaussVec);

            for (int j=1; j<=n; j++)
            {
                gaussPoints(1,j) = (*gaussVec)(j);
            }

            weights = (*weightVec);

            delete gaussVec;
            delete weightVec;
        } break;

    case 1:
        {
            assert(gaussPoints.GetNumberOfRows() == 2);

            n = sqrt(n_q);

            Vector* gaussVec = new Vector(n);
            Vector* weightVec = new Vector(n);
            ComputeGQPoints(*gaussVec);
            ComputeGQWeights(*weightVec, *gaussVec);

            for (int i=1; i<=n; i++)
            {
                for (int j=1; j<=n; j++)
                {
                    gaussPoints(1,((i-1)*n)+j) = (*gaussVec)(i);
                    gaussPoints(2,((i-1)*n)+j) = (*gaussVec)(j);
                    weights((i-1)*n+j) = (*weightVec)(i)*(*weightVec)(j);
                }
            }

            delete gaussVec;
            delete weightVec;

            Matrix* quadrilateralPoints = new Matrix(gaussPoints);
            MapReferenceQuadrilateralToTriangle(*quadrilateralPoints, gaussPoints);
            delete quadrilateralPoints;

            for (int j=1; j<=weights.GetSize(); j++)
            {
                weights(j) = weights(j)*((1-gaussPoints(2,j))*0.5);
            }
        } break;

    case 2:
        {
            assert(gaussPoints.GetNumberOfRows() == 2);

            n = sqrt(n_q);

            Vector* gaussVec = new Vector(n);
            Vector* weightVec = new Vector(n);
            ComputeGQPoints(*gaussVec);
            ComputeGQWeights(*weightVec, *gaussVec);

            for (int i=1; i<=n; i++)
            {
                for (int j=1; j<=n; j++)
                {
                    gaussPoints(1,((i-1)*n)+j) = (*gaussVec)(i);
                    gaussPoints(2,((i-1)*n)+j) = (*gaussVec)(j);
                    weights((i-1)*n+j) = (*weightVec)(i)*(*weightVec)(j);
                }
            }

            delete gaussVec;
            delete weightVec;
        } break;
    }
}


#ifndef QUADRATURELIBRARYHEADERDEF
#define QUADRATURELIBRARYHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class QuadratureLibrary
{
private:

    void ComputeGQWeights(Vector& GQWeights, Vector& gaussPoints);
    void ComputeGQPoints(Vector& gaussPoints);

    void NewtonForLegendre(double tolerance, int maxIterations, Vector& initialGuess, Vector& legendreRoots);

    void EvaluateNthLegendrePolynomial(Vector& pointsToEvaluate, Vector& legendrePoints);
    void EvaluateNthLegendrePolynomialFirstDerivative(Vector& pointsToEvaluate, Vector& legendrePoints);

    void MapReferenceQuadrilateralToTriangle(Matrix& quadrilateralPoints, Matrix& trianglePoints);

public:
    // Default Constructor
    QuadratureLibrary();

    void GetQuadrature(const int elementType, const int n_q, Vector& weights, Matrix& gaussPoints);
};

#endif

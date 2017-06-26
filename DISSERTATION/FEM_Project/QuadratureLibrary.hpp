#ifndef QUADRATURELIBRARYHEADERDEF
#define QUADRATURELIBRARYHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class QuadratureLibrary
{
private:

    void TransformGQPoints();
    void ComputeGQWeights(Vector& GQWeights, Vector& gaussPoints);

    void NewtonForLegendre(double tolerance, int maxIterations, Vector& initialGuess, Vector& legendreRoots);

    void EvaluateNthLegendrePolynomial(Vector& pointsToEvaluate, Vector& legendrePoints);
    void EvaluateNthLegendrePolynomialFirstDerivative(Vector& pointsToEvaluate, Vector& legendrePoints);

public:
    // Default Constructor
    QuadratureLibrary();

    void GetGaussQuadraturePoints(Vector& gaussPoints);
    double GaussQuadrature(int elementType, Vector& gaussPoints, Vector& functionPoints, double mappingJacobianDeterminant);

    double TrapezoidRule1D(Vector& gridPoints, Vector& functionPoints);

};

#endif

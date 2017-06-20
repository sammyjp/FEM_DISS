#ifndef QUADRATUREHEADERDEF
#define QUADRATUREHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class Quadrature
{
private:

    double mStartPoint;
    double mEndPoint;
    int mNumPoints;
    int mElementType;
    Matrix* mVertices;

    Vector* mGridPoints;

    double mUniformStepSize;

    void TransformGQPoints();
    void ComputeGQWeights(Vector& GQWeights);

    void NewtonForLegendre(double tolerance, int maxIterations, Vector& initialGuess, Vector& legendreRoots);

    void EvaluateNthLegendrePolynomial(Vector& pointsToEvaluate, Vector& legendrePoints);
    void EvaluateNthLegendrePolynomialFirstDerivative(Vector& pointsToEvaluate, Vector& legendrePoints);

public:
    // Specialised Constructor
    Quadrature(int numPoints, Matrix& ElementVertices, int ElementType);
    Quadrature(double startPoint, double endPoint, int numNodes);

    // Copy Constructor
    Quadrature(Quadrature& otherQuadrature);

    // Destructor
    ~Quadrature();

    Vector GetUniformGridPoints();

    Vector GetGQPoints();
    double GaussQuadrature(Vector& functionPoints);

    double TrapezoidRule(Vector& functionPoints);


};

#endif

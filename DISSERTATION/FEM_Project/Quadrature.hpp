#ifndef QUADRATUREHEADERDEF
#define QUADRATUREHEADERDEF

#include "Vector.hpp"

class Quadrature
{
private:

    double mStartPoint;
    double mEndPoint;
    int mNumNodes;

    Vector* mGridPoints;

    double mUniformStepSize;

    void TransformGQPoints();
    void ComputeGQWeights(Vector& GQWeights);
    void EvaluateNthLegendrePolynomial(int n, Vector& pointsToEvaluate, Vector& legendrePoints);
    void NewtonForLegendre();

public:

    // Specialised Constructor
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

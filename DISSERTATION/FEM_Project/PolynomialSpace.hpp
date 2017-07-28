#ifndef POLYNOMIALSPACEHEADERDEF
#define POLYNOMIALSPACEHEADERDEF

#include "Element.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

class PolynomialSpace
{
private:
    int mPolynomialDegree;
    int mElementType;
    int mNumElementDofs;

    double EvaluateNthLegendrePolynomial(int n, double pointToEvaluate);
    double EvaluateNthLobattoPolynomial(int n, double pointToEvaluate);
    void EvaluateNthLegendrePolynomial(int n, Vector& pointsToEvaluate, Vector& legendrePoints);
    void EvaluateNthLobattoPolynomial(int n, Vector& pointsToEvaluate, Vector& lobattoPoints);

public:

    PolynomialSpace(int polynomialDegree, int elementType);

    int GetNumElementDofs();

    void ComputeBasis(double localGridPoint, Vector& basisValues);

    void ComputeGradBasis(double localGridPoint, Matrix& gradBasisValues);


};

#endif

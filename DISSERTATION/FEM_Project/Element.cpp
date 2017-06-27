#include <cmath>
#include <cassert>

#include "Element.hpp"

Element::~Element()
{
    delete mElementConnectivityArray;
}

Vector Element::GetElementConnectivityArray() const
{
    return *mElementConnectivityArray;
}

void Element::ComputeElementQuadraturePoints(Vector& quadraturePoints)
{
    QuadratureLibrary().GetGaussQuadraturePoints(quadraturePoints);
}

double Element::PerformElementQuadrature(Vector& quadraturePoints, Vector& functionPoints, Matrix& mappingJacobian)
{
    return QuadratureLibrary().GaussQuadrature(GetElementType(), quadraturePoints, functionPoints, mappingJacobian.CalculateDeterminant());
}

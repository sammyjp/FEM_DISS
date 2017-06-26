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

void Element::PerformElementQuadrature(Vector& quadraturePoints, Vector& functionPoints)
{

}

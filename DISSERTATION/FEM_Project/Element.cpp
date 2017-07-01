#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "Mesh.hpp"

Element::~Element()
{
    delete mElementConnectivityArray;
}

Matrix Element::GetElementCoordinates() const
{
    return mMeshReference->GetGridPoints((*mElementConnectivityArray));
}

Vector Element::GetElementConnectivityArray() const
{
    return *mElementConnectivityArray;
}

void Element::GetQuadrature(const int n_q, Vector& quadratureWeights, Matrix& quadraturePoints)
{
    QuadratureLibrary().Quadrature(GetElementType(), n_q, quadratureWeights, quadraturePoints);

    Matrix* localPoints = new Matrix (quadraturePoints);
    MapLocalToGlobal(*localPoints, quadraturePoints);
    delete localPoints;
}

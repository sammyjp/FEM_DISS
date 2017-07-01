#include <cmath>
#include <cassert>

#include "Quadrilateral.hpp"

Quadrilateral::Quadrilateral(Vector& elementConnectivity, Mesh& meshReference)
{
    mElementConnectivityArray = new Vector(elementConnectivity);
    mMeshReference = &meshReference;
}

Quadrilateral::Quadrilateral(Quadrilateral& otherQuadrilateral)
{
    mElementConnectivityArray = new Vector (*otherQuadrilateral.mElementConnectivityArray);
}

int Quadrilateral::GetElementType() const
{
    return ElementType::Quadrilateral;
}

void Quadrilateral::MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords)
{

}

void Quadrilateral::MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords)
{

}

void Quadrilateral::ComputeMappingJacobian(Vector& pointToEval, Matrix& Jacobian)
{

}

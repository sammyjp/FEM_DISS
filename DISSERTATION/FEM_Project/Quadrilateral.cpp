#include <cmath>
#include <cassert>

#include "Quadrilateral.hpp"

Quadrilateral::Quadrilateral(Vector& elementConnectivity)
{
    mElementConnectivityArray = new Vector(elementConnectivity);
}

Quadrilateral::Quadrilateral(Quadrilateral& otherQuadrilateral)
{
    mElementConnectivityArray = new Vector (*otherQuadrilateral.mElementConnectivityArray);
}

int Quadrilateral::GetElementType() const
{
    return ElementType::Quadrilateral;
}

void Quadrilateral::MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords)
{

}

void Quadrilateral::MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords)
{

}

void Quadrilateral::ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian)
{

}



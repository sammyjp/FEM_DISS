#include <cmath>
#include <cassert>

#include "Quadrilateral.hpp"

Quadrilateral::Quadrilateral(Vector& elementConnectivity)
{
    mElementConnectivity = new Vector(elementConnectivity);
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



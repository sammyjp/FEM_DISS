#include <cmath>
#include <cassert>

#include "Interval.hpp"

Interval::Interval(Vector& elementConnectivity)
{
    mElementConnectivityArray = new Vector(elementConnectivity);
}

Interval::Interval(Interval& otherInterval)
{
    mElementConnectivityArray = new Vector (*otherInterval.mElementConnectivityArray);
}

int Interval::GetElementType() const
{
    return ElementType::Interval;
}

void Interval::MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords)
{
    assert(localCoords.GetNumberOfRows() == globalCoords.GetNumberOfRows());
    assert(localCoords.GetNumberOfColumns() == globalCoords.GetNumberOfColumns());
    assert(localCoords.GetNumberOfRows() == 1);
    assert(nodes.GetNumberOfRows() == 1);
    assert(nodes.GetNumberOfColumns() == 2);

    for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
    {
        globalCoords(1,j) = nodes(1,1) + 0.5*(nodes(1,2)-nodes(1,1))*(localCoords(1,j)+1);
    }
}

void Interval::MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords)
{
    assert(globalCoords.GetNumberOfRows() == localCoords.GetNumberOfRows());
    assert(globalCoords.GetNumberOfColumns() == localCoords.GetNumberOfColumns());
    assert(globalCoords.GetNumberOfRows() == 1);
    assert(nodes.GetNumberOfRows() == 1);
    assert(nodes.GetNumberOfColumns() == 2);

    for (int j=1; j<=localCoords.GetNumberOfColumns(); j++)
    {
        localCoords(1,j) = (nodes(1,1) + nodes(1,2) - 2*globalCoords(1,j))/(nodes(1,1) - nodes(1,2));
    }
}

void Interval::ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian)
{
    assert(Jacobian.GetNumberOfRows() == Jacobian.GetNumberOfColumns());
    assert(Jacobian.GetNumberOfRows() == 1);
    assert(nodes.GetNumberOfRows() == 1);
    assert(nodes.GetNumberOfColumns() == 2);

    Jacobian(1,1) = 0.5*(nodes(1,2) - nodes(1,1));
}

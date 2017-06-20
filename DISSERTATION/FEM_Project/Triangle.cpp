#include <cmath>
#include <cassert>

#include "Triangle.hpp"

Triangle::Triangle(Vector& elementConnectivity)
{
    mElementConnectivityArray = new Vector(elementConnectivity);
}

Triangle::Triangle(Triangle& otherTriangle)
{
    mElementConnectivityArray = new Vector (*otherTriangle.mElementConnectivityArray);
}

int Triangle::GetElementType() const
{
    return ElementType::Triangle;
}

void Triangle::MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords)
{
    assert(localCoords.GetNumberOfRows() == globalCoords.GetNumberOfRows());
    assert(localCoords.GetNumberOfColumns() == globalCoords.GetNumberOfColumns());
    assert(localCoords.GetNumberOfRows() == 2);
    assert(nodes.GetNumberOfRows() == 2);
    assert(nodes.GetNumberOfColumns() == 3);

    for (int i=1; i<=globalCoords.GetNumberOfRows(); i++)
    {
        for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
        {
            globalCoords(i,j) = -((localCoords(1,j) + localCoords(2,j))/2.0)*nodes(i,1) + ((localCoords(1,j) + 1)/2)*nodes(i,2) + ((localCoords(2,j) + 1)/2.0)*nodes(i,3);
        }
    }
}

void Triangle::MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords)
{
    assert(globalCoords.GetNumberOfRows() == localCoords.GetNumberOfRows());
    assert(globalCoords.GetNumberOfColumns() == localCoords.GetNumberOfColumns());
    assert(globalCoords.GetNumberOfRows() == 2);
    assert(nodes.GetNumberOfRows() == 2);
    assert(nodes.GetNumberOfColumns() == 3);

    for (int j=1; j<=localCoords.GetNumberOfColumns(); j++)
    {
        localCoords(1,j) = (nodes(1,3)*(-2.0*globalCoords(2,j)+nodes(2,1)+nodes(2,2)) + nodes(1,1)*(2.0*globalCoords(2,j)-nodes(2,2)-nodes(2,3)) - (2.0*globalCoords(1,j)-nodes(1,2))*(nodes(2,1)-nodes(2,3)))/
                           (nodes(1,3)*(nodes(2,1)-nodes(2,2)) + nodes(1,1)*(nodes(2,2)-nodes(2,3)) + nodes(1,2)*(nodes(2,3)-nodes(2,1)));

        localCoords(2,j) = (nodes(1,1)*(-2.0*globalCoords(2,j)+nodes(2,2)+nodes(2,3)) + nodes(1,2)*(2.0*globalCoords(2,j)-nodes(2,1)-nodes(2,3)) - (2.0*globalCoords(1,j)-nodes(1,3))*(nodes(2,1)-nodes(2,2)))/
                           (nodes(1,3)*(nodes(2,1)-nodes(2,2)) + nodes(1,1)*(nodes(2,2)-nodes(2,3)) + nodes(1,2)*(nodes(2,3)-nodes(2,1)));
    }
}

void Triangle::ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian)
{
    assert(Jacobian.GetNumberOfRows() == Jacobian.GetNumberOfColumns());
    assert(Jacobian.GetNumberOfRows() == 2);
    assert(nodes.GetNumberOfRows() == 2);
    assert(nodes.GetNumberOfColumns() == 3);

    Jacobian(1,1) = (nodes(1,2) - nodes(1,1))/2.0;
    Jacobian(1,2) = (nodes(1,3) - nodes(1,1))/2.0;
    Jacobian(2,1) = (nodes(2,2) - nodes(2,1))/2.0;
    Jacobian(2,2) = (nodes(2,3) - nodes(2,1))/2.0;
}


#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

Element::Element(int elementType)
{
    assert(elementType >= 1);
    mElementType = elementType;
}


void Element::MapLocalToGlobal(Vector& nodes, Vector& localCoords, Vector& globalCoords)
{
    assert(mElementType == 1);
    assert(localCoords.GetSize() == globalCoords.GetSize());

    for (int i=0; i<globalCoords.GetSize(); i++)
    {
        globalCoords[i] = nodes[0] + 0.5*(nodes[1]-nodes[0])*(localCoords[i]+1);
    }
}

void Element::MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords)
{
    assert(localCoords.GetNumberOfRows() == globalCoords.GetNumberOfRows());
    assert(localCoords.GetNumberOfColumns() == globalCoords.GetNumberOfColumns());

    switch(mElementType)
    {
        case 1:
        {
            assert(localCoords.GetNumberOfRows() == 1);
            assert(nodes.GetNumberOfRows() == 1);
            assert(nodes.GetNumberOfColumns() == 2);

            for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
            {
                globalCoords(1,j) = nodes(1,1) + 0.5*(nodes(1,2)-nodes(1,1))*(localCoords(1,j)+1);
            }
        } break;

        case 2:
        {
            assert(localCoords.GetNumberOfRows() == 2);
            assert(nodes.GetNumberOfRows() == mElementType);
            assert(nodes.GetNumberOfColumns() == 3);

            for (int i=1; i<=globalCoords.GetNumberOfRows(); i++)
            {
                for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
                {
                    globalCoords(i,j) = (1.0 - localCoords(1,j) - localCoords(2,j))*nodes(i,1) + localCoords(1,j)*nodes(i,2) + localCoords(2,j)*nodes(i,3);
                }
            }
        } break;

    }
}

void Element::MapGlobalToLocal(Vector& nodes, Vector& globalCoords, Vector& localCoords)
{
    assert(mElementType == 1);
    assert(localCoords.GetSize() == globalCoords.GetSize());

    for (int i=0; i<localCoords.GetSize(); i++)
    {
        localCoords[i] = (nodes[0] + nodes[1] - 2*globalCoords[i])/(nodes[0] - nodes[1]);
    }
}

void Element::MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords)
{
    assert(globalCoords.GetNumberOfRows() == localCoords.GetNumberOfRows());
    assert(globalCoords.GetNumberOfColumns() == localCoords.GetNumberOfColumns());

    switch(mElementType)
    {
        case 1:
        {
            assert(globalCoords.GetNumberOfRows() == 1);
            assert(nodes.GetNumberOfRows() == 1);
            assert(nodes.GetNumberOfColumns() == 2);

            for (int j=1; j<=localCoords.GetNumberOfColumns(); j++)
            {
                localCoords(1,j) = (nodes(1,1) + nodes(1,2) - 2*globalCoords(1,j))/(nodes(1,1) - nodes(1,2));
            }
        } break;

        case 2:
        {
            assert(globalCoords.GetNumberOfRows() == 2);
            assert(nodes.GetNumberOfRows() == 2);
            assert(nodes.GetNumberOfColumns() == 3);

            for (int j=1; j<=localCoords.GetNumberOfColumns(); j++)
            {
                localCoords(1,j) = (nodes(1,3)*(nodes(2,1)-globalCoords(2,j)) + nodes(1,1)*(globalCoords(2,j)-nodes(2,3)) + globalCoords(1,j)*(nodes(2,3)-nodes(2,1)))/
                                   (nodes(1,3)*(nodes(2,1)-nodes(2,2)) + nodes(1,1)*(nodes(2,2)-nodes(2,3)) + nodes(1,2)*(nodes(2,3)-nodes(2,1)));

                localCoords(2,j) = (globalCoords(2,j)*(nodes(1,1)-nodes(1,2)) + nodes(2,1)*(nodes(1,2)-globalCoords(1,j)) + nodes(2,2)*(globalCoords(1,j)-nodes(1,1)))/
                                   (nodes(2,1)*(nodes(1,2)-nodes(1,3)) + nodes(2,2)*(nodes(1,3)-nodes(1,1)) + nodes(2,3)*(nodes(1,1)-nodes(1,2)));
            }
        } break;

    }
}

void Element::GetMappingJacobian(Matrix& nodes, Matrix& Jacobian)
{
    assert(Jacobian.GetNumberOfRows() == Jacobian.GetNumberOfColumns());

    switch(mElementType)
    {
        case 1:
        {
            assert(Jacobian.GetNumberOfRows() == 1);
            assert(nodes.GetNumberOfRows() == 1);
            assert(nodes.GetNumberOfColumns() == 2);

            Jacobian(1,1) = 0.5*(nodes(1,2) - nodes(1,1));
        } break;

        case 2:
        {
            assert(Jacobian.GetNumberOfRows() == 2);
            assert(nodes.GetNumberOfRows() == 2);
            assert(nodes.GetNumberOfColumns() == 3);

            Jacobian(1,1) = nodes(1,2) - nodes(1,1);
            Jacobian(1,2) = nodes(1,3) - nodes(1,1);
            Jacobian(2,1) = nodes(2,2) - nodes(2,1);
            Jacobian(2,2) = nodes(2,3) - nodes(2,1);
        } break;

    }
}

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

    if (mElementType == 1)
    {
        assert(localCoords.GetNumberOfRows() == 1);
        assert(nodes.GetNumberOfRows() == 1);
        assert(nodes.GetNumberOfColumns() == 2);

        for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
        {
            globalCoords(1,j) = nodes(1,1) + 0.5*(nodes(1,2)-nodes(1,1))*(localCoords(1,j)+1);
        }
    }

    if (mElementType == 2)
    {
        assert(localCoords.GetNumberOfRows() == 2);
        assert(nodes.GetNumberOfRows() == 2);
        assert(nodes.GetNumberOfColumns() == 3);

        for (int i=1; i<=globalCoords.GetNumberOfRows(); i++)
        {
            for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
            {
                globalCoords(i,j) = (1.0 - localCoords(1,j) - localCoords(2,j))*nodes(i,1) + localCoords(1,j)*nodes(i,2) + localCoords(2,j)*nodes(i,3);
            }
        }
    }
}

void Element::MapGlobalToLocal(Vector& globalCoords, Vector& localCoords)
{
    if (mElementType == 1)
    {

    }

    if (mElementType == 2)
    {

    }
}

void Element::MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords)
{
    if (mElementType == 1)
    {

    }

    if (mElementType == 2)
    {

    }
}

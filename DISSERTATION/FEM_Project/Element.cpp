#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

// Default Constructor
Element::Element()
{
}

Element::Element(Vector& Connectivity)
{
    mConnectivity = new Vector(Connectivity);
}

// Copy Constructor
Element::Element(const Element& otherElement)
{
    mConnectivity = otherElement.mConnectivity;
}

// Destructor
Element::~Element()
{
    delete[] mConnectivity;
}

Element& Element::operator=(const Element& otherElement)
{
    mConnectivity = otherElement.mConnectivity;
}


void MapLocalToGlobal(Vector& localCoords, Vector& globalCoords)
{

}

void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords)
{

}

void MapGlobalToLocal(Vector& globalCoords, Vector& localCoords)
{

}

void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords)
{

}

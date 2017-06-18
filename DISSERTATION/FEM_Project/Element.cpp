#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

Element::~Element()
{
    delete mElementConnectivity;
}

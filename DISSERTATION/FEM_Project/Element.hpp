#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include "Matrix.hpp"
#include "Quadrature.hpp"
#include "Vector.hpp"

class Element
{
public:

    enum class ElementType
    {
        Interval,
        Triangle,
        Quad,
    };

    // Specialised Constructor
    Element(ElementType elementType);

    void MapLocalToGlobal(Vector& nodes, Vector& localCoords, Vector& globalCoords);
    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Vector& nodes, Vector& globalCoords, Vector& localCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Vector& nodes, Matrix& Jacobian);
    void ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian);

protected:

    ElementType mElementType;


};

#endif

#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class Element
{
protected:

    int mElementType;

public:

    // Specialised Constructor
    Element(int elementType);
    // 1 for interval, 2 for triangle,... corresponding to dimension.

    void MapLocalToGlobal(Vector& nodes, Vector& localCoords, Vector& globalCoords);
    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Vector& nodes, Vector& globalCoords, Vector& localCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);





};

#endif

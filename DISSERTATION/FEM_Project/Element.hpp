#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class Element
{
protected:

public:


    void MapLocalToGlobal(Vector& localCoords, Vector& globalCoords);
    void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Vector& globalCoords, Vector& localCoords);
    void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords);





};

#endif

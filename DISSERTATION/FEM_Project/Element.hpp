#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include "Matrix.hpp"
#include "Vector.hpp"

class Element
{
private:

    Vector* mConnectivity;


public:

    // Default Constructor
    Element();

    Element(Vector& Connectivity);

    // Copy Constructor
    Element(const Element& otherElement);

    // Destructor
    ~Element();

    Element& operator=(const Element& otherElement);


    void MapLocalToGlobal(Vector& localCoords, Vector& globalCoords);
    void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Vector& globalCoords, Vector& localCoords);
    void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords);




};

#endif

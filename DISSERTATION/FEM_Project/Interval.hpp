#ifndef INTERVALHEADERDEF
#define INTERVALHEADERDEF

#include "Element.hpp"
#include "Quadrature.hpp"

class Interval: public Element
{
private:

public:

    Interval(Vector elementConnectivity);

    Interval(Interval& otherInterval);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian);

};

#endif

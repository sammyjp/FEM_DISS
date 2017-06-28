#ifndef INTERVALHEADERDEF
#define INTERVALHEADERDEF

#include "Element.hpp"
#include "QuadratureLibrary.hpp"

class Interval: public Element
{
private:

public:

    Interval(Vector& elementConnectivity);

    Interval(Interval& otherInterval);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapLocalToGlobal(Vector& nodes, Vector& localCoords, Vector& globalCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Matrix& nodes, Vector& pointsToEval, Matrix& Jacobian);

};

#endif

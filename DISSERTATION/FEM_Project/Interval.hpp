#ifndef INTERVALHEADERDEF
#define INTERVALHEADERDEF

#include "Element.hpp"
#include "QuadratureLibrary.hpp"

class Interval: public Element
{
private:

public:

    Interval(Vector& elementConnectivity, Mesh& meshReference);

    Interval(Interval& otherInterval);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Vector& pointToEval, Matrix& Jacobian);

};

#endif

#ifndef TRIANGLEHEADERDEF
#define TRIANGLEHEADERDEF

#include "Element.hpp"
#include "QuadratureLibrary.hpp"

class Triangle: public Element
{
private:

public:

    Triangle(Vector& elementConnectivity, Mesh& meshReference);

    Triangle(const Triangle& otherTriangle);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Vector& pointToEval, Matrix& Jacobian);
};

#endif

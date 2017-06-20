#ifndef TRIANGLEHEADERDEF
#define TRIANGLEHEADERDEF

#include "Element.hpp"
#include "Quadrature.hpp"

class Triangle: public Element
{
private:

public:

    Triangle(Vector& elementConnectivity);

    Triangle(Triangle& otherTriangle);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian);

};

#endif

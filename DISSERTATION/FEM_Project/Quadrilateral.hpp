#ifndef QUADRILATERALHEADERDEF
#define QUADRILATERALHEADERDEF

#include "Element.hpp"
#include "Quadrature.hpp"

class Quadrilateral: public Element
{
private:

public:

    Quadrilateral(Vector& elementConnectivity);

    Quadrilateral(Quadrilateral& otherQuadrilateral);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian);

};

#endif

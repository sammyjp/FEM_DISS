#ifndef QUADRILATERALHEADERDEF
#define QUADRILATERALHEADERDEF

#include "Element.hpp"
#include "QuadratureLibrary.hpp"

class Quadrilateral: public Element
{
private:

public:

    Quadrilateral(Vector& elementConnectivity);

    Quadrilateral(Quadrilateral& otherQuadrilateral);

    int GetElementType() const;

    void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Matrix& nodes, Vector& pointsToEval, Matrix& Jacobian);
};

#endif

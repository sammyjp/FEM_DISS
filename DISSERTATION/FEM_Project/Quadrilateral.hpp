#ifndef QUADRILATERALHEADERDEF
#define QUADRILATERALHEADERDEF

#include "Element.hpp"
#include "QuadratureLibrary.hpp"

class Quadrilateral: public Element
{
private:

public:

    Quadrilateral(Vector& elementConnectivity, Mesh& meshReference);

    Quadrilateral(const Quadrilateral& otherQuadrilateral);

    int GetElementType() const;
    int GetNumElementNodes() const;

    void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords);
    void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords);

    void ComputeMappingJacobian(Vector& pointToEval, Matrix& Jacobian);
};

#endif

#ifndef ELEMENTHEADERDEF
#define ELEMENTHEADERDEF

#include "Matrix.hpp"
#include "QuadratureLibrary.hpp"
#include "Vector.hpp"

class Element
{
protected:

    Vector* mElementConnectivityArray;

public:

    enum ElementType
    {
        Interval,
        Triangle,
        Quadrilateral,
    };

    virtual ~Element();

    virtual int GetElementType() const = 0;
    virtual Vector GetElementConnectivityArray() const;

    virtual void MapLocalToGlobal(Matrix& nodes, Matrix& localCoords, Matrix& globalCoords) = 0;
    virtual void MapGlobalToLocal(Matrix& nodes, Matrix& globalCoords, Matrix& localCoords) = 0;
    virtual void ComputeMappingJacobian(Matrix& nodes, Matrix& Jacobian) = 0;

    virtual void ComputeElementQuadraturePoints(Vector& quadraturePoints);
    virtual double PerformElementQuadrature(Vector& quadraturePoints, Vector& functionPoints, Matrix& mappingJacobian);
};

#endif

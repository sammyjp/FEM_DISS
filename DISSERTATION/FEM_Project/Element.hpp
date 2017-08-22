#ifndef ABSTRACTELEMENTHEADERDEF
#define ABSTRACTELEMENTHEADERDEF

#include "Matrix.hpp"
#include "QuadratureLibrary.hpp"
#include "Vector.hpp"

class Mesh;

class Element
{
protected:

    Vector* mElementConnectivityArray;
    Mesh* mMeshReference;

public:

    enum ElementType
    {
        Interval,
        Triangle,
        Quadrilateral,
    };

    virtual ~Element();

    virtual int GetElementType() const = 0;
    virtual int GetNumElementNodes() const = 0;
    virtual Matrix GetElementCoordinates() const;
    virtual Vector GetElementConnectivityArray() const;

    virtual void MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords) = 0;
    virtual void MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords) = 0;
    virtual void ComputeMappingJacobian(Vector& pointToEval, Matrix& Jacobian) = 0;

    virtual void GetQuadrature(const int n_q, Vector& quadratureWeights, Matrix& localQuadraturePoints, Matrix& globalQuadraturePoints);
    virtual double PerformElementQuadrature(const int n_q, Vector& quadratureWeights, Matrix& localQuadraturePoints, Vector& functionPoints);
};

#endif

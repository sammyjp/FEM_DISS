#include <cmath>
#include <cassert>

#include "Quadrilateral.hpp"

Quadrilateral::Quadrilateral(Vector& elementConnectivity, Mesh& meshReference)
{
    mElementConnectivityArray = new Vector(elementConnectivity);
    mMeshReference = &meshReference;
}

Quadrilateral::Quadrilateral(const Quadrilateral& otherQuadrilateral)
{
    mElementConnectivityArray = new Vector (*otherQuadrilateral.mElementConnectivityArray);
}

int Quadrilateral::GetElementType() const
{
    return ElementType::Quadrilateral;
}

void Quadrilateral::MapLocalToGlobal(Matrix& localCoords, Matrix& globalCoords)
{
    Matrix* nodes = new Matrix (GetElementCoordinates());
    assert(localCoords.GetNumberOfRows() == globalCoords.GetNumberOfRows());
    assert(localCoords.GetNumberOfColumns() == globalCoords.GetNumberOfColumns());
    assert(localCoords.GetNumberOfRows() == 2);
    assert(nodes->GetNumberOfRows() == 2);
    assert(nodes->GetNumberOfColumns() == 4);

    for (int i=1; i<=globalCoords.GetNumberOfRows(); i++)
    {
        for (int j=1; j<=globalCoords.GetNumberOfColumns(); j++)
        {
            globalCoords(i,j) = 0.25*((*nodes)(i,1)*((1-localCoords(1,j))*(1-localCoords(2,j))) + (*nodes)(i,2)*((1+localCoords(1,j))*(1-localCoords(2,j)))+
                                (*nodes)(i,3)*((1+localCoords(1,j))*(1+localCoords(2,j))) + (*nodes)(i,4)*((1-localCoords(1,j))*(1+localCoords(2,j))));
        }
    }

    delete nodes;
}

void Quadrilateral::MapGlobalToLocal(Matrix& globalCoords, Matrix& localCoords)
{

}

void Quadrilateral::ComputeMappingJacobian(Vector& pointToEval, Matrix& Jacobian)
{
    Matrix* nodes = new Matrix (GetElementCoordinates());
    assert(Jacobian.GetNumberOfRows() == Jacobian.GetNumberOfColumns());
    assert(Jacobian.GetNumberOfRows() == 2);
    assert(nodes->GetNumberOfRows() == 2);
    assert(nodes->GetNumberOfColumns() == 4);

    Jacobian(1,1) = 0.25*(-(*nodes)(1,1)*(1-pointToEval(2)) + (*nodes)(1,2)*(1-pointToEval(2)) + (*nodes)(1,3)*(1+pointToEval(2)) - (*nodes)(1,4)*(1+pointToEval(2)));
    Jacobian(1,2) = 0.25*(-(*nodes)(1,1)*(1-pointToEval(1)) - (*nodes)(1,2)*(1+pointToEval(1)) + (*nodes)(1,3)*(1+pointToEval(1)) + (*nodes)(1,4)*(1-pointToEval(1)));
    Jacobian(2,1) = 0.25*(-(*nodes)(2,1)*(1-pointToEval(2)) + (*nodes)(2,2)*(1-pointToEval(2)) + (*nodes)(2,3)*(1+pointToEval(2)) - (*nodes)(2,4)*(1+pointToEval(2)));
    Jacobian(2,2) = 0.25*(-(*nodes)(2,1)*(1-pointToEval(1)) - (*nodes)(2,2)*(1+pointToEval(1)) + (*nodes)(2,3)*(1+pointToEval(1)) + (*nodes)(2,4)*(1-pointToEval(1)));

    delete nodes;
}

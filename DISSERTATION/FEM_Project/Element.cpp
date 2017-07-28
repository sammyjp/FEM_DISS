#include <cmath>
#include <cassert>

#include "Element.hpp"
#include "Mesh.hpp"

Element::~Element()
{
    delete mElementConnectivityArray;
}

Matrix Element::GetElementCoordinates() const
{
    return mMeshReference->GetGridPoints((*mElementConnectivityArray));
}

Vector Element::GetElementConnectivityArray() const
{
    return *mElementConnectivityArray;
}

void Element::GetQuadrature(const int n_q, Vector& quadratureWeights, Matrix& localQuadraturePoints, Matrix& globalQuadraturePoints)
{
    QuadratureLibrary().GetQuadrature(GetElementType(), n_q, quadratureWeights, localQuadraturePoints);
    MapLocalToGlobal(localQuadraturePoints, globalQuadraturePoints);
}

double Element::PerformElementQuadrature(const int n_q, Vector& quadratureWeights, Matrix& localQuadraturePoints, Vector& functionPoints)
{
    double I = 0;

    Matrix* jacobian = new Matrix (mMeshReference->GetDimension(), mMeshReference->GetDimension());

    for (int i=1; i<=n_q; i++)
    {
        Vector* pointToEvalJacobian = new Vector (localQuadraturePoints.GetColumnAsVector(i));
        ComputeMappingJacobian(*pointToEvalJacobian, *jacobian);

        I += quadratureWeights(i)*functionPoints(i)*(jacobian->CalculateDeterminant());

        delete pointToEvalJacobian;
    }

    delete jacobian;

    return I;
}

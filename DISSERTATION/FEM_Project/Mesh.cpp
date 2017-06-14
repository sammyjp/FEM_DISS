#include <cmath>
#include <cassert>
#include "Mesh.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

// Specialised Constructor
Mesh::Mesh(Matrix& gridPoints, int numElements, Matrix& connectivity)
{
    assert(numElements > 0);

    mDimension = gridPoints.GetNumberOfRows();
    mNumNodes = gridPoints.GetNumberOfRows();
    mNumElements = numElements;

    mGridPoints = new Matrix(gridPoints);
    mConnectivity = new Matrix(connectivity);
    mElement = new Element();
}


// Copy Constructor
Mesh::Mesh(const Mesh& otherMesh)
{
    mDimension = otherMesh.mDimension;
    mNumElements = otherMesh.mNumElements;
    mNumNodes = otherMesh.mNumNodes;

    mGridPoints = otherMesh.mGridPoints;
    mConnectivity = otherMesh.mConnectivity;
}

// Destructor
Mesh::~Mesh()
{
    delete mGridPoints;
    delete mConnectivity;
    delete mElement;
}

int Mesh::GetDimension() const
{
    return mDimension;
}

int Mesh::GetNumElements() const
{
    return mNumElements;
}

int Mesh::GetNumNodes() const
{
    return mNumNodes;
}

Matrix Mesh::GetGridPoints() const
{
    return *mGridPoints;
}

Matrix Mesh::GetConnectivityArray() const
{
    return *mConnectivity;
}

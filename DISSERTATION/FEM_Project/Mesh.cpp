#include <cmath>
#include <cassert>
#include "Mesh.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

// Specialised Constructor
Mesh::Mesh(int dimension, int numElements)
{
    assert(numElements > 0);
    assert(dimension > 0);

    mDimension = dimension;
    mNumElements = numElements;

    //
    mNumNodes = numElements + 1;
    //

    mGridPoints = new Matrix(mDimension, mNumNodes);

    mElementsArray = new Element [mNumElements];
    for (int i=0; i<mNumElements; i++)
    {
        mElementsArray[i] = new Element();
    }
}


// Copy Constructor
Mesh::Mesh(const Mesh& otherMesh)
{
    mDimension = otherMesh.mDimension;
    mNumElements = otherMesh.mNumElements;
    mNumNodes = otherMesh.mNumNodes;
    mNumEdges = otherMesh.mNumEdges;
    mNumFaces = otherMesh.mNumFaces;

    mGridPoints = otherMesh.mGridPoints;
    mElementsArray = otherMesh.mElementsArray;
}

// Destructor
Mesh::~Mesh()
{
    delete mGridPoints;
    delete mElementsArray;
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

int Mesh::GetNumEdges() const
{
    return mNumEdges;
}

int Mesh::GetNumFaces() const
{
    return mNumFaces;
}

Matrix Mesh::GetGridPoints() const
{
    return *mGridPoints;
}

#include <cmath>
#include <cassert>
#include "Mesh.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

// Specialised Constructor
Mesh::Mesh(int numElements, int dimension)
{
    assert(numElements > 0);
    assert(dimension > 0);

    mNumElements = numElements;
    mDimension = dimension;

    //
    mNumNodes = numElements + 1;
    //

    if (dimension == 1)
    {
        mNumNodes = numElements + 1;
        mXGridPoints = new Vector(mNumNodes);
    }
}


// Copy Constructor
Mesh::Mesh(const Mesh& otherMesh)
{
    mNumElements = otherMesh.mNumElements;
    mNumNodes = otherMesh.mNumNodes;
    mNumEdges = otherMesh.mNumEdges;
    mNumFaces = otherMesh.mNumFaces;
    mDimension = otherMesh.mDimension;

    mXGridPoints = otherMesh.mXGridPoints;
    mYGridPoints = otherMesh.mYGridPoints;
    mZGridPoints = otherMesh.mZGridPoints;
}

// Destructor
Mesh::~Mesh()
{
    delete mXGridPoints;
    if (mDimension > 1)
    {
        delete mYGridPoints;
    }
    if (mDimension > 2)
    {
        delete mZGridPoints;
    }
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

int Mesh::GetDimension() const
{
    return mDimension;
}

Vector Mesh::GetXGridPoints() const
{
    return *mXGridPoints;
}

Vector Mesh::GetYGridPoints() const
{
    assert(mDimension > 1);
    return *mYGridPoints;
}

Vector Mesh::GetZGridPoints() const
{
    assert(mDimension > 2);
    return *mZGridPoints;
}

void Mesh::GenerateUniformMesh()
{
    mStepSize = 1.0/(mNumNodes - 1);
    for (int i=0; i<mXGridPoints->GetSize(); i++)
    {
        (*mXGridPoints)[i] = i*mStepSize;
    }
}

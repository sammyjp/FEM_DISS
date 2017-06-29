#include <cmath>
#include <cassert>

#include "Mesh.hpp"

// Specialised Constructor
Mesh::Mesh(Matrix& gridPoints, int numElements)
{
    mDimension = gridPoints.GetNumberOfRows();
    mNumNodes = gridPoints.GetNumberOfColumns();
    mNumElements = numElements;

    mGridPoints = new Matrix(gridPoints);
    mElementsArray = new Element* [mNumElements];
}

// Copy Constructor
Mesh::Mesh(const Mesh& otherMesh)
{
    mDimension = otherMesh.mDimension;
    mNumElements = otherMesh.mNumElements;
    mNumNodes = otherMesh.mNumNodes;

    mGridPoints = new Matrix (*otherMesh.mGridPoints);
    mElementsArray = new Element* [mNumElements];

    for (int i=0; i<mNumElements; i++)
    {
        mElementsArray[i] = otherMesh.mElementsArray[i];
    }
}

// Destructor
Mesh::~Mesh()
{
    for (int i=0; i<mNumElements; i++)
    {
        delete mElementsArray[i];
    }
    delete[] mElementsArray;

    delete mGridPoints;
}

void Mesh::InitialiseElement(int elementNumber, Vector& connectvity, int elementType)
{
    switch(elementType)
    {
    case 0:
        {
            mElementsArray[elementNumber - 1] = new Interval(connectvity, *this);
        } break;
    case 1:
        {
            mElementsArray[elementNumber - 1] = new Triangle(connectvity, *this);
        } break;
    case 2:
        {
            mElementsArray[elementNumber - 1] = new Quadrilateral(connectvity, *this);
        } break;
    }
}

Element* Mesh::GetElement(int elementArrayIndex) const
{
    assert(elementArrayIndex > 0 && elementArrayIndex <= mNumElements);
    return mElementsArray[elementArrayIndex-1];
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

Matrix Mesh::GetAllGridPoints() const
{
    return *mGridPoints;
}

Matrix Mesh::GetGridPoints(Vector& elementConnectivityArray) const
{
    Matrix mat(mDimension, elementConnectivityArray.GetSize());

    for (int i=1; i<=mat.GetNumberOfRows(); i++)
    {
        for (int j=1; j<=mat.GetNumberOfColumns(); j++)
        {
            mat(i,j) = (*mGridPoints)(i,elementConnectivityArray[j-1]);
        }
    }
    return mat;
}

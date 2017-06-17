#include <cmath>
#include <cassert>
#include "Mesh.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

// Specialised Constructor
Mesh::Mesh(Matrix& gridPoints, Matrix& connectivity)
{
    mDimension = gridPoints.GetNumberOfRows();
    mNumNodes = gridPoints.GetNumberOfRows();
    mNumElements = connectivity.GetNumberOfRows();

    mGridPoints = new Matrix(gridPoints);
    mConnectivity = new Matrix(connectivity);
    mElementsArray = new Element* [mNumElements];

    switch(mDimension)
    {
        case 1:
        {
            for (int i=0; i<mNumElements; i++)
            {
                mElementsArray[i] = new Element(Element::ElementType::Interval);
            }
        } break;

        case 2:
        {
            for (int i=1; i<=connectivity.GetNumberOfRows(); i++)
            {
                int nonZeros = 0;
                for (int j=1; j<=connectivity.GetNumberOfColumns(); j++)
                {
                    if (connectivity(i,j) != 0)
                    {
                        nonZeros++;
                    }
                }
                switch(nonZeros)
                {
                    case 3:
                    {
                        mElementsArray[i] = new Element(Element::ElementType::Triangle);
                    }
                    case 4:
                    {
                        mElementsArray[i] = new Element(Element::ElementType::Quad);
                    }
                }
            }
        } break;
    }


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
    for (int i=0; i<mNumElements; i++)
    {
        delete[] mElementsArray[i];
    }
    delete[] mElementsArray;

    delete mGridPoints;
    delete mConnectivity;
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

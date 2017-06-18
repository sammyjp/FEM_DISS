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
            for (int i=1; i<=mNumElements; i++)
            {
                Vector* elementConnectivity =  new Vector(mConnectivity->GetRowAsVector(i));
                mElementsArray[i-1] = new Interval(*elementConnectivity);
                delete elementConnectivity;
            }
        } break;

        case 2:
        {
            for (int i=1; i<=mNumElements; i++)
            {
                int nonZeros = 0;
                for (int j=1; j<=connectivity.GetNumberOfColumns(); j++)
                {
                    if (connectivity(i,j) != 0)
                    {
                        nonZeros++;
                    }
                }

                Vector* elementConnectivity =  new Vector(mConnectivity->GetRowAsVector(i));

                switch(nonZeros)
                {
                    case 3:
                    {
                        mElementsArray[i-1] = new Triangle(*elementConnectivity);
                    } break;

                    case 4:
                    {
                        mElementsArray[i-1] = new Quadrilateral(*elementConnectivity);
                    } break;
                }
                delete elementConnectivity;
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
        delete mElementsArray[i];
    }
    delete[] mElementsArray;

    delete mGridPoints;
    delete mConnectivity;
}

Element* Mesh::GetElement(int elementArrayIndex) const
{
    assert(elementArrayIndex >= 0 && elementArrayIndex < mNumElements);
    return mElementsArray[elementArrayIndex];
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

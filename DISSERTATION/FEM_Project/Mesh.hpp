#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include "Element.hpp"
#include "Interval.hpp"
#include "Triangle.hpp"
#include "Quadrilateral.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"


class Mesh
{
private:

    int mDimension;
    int mNumElements;
    int mNumNodes;

    Matrix* mGridPoints;
    Element** mElementsArray;

public:

    // Specialised Constructor
    Mesh(Matrix& gridPoints, int numElements);

    // Copy Constructor
    Mesh(const Mesh& otherMesh);

    // Destructor
    ~Mesh();

    void InitialiseElement(int elementNumber, Vector& connectvity, int elementType);

    int GetDimension() const;
    int GetNumElements() const;
    int GetNumNodes() const;

    Element* GetElement(int elementArrayIndex) const;

    Matrix GetAllGridPoints() const;
    Matrix GetGridPoints(Vector& elementConnectivityArray) const;
};

#endif

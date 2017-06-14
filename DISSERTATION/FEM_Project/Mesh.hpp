#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include "Element.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

class Mesh
{
private:

    int mDimension;
    int mNumElements;
    int mNumNodes;

    Matrix* mGridPoints;
    Matrix* mConnectivity;
    Element* mElement;

public:

    // Specialised Constructor
    Mesh(Matrix& gridPoints, int numElements, Matrix& connectivity);

    // Copy Constructor
    Mesh(const Mesh& otherMesh);

    // Destructor
    ~Mesh();

    int GetDimension() const;
    int GetNumElements() const;
    int GetNumNodes() const;

    Matrix GetGridPoints() const;
    Matrix GetConnectivityArray() const;
};

#endif

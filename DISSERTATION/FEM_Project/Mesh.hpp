#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include "Element.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

class Mesh
{
private:

    int mNumElements;
    int mNumNodes;
    int mNumEdges;
    int mNumFaces;
    int mDimension;

    Matrix* mGridPoints;
    Element* mElementsArray;

public:

    // Specialised Constructor
    Mesh(int dimension, int numElements);

    // Copy Constructor
    Mesh(const Mesh& otherMesh);

    // Destructor
    ~Mesh();

    int GetNumElements() const;
    int GetNumNodes() const;
    int GetNumEdges() const;
    int GetNumFaces() const;
    int GetDimension() const;

    Matrix GetGridPoints() const;
};

#endif

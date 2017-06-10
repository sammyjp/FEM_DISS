#ifndef MESHHEADERDEF
#define MESHHEADERDEF

#include "Vector.hpp"
#include "Element.hpp"

class Mesh
{
private:

    int mNumElements;
    int mNumNodes;
    int mNumEdges;
    int mNumFaces;
    int mDimension;

    double mStepSize;

    Vector* mXGridPoints;
    Vector* mYGridPoints;
    Vector* mZGridPoints;

public:

    // Specialised Constructor
    Mesh(int numElements, int dimension);

    // Copy Constructor
    Mesh(const Mesh& otherMesh);

    // Destructor
    ~Mesh();

    int GetNumElements() const;
    int GetNumNodes() const;
    int GetNumEdges() const;
    int GetNumFaces() const;
    int GetDimension() const;

    Vector GetXGridPoints() const;
    Vector GetYGridPoints() const;
    Vector GetZGridPoints() const;

    void GenerateUniformMesh();
};

#endif

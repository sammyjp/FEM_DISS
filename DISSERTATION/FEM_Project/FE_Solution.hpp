#ifndef FE_SOLUTIONHEADERDEF
#define FE_SOLUTIONHEADERDEF

#include "Mesh.hpp"
#include "PolynomialSpace.hpp"
#include "QuadratureLibrary.hpp"
#include "Vector.hpp"

class FE_Solution
{
private:

    Mesh* mMesh;
    int mPolynomialDegree;

    Vector* dofStart;

    PolynomialSpace** mElementPolySpaceArray;

    void InitialiseElementDofs();

public:

    FE_Solution(Mesh& mesh, int polynomialDegree);

    ~FE_Solution();

    PolynomialSpace* GetElementPolynomialSpace(int elementNumber) const;


};

#endif

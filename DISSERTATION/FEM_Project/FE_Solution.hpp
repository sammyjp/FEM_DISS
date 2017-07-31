#ifndef FE_SOLUTIONHEADERDEF
#define FE_SOLUTIONHEADERDEF

#include "Mesh.hpp"
#include "PolynomialSpace.hpp"
#include "QuadratureLibrary.hpp"
#include "Vector.hpp"

class FE_Solution
{
private:

    Mesh* mMeshReference;
    int mPolynomialDegree;
    int mNumDofs;
    Vector* mSolutionVector;

    Vector* dofStart;

    PolynomialSpace** mElementPolySpaceArray;

    void InitialiseElementDofs();

public:

    FE_Solution(Mesh& mesh, int polynomialDegree);

    ~FE_Solution();

    PolynomialSpace* GetElementPolynomialSpace(int elementNumber) const;

    Vector GetElementDofs(int elementNumber);
    int GetNumberOfDofs();
    int GetNumElementDofs(int elementNumber);

    void SetSolutionVector();
    Vector& GetSolutionVector();

    void ComputeBasis(int elementNumber, double localGridPoint, Vector& basisValues);
    void ComputeBasis(int elementNumber, Vector& localGridPoint, Vector& basisValues);

    void ComputeGradBasis(int elementNumber, double localGridPoint, Matrix& gradBasisValues);
    void ComputeGradBasis(int elementNumber, Vector& localGridPoint, Matrix& gradBasisValues);
};

#endif

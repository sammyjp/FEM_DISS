#ifndef FE_SOLUTIONHEADERDEF
#define FE_SOLUTIONHEADERDEF

#include "Mesh.hpp"
#include "Vector.hpp"

class FE_Solution
{
private:

    Mesh* mMesh;
    int mPolynomialDegree;

public:

    FE_Solution(Mesh& mesh, int polynomialDegree);

    ~FE_Solution();

    void ComputeBasisFunctionValues(int dofNumber, Vector& functionValues);
    void ComputeBasisFunctionGrad(int dofNumber, Matrix& gradValues);

    void ComputeLinearBasisFunctionValues(int i, Vector& functionValues, Matrix& x);

    void ComputeLinearBasisFunctionDerivativeValues(int i, Vector& derivativeValues, Matrix& x);

};

#endif

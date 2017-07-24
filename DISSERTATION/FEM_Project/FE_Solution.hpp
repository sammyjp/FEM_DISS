#ifndef FE_SOLUTIONHEADERDEF
#define FE_SOLUTIONHEADERDEF

#include "Mesh.hpp"
#include "Vector.hpp"

class FE_Solution
{
private:

    Mesh* mMesh;

public:

    FE_Solution(Mesh& mesh);

    void ComputeLinearBasisFunctionValues(int i, Vector& functionValues, Matrix& x);

    void ComputeLinearBasisFunctionDerivativeValues(int i, Vector& derivativeValues, Matrix& x);

};

#endif

#ifndef FE_SOLUTIONHEADERDEF
#define FE_SOLUTIONHEADERDEF

#include "Mesh.hpp"
#include "Vector.hpp"

class FE_Solution
{
private:

    Mesh* mMesh;


public:

    // Specialised Constructor
    FE_Solution(Mesh& mesh);

    // Destructor
    ~FE_Solution();

    void ComputeLinearBasisFunctionValues(int i, Vector& functionValues, Vector& x);



};

#endif

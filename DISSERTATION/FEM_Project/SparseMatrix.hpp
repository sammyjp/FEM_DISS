#ifndef SPARSEMATRIXHEADERDEF
#define SPARSEMATRIXHEADERDEF

#include "FE_Solution.hpp"
#include "Matrix.hpp"
#include "Vector.hpp"

class SparseMatrix
{
private:
    Vector* mData;
    Vector* mRow_ptr;
    Vector* mCol_index;

    int mNumNonZeros;
    int mNumRows;
    int mNumCols;

    double ReadDataArray(int i) const;
    int ReadRowPointerArray(int i) const;
    int ReadColumnIndexArray(int i) const;

public:

    SparseMatrix(FE_Solution& FE, int numElements);

    SparseMatrix(SparseMatrix& otherSparseMatrix);

    ~SparseMatrix();

    double Read(int i, int j) const;

    Vector GetDataArray() const;
    Vector GetRowPointerArray() const;
    Vector GetColumnIndexArray() const;
    int GetNumberOfRows() const;
    int GetNumberOfColumns() const;

    void AddValue(double value, int i, int j);

    friend std::ostream& operator<<(std::ostream& output, const SparseMatrix& m);

    friend Vector operator*(const SparseMatrix& m, const Vector& v);

    void CGSolveSystem(const Vector& rightHandVector, Vector& solutionVector, double tolerance, int maxIterations);

};

#endif

#ifndef SPARSEMATRIXHEADERDEF
#define SPARSEMATRIXHEADERDEF

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

public:

    SparseMatrix(Matrix& otherMatrix);

    SparseMatrix(SparseMatrix& otherSparseMatrix);

    ~SparseMatrix();

    Vector GetDataArray() const;
    Vector GetRowPointerArray() const;
    Vector GetColumnIndexArray() const;
    int GetNumberOfRows() const;
    int GetNumberOfColumns() const;

    double Read(int i, int j) const;

    double ReadDataArray(int i) const;
    int ReadRowPointerArray(int i) const;
    int ReadColumnIndexArray(int i) const;

    friend std::ostream& operator<<(std::ostream& output, const SparseMatrix& m);

    friend Vector operator*(const SparseMatrix& m, const Vector& v);

    void CGSolveSystem(const Vector& rightHandVector, Vector& solutionVector, double tolerance, int maxIterations);

};

#endif

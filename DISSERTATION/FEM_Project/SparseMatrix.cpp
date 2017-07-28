#include <cmath>
#include <cassert>
#include <iomanip>

#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(FE_Solution& FE, int numElements)
{
    mNumRows = FE.GetNumberOfDofs();
    mNumCols = FE.GetNumberOfDofs();
    mNumNonZeros = 0;

    for (int k=1; k<=numElements; k++)
    {
        FE.GetElementDofs(k);

        for (int i=1; i<=FE.GetElementPolynomialSpace(k)->GetNumElementDofs(); i++)
        {
            for (int j=1; j<=FE.GetElementPolynomialSpace(k)->GetNumElementDofs(); j++)
            {
            }
        }
    }
}

SparseMatrix::SparseMatrix(Matrix& otherMatrix)
{
    mNumRows = otherMatrix.GetNumberOfRows();
    mNumCols = otherMatrix.GetNumberOfColumns();

    for (int i=1; i<=mNumRows; i++)
    {
        for (int j=1; j<=mNumCols; j++)
        {
            if (otherMatrix(i,j) != 0)
            {
                mNumNonZeros++;
            }
        }
    }
    assert(mNumNonZeros != 0);

    mData = new Vector (mNumNonZeros);
    mRow_ptr = new Vector (mNumRows+1);
    mCol_index = new Vector (mNumNonZeros);

    int dataIndex = 0;

    for (int i=1; i<=mNumRows; i++)
    {
        (*mRow_ptr)[i-1] = dataIndex;

        for (int j=1; j<=mNumCols; j++)
        {
            if (otherMatrix(i,j) != 0)
            {
                (*mData)[dataIndex] = otherMatrix(i,j);
                (*mCol_index)[dataIndex] = j-1;
                dataIndex++;
            }
        }
    }

    (*mRow_ptr)[otherMatrix.GetNumberOfRows()] = mNumNonZeros;
}

SparseMatrix::SparseMatrix(SparseMatrix& otherSparseMatrix)
{
    mNumRows = otherSparseMatrix.mNumRows;
    mNumCols = otherSparseMatrix.mNumCols;
    mNumNonZeros = otherSparseMatrix.mNumNonZeros;

    mData = new Vector (mNumNonZeros);
    mRow_ptr = new Vector (mNumRows+1);
    mCol_index = new Vector (mNumNonZeros);

    for (int i=0; i<mNumNonZeros; i++)
    {
        mData[i] = otherSparseMatrix.mData[i];
        mCol_index[i] = otherSparseMatrix.mCol_index[i];
    }
    for (int i=0; i<mNumRows+1; i++)
    {
        mRow_ptr[i] = otherSparseMatrix.mRow_ptr[i];
    }
}

SparseMatrix::~SparseMatrix()
{
    delete mData;
    delete mRow_ptr;
    delete mCol_index;
}

double SparseMatrix::Read(int i, int j) const
{
    assert(i > 0);
    assert(i < mNumRows+1);
    assert(j > 0);
    assert(j < mNumCols+1);

    for (int k=(*mRow_ptr)[i-1]; k<(*mRow_ptr)[i]; k++)
    {
        if((*mCol_index)[k]==j-1)
        {
            return (*mData)[k];
        }
    }
    return 0.0;
}

double SparseMatrix::ReadDataArray(int i) const
{
    return (*mData)[i-1];
}

int SparseMatrix::ReadRowPointerArray(int i) const
{
    return (*mRow_ptr)[i-1];
}

int SparseMatrix::ReadColumnIndexArray(int i) const
{
    return (*mCol_index)[i-1];
}

Vector SparseMatrix::GetDataArray() const
{
    return *mData;
}

Vector SparseMatrix::GetRowPointerArray() const
{
    return *mRow_ptr;
}

Vector SparseMatrix::GetColumnIndexArray() const
{
    return *mCol_index;
}

int SparseMatrix::GetNumberOfRows() const
{
    return mNumRows;
}

int SparseMatrix::GetNumberOfColumns() const
{
    return mNumCols;
}

std::ostream& operator<<(std::ostream& output, const SparseMatrix& m)
{
  for (int i=1; i<=m.GetNumberOfRows(); i++)
  {
    for (int j=1; j<=m.GetNumberOfColumns(); j++)
    {
      output << std::setw(14)
             << std::setprecision(5)
	     << std::scientific
	     << m.Read(i,j);
    }
    output << std::endl;
  }
  output << std::endl;

  return output;

}

Vector operator*(const SparseMatrix& m, const Vector& v)
{
    int original_vector_size = v.GetSize();
    assert(m.GetNumberOfColumns() == original_vector_size);
    int new_vector_length = m.GetNumberOfRows();
    Vector new_vector(new_vector_length);

    for (int i=0; i<new_vector_length; i++)
    {
        for (int j=m.ReadRowPointerArray(i+1); j<m.ReadRowPointerArray(i+2); j++)
        {
            new_vector[i] += m.ReadDataArray(j+1)*v.Read(m.ReadColumnIndexArray(j+1));
        }
    }

    return new_vector;
}

void SparseMatrix::CGSolveSystem(const Vector& rightHandVector, Vector& solutionVector, double tolerance, int maxIterations)
{
    assert(mNumCols == solutionVector.GetSize());
    assert(mNumRows == rightHandVector.GetSize());

    int iterationCounter = 0;
    Vector* r = new Vector (rightHandVector.GetSize());
    Vector* r_prev = new Vector (rightHandVector.GetSize());
    Vector* p = new Vector (r->GetSize());
    Vector* w = new Vector (mNumRows);
    double alpha;
    double beta;

    (*r) = rightHandVector - ((*this)*solutionVector);
    (*p) = (*r);

    do
    {
        (*w) = (*this)*(*p);
        alpha = (r->ScalarProduct((*r)))/(p->ScalarProduct((*w)));
        solutionVector = solutionVector + ((*p)*alpha);
        (*r_prev) = (*r);
        (*r) = (*r_prev) - ((*w)*alpha);
        beta = (r->ScalarProduct((*r)))/(r_prev->ScalarProduct((*r_prev)));
        (*p) = (*r) + ((*p)*beta);

        iterationCounter++;
    } while (r->CalculateNorm(2) > tolerance && iterationCounter <= maxIterations);

    delete r;
    delete r_prev;
    delete p;
    delete w;
}

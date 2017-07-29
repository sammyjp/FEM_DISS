#include <cmath>
#include <cassert>
#include <iomanip>

#include "SparseMatrix.hpp"

SparseMatrix::SparseMatrix(FE_Solution& FE, int numElements)
{
    mNumRows = FE.GetNumberOfDofs();
    mNumCols = FE.GetNumberOfDofs();

    mRow_ptr = new Vector (mNumRows+1);
    (*mRow_ptr) = 1;

    Vector* NNZ = new Vector (FE.GetNumberOfDofs());

    mNumNonZeros = 0;

    for (int k=1; k<=numElements; k++)
    {
        for (int i=1; i<=FE.GetElementPolynomialSpace(k)->GetNumElementDofs(); i++)
        {
            for (int j=1; j<=FE.GetElementPolynomialSpace(k)->GetNumElementDofs(); j++)
            {
                (*NNZ)((FE.GetElementDofs(k))(i)) ++;
            }
        }
    }

    for (int i=0; i<NNZ->GetSize(); i++)
    {
        mNumNonZeros += (*NNZ)[i];
    }

    for (int i=1; i<mRow_ptr->GetSize(); i++)
    {
        (*mRow_ptr)[i] = (*mRow_ptr)[i-1] + (*NNZ)[i-1];
    }

    delete NNZ;

    mData = new Vector (mNumNonZeros);
    mCol_index = new Vector (mNumNonZeros);

    for (int k=1; k<=numElements; k++)
    {
        for (int i=1; i<=FE.GetElementPolynomialSpace(k)->GetNumElementDofs(); i++)
        {
            for (int j=1; j<=FE.GetElementPolynomialSpace(k)->GetNumElementDofs(); j++)
            {
                for (int l=0; l<(*mRow_ptr)((FE.GetElementDofs(k))(i)+1) - (*mRow_ptr)((FE.GetElementDofs(k))(i)); l++)
                {
                    if ((*mCol_index)((*mRow_ptr)((FE.GetElementDofs(k))(i)) + l) == 0)
                    {
                        (*mCol_index)((*mRow_ptr)((FE.GetElementDofs(k))(i)) + l) = (FE.GetElementDofs(k))(j);
                        break;
                    }
                }
            }
        }
    }
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

    for (int k=(*mRow_ptr)(i); k<(*mRow_ptr)(i+1); k++)
    {
        if((*mCol_index)(k)==j)
        {
            return (*mData)(k);
        }
    }
    return 0.0;
}

double SparseMatrix::ReadDataArray(int i) const
{
    assert(i > 0 && i <= mData->GetSize());
    return (*mData)(i);
}

int SparseMatrix::ReadRowPointerArray(int i) const
{
    assert(i > 0 && i <= mRow_ptr->GetSize());
    return (*mRow_ptr)(i);
}

int SparseMatrix::ReadColumnIndexArray(int i) const
{
    assert(i > 0 && i <= mCol_index->GetSize());
    return (*mCol_index)(i);
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

void SparseMatrix::AddValue(double value, int i, int j)
{
    assert(i > 0);
    assert(i < mNumRows+1);
    assert(j > 0);
    assert(j < mNumCols+1);

    for (int k=(*mRow_ptr)(i); k<(*mRow_ptr)(i+1); k++)
    {
        if((*mCol_index)(k)==j)
        {
            (*mData)(k) += value;
            return;
        }
    }
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

    for (int i=1; i<=new_vector_length; i++)
    {
        for (int j=m.ReadRowPointerArray(i); j<m.ReadRowPointerArray(i+1); j++)
        {
            std::cout << "i = " << i << ". j = " << j << ".\n";
            std::cout << "Column index = " << m.ReadColumnIndexArray(j) << ".\n";
            new_vector(i) += m.ReadDataArray(j)*v.Read(m.ReadColumnIndexArray(j));
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

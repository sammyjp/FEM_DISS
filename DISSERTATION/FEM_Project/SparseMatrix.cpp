#include <cmath>
#include <cassert>

#include "SparseMatrix.hpp"

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

Vector operator*(const Vector& v, const SparseMatrix& m)
{
    int original_vector_size = v.GetSize();
    assert(m.GetNumberOfRows() == original_vector_size);
    int new_vector_length = m.GetNumberOfColumns();
    Vector new_vector(new_vector_length);

    ////TODO////

//    for (int i=0; i<new_vector_length; i++)
//    {
//        for (int j=0; j<original_vector_size; j++)
//        {
//            new_vector[i] += v.Read(j)*m.mData[j][i];
//        }
//    }

    return new_vector;
}

#ifndef MATRIXHEADERDEF
#define MATRIXHEADERDEF
#include "Vector.hpp"

// Class adapted from Whiteley, Pitt-Francis
class Matrix
{

private:

  // Entries of matrix
  double** mData;

  // Dimensions
  int mNumRows, mNumCols;

public:

  // Copy constructor
  Matrix(const Matrix& otherMatrix);

  // Specialised Constructor
  Matrix(int numRows, int numCols);

  // Destructor
  ~Matrix();

  // Accessors
  int GetNumberOfRows() const;
  int GetNumberOfColumns() const;

  // 1-based indexing
  double& operator()(int i, int j);

  //Overloaded assignment operator
  Matrix& operator=(const Matrix& otherMatrix);

  // Unary +
  Matrix operator+() const;

  // Unary -
  Matrix operator-() const;

  // Binary +
  Matrix operator+(const Matrix& m1) const;

  // Binary -
  Matrix operator-(const Matrix& m1) const;

  // Scalar multiplication
  Matrix operator*(double a) const;

  // Determinant
  double CalculateDeterminant() const;

  // Declare vector multiplication friendship
  friend Vector operator*(const Matrix& m,
                          const Vector& v);
  friend Vector operator*(const Vector& v,
                          const Matrix& m);

  // Overridden << operator
  friend std::ostream& operator<<(std::ostream& output,
                        const Matrix& m);

  // Matrix multiplication
  Matrix operator*(const Matrix& m1) const;

  // Kronecker Product
  Matrix KroneckerProduct(const Matrix& m1) const;

  // Sets the matrix to the identity matrix
  void MakeIdentity();

  // Computes transpose
  Matrix Transpose();

};

// prototype signatures for friend operators
Vector operator*(const Matrix& m, const Vector& v);
Vector operator*(const Vector& v, const Matrix& m);

#endif

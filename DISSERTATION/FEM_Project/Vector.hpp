#ifndef VECTORHEADERDEF
#define VECTORHEADERDEF

#include <iostream>

// Class adapted from Whiteley, Pitt-Francis
class Vector
{

private:

  // Data stored in vector
  double* mData;

  // Size of vector
  int mSize;

public:

  // Copy constructor
  Vector(const Vector& otherVector);

  // Specialised constructor
  Vector(int size);

  // Destructor
  ~Vector();

  // Accessor
  int GetSize() const;

  // Zero-based indexing
  double& operator[](int i);

  // Read-only zero-based indexing
  double Read(int i) const;

  // One-based indexing
  double& operator()(int i);

  // Assignment
  Vector& operator=(const Vector& otherVector);
  Vector& operator=(const double constant);

  // Unary +
  Vector operator+() const;

  // Unary -
  Vector operator-() const;

  // Binary +
  Vector operator+(const Vector& v1) const;

  // Binary -
  Vector operator-(const Vector& v1) const;

  // Scalar multiplication
  Vector operator*(double a) const;

  // p-norm method
  double CalculateNorm(int p=2) const;

  // infinity-norm method
  double CalculateInfinityNorm() const;

  // Scalar product with another vector
  double ScalarProduct(const Vector& v) const;

  // Declare length function as a friend
  friend int length(const Vector& v);

  // Override << operator
  friend std::ostream& operator<<(std::ostream& output,
                       const Vector& v);

  // Find maximum value in vector // ADDED BY SAMMY PETROS
  double FindMax() const;

  // Equivalence // ADDED BY SAMMY PETROS
  bool operator==(const Vector& v1) const;

  // Non-Equivalence // ADDED BY SAMMY PETROS
  bool operator!=(const Vector& v1) const;

  // Bubble Sort // ADDED BY SAMMY PETROS
  void BubbleSort();

};

// Prototype signature of length() friend function
int length(const Vector& v);

#endif

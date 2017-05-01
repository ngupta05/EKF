#include "Matrix.h"
#include <math.h>
#include <stdexcept>

Matrix::Matrix(int n_rows, int n_cols): m_rows(n_rows), m_cols(n_cols) {
  m_matrix.resize(m_rows, std::vector<float>(m_cols, 0));
}

Matrix Matrix::operator+(const Matrix& rhs) const {
  if (m_rows != rhs.getNumRows() || m_cols != rhs.getNumCols()) {
   throw std::invalid_argument("Matrix dimensions not same");
  } 
  Matrix result(m_rows, m_cols);
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < m_cols; j++)
      result(i, j) = (*this).get(i, j) + rhs.get(i, j);
  return result;
}

Matrix Matrix::operator-(const Matrix& rhs) const {
  if (m_rows != rhs.getNumRows() || m_cols != rhs.getNumCols()) {
   throw std::invalid_argument("Matrix dimensions not same");
  } 
  Matrix result(m_rows, m_cols);
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < m_cols; j++)
      result(i, j) = (*this).get(i, j) - rhs.get(i, j);
  return result;
}

Matrix Matrix::operator*(const Matrix& rhs) const {
  if (m_cols != rhs.getNumRows()) {
   throw std::invalid_argument("Matrix dimensions not same");
  }

  Matrix result(m_rows, rhs.getNumCols());
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < rhs.getNumCols(); j++) {
      float res = 0;
      for (int k = 0; k < m_cols; k++)
        res += (*this).get(i, k) * rhs.get(k, j);
      result(i, j) = res;
    }
  return result;
}

Matrix Matrix::transpose() const {
  Matrix result(m_cols, m_rows);
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < m_cols; j++)
      result(j, i) = this->get(i, j);
  return result;
}

float& Matrix::operator()(int row, int col) {
  if (row >= m_rows || col >= m_cols)
    throw std::invalid_argument("Invalid index");
  return m_matrix[row][col];
}

#define N 10
// Function to get cofactor of A[p][q] in temp[][]. n is current
// dimension of A[][]
void getCofactor(float A[N][N], float temp[10][10], int p, int q, int n) {
  int i = 0, j = 0;
  // Looping for each element of the matrix
  for (int row = 0; row < n; row++) {
    for (int col = 0; col < n; col++) {
      //  Copying into temporary matrix only those element
      //  which are not in given row and column
      if (row != p && col != q) {
        temp[i][j++] = A[row][col];
 
        // Row is filled, so increase row index and
        // reset col index
        if (j == n - 1) {
          j = 0;
          i++;
        }
      }
    }
  }
}
 
/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
float determinant(float A[N][N], int n) {
  float D = 0; // Initialize result
 
  //  Base case : if matrix contains single element
  if (n == 1)
    return A[0][0];
 
  float temp[N][N]; // To store cofactors
 
  int sign = 1;  // To store sign multiplier
 
  // Iterate for each element of first row
  for (int f = 0; f < n; f++) {
    // Getting Cofactor of A[0][f]
    getCofactor(A, temp, 0, f, n);
    D += sign * A[0][f] * determinant(temp, n - 1);
 
    // terms are to be added with alternate sign
    sign = -sign;
  }
  return D;
}
 
// Function to get adjoint of A[N][N] in adj[N][N].
void adjoint(float A[N][N], float adj[N][N], int n) {
  if (n == 1) {
    adj[0][0] = 1;
    return;
  }
 
  // temp is used to store cofactors of A[][]
  int sign = 1;
  float temp[N][N];
 
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      // Get cofactor of A[i][j]
      getCofactor(A, temp, i, j, n);
 
      // sign of adj[j][i] positive if sum of row
      // and column indexes is even.
      sign = ((i+j)%2==0)? 1: -1;
 
      // Interchanging rows and columns to get the
      // transpose of the cofactor matrix
      adj[j][i] = (sign)*(determinant(temp, n-1));
    }
  }
}
 
// Function to calculate and store inverse, returns false if
// matrix is singular
Matrix Matrix::inverse() const {
  if (m_rows != m_cols) {
    throw std::invalid_argument("Cannot invert non-sq matrix");
  }

  float A[N][N];
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < m_cols; j++)
      A[i][j] = m_matrix[i][j];

  Matrix inverse(m_rows, m_cols);
  // Find determinant of A[][]
  float det = determinant(A, m_rows);
  if (fabs(det) < 1e-8) {
    throw std::invalid_argument("Singular matrix, can't find its inverse");
  }
 
  // Find adjoint
  float adj[N][N];
  adjoint(A, adj, m_rows);
 
  // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < m_cols; j++)
      inverse(i, j) = adj[i][j]/ det;
 
  return inverse;
}

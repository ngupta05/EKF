#pragma once

#include <vector>

class Matrix {
  int m_rows;
  int m_cols;

  std::vector<std::vector<float> > m_matrix;

  public:
  Matrix(int n_rows, int n_cols);
  Matrix operator+(const Matrix& rhs) const;
  Matrix operator-(const Matrix& rhs) const;
  Matrix operator*(const Matrix& rhs) const;
  Matrix inverse() const;
  Matrix inverse_cofactors() const; // minors, cofactor method
  Matrix transpose() const;
  float& operator()(int row, int col);
 
  float get(int row, int col) const { return m_matrix[row][col]; } 
  int getNumRows() const { return m_rows; }
  int getNumCols() const { return m_cols; }
};



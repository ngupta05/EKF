#include "Matrix.h"

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

Matrix Matrix::inverse() const {
  return *this;
}

Matrix Matrix::transpose() const {
  Matrix result(m_cols, m_rows);
  for (int i = 0; i < m_rows; i++)
    for (int j = 0; j < m_cols; j++)
      result(j, i) = this->get(i, j);
  return result;
}

float& Matrix::operator()(int row, int col) {
  return m_matrix[row][col];
}

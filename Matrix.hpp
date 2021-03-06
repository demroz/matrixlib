/* Matrix Operations library for C++
   Requires Vector.hpp for matrix/vector multiplication*/

#ifndef MATRIX_HPP
#define MATRIX_HPP

#include "Vector.hpp"

class Matrix
{
private:
  int rows_;
  int cols_;
  //  Vector data_;
  std::vector<std::vector<double> > m;

public:
  //Constructor
  Matrix(int rows, int cols);

  //get row/col
  int numRows() const ;
  int numCols() const ;
  //access indiidual matrix elements
  const double& operator()(const int& row, const int& col) const;
  double& operator()(const int& row,const int& col);
  //matrix matrix operators
  Matrix operator+(const Matrix& m_RHS);
  Matrix operator-(const Matrix& m_RHS);
  Matrix operator*(const Matrix& m_RHS);

  //matrix/scalar
  Matrix operator+(const double& a);
  Matrix operator-(const double& a); // get rid of
  Matrix operator*(const double& a);
  Matrix operator/(const double& a);
  //Matrix/vector operations
  Vector operator*(const Vector& v);
  // Gaussian elimination
  Vector Solve();  
};

inline std::ostream& operator<<(std::ostream& os, const Matrix& m_RHS)
{
  os << "[";
  for (int i=0; i<m_RHS.numRows(); i++)
    {
      for(int j=0; j<m_RHS.numCols(); j++)
	{
	  if(j!=0) os << ", ";
	  os << m_RHS(i,j);
	}
      if(i!=m_RHS.numRows()-1) os << ";\n";
    }
  os << "]";
  return os;	
  

}


inline Matrix I(int rows, int cols)
{
  if(rows!=cols)
    {
      throw std::domain_error("matrix must be square");
      
    }
  Matrix rtn(rows,cols);
  int j = 0;
  for(int i = 0; i<rows;i++)
    {
      rtn(i,j) = 1;
      j++;      
    }

  return rtn;
}


#endif

#include "Matrix.hpp"
#define TINY 1.0e-20

//Constructor - Matrix Initializer
Matrix::Matrix(int rows, int cols)
{
  std::vector<double> temp_vec(cols);
  m.resize(rows,temp_vec);
  rows_ = rows;
  cols_ = cols;

  temp_vec.clear();
}

// number of rows/colums
int Matrix::numRows() const
{
  return rows_;
}

int Matrix::numCols() const
{
  return cols_;
}

// acess individual matrix elements
const double& Matrix::operator()(const int& row, const int& col) const
{
  return m[row][col];
}
double& Matrix::operator()(const int& row,const int& col)
{
  return m[row][col];
}
//matrix matrix operations
//Addition
Matrix Matrix::operator+(const Matrix& m_RHS)
{
  int rows = m_RHS.numRows();
  int cols = m_RHS.numCols();
  
  Matrix rtn(rows, cols);

  for (int i =0; i<rows; i++)
    {
      for(int j=0; j<cols; j++)
	{
	  rtn(i,j) = m[i][j] + m_RHS(i,j);
	}
    }

  return rtn;
}
//Subtraction
Matrix Matrix::operator-(const Matrix& m_RHS)
{
  int rows = m_RHS.numRows();
  int cols = m_RHS.numCols();
  Matrix rtn(rows, cols);

  for (int i =0; i<rows; i++)
    {
      for(int j=0; j<cols; j++)
	{
	  rtn(i,j) = m[i][j] - m_RHS(i,j);
	}
    }

  return rtn;
}
//Multiplication (AB =/= BA)
Matrix Matrix::operator*(const Matrix& m_RHS)
{
  int rows_RHS = m_RHS.numRows();
  int cols_RHS = m_RHS.numCols();
  Matrix rtn(rows_RHS, cols_RHS);

  //fill rtn with zeros
  for(int i=0; i < rows_RHS; i++)
    {
      for(int j=0; j < cols_RHS;j++)
	{
	  rtn(i,j) = 0;
	}
    }
  
  // Matrix Multiplication O(n^3)
  for (int i=0; i<rows_RHS; i++) {
    for (int j=0; j<cols_RHS; j++) {
      for (int k=0; k<rows_RHS; k++) {
        rtn(i,j) += m[i][k] * m_RHS(k,j);
      }
    }
  }
  //  rows_RHS.clear();
  //  cols_RHS.clear();
  
  return rtn;

}

// Matrix inversion

// Matrix/Scalar operation
// Addition
Matrix Matrix::operator+(const double& a)
{
  Matrix rtn(rows_, cols_);

  for (int i = 0; i< rows_; i++)
    {
      for(int j = 0; j<cols_; j++)
	{
	  rtn(i,j) = m[i][j] + a;
	}
    }

  return rtn;  
}
//Subtration 
Matrix Matrix::operator-(const double& a)
{
  Matrix rtn(rows_, cols_);

  for (int i = 0; i< rows_; i++)
    {
      for(int j = 0; j<cols_; j++)
	{
	  rtn(i,j) = m[i][j] - a;
	}
    }

  return rtn;  
}
// Multiplication
Matrix Matrix::operator*(const double& a)
{
  Matrix rtn(rows_, cols_);

  for (int i = 0; i< rows_; i++)
    {
      for(int j = 0; j<cols_; j++)
	{
	  rtn(i,j) = m[i][j] * a;
	}
    }

  return rtn;  
}

//Division
Matrix Matrix::operator/(const double& a)
{
  Matrix rtn(rows_, cols_);

  for (int i = 0; i< rows_; i++)
    {
      for(int j = 0; j<cols_; j++)
	{
	  rtn(i,j) = m[i][j] / a;
	}
    }

  return rtn;  
}

// Matrix/Vector multiplication
// Column Vector
// test if algorithm is correct?
Vector Matrix::operator*(const Vector& v)
{
  Vector rtn(rows_);
  
  for(int i = 0; i<rows_; i++)
    {
      for(int j = 0; j < cols_; j++)
	{
	  rtn[i] += m[i][j] * v[j];
	}
    }
  return rtn;
}

//Gaussian Elimination

Vector Matrix::Solve()
{
  double maxEl,maxRow,tmp;
  Vector rtn(rows_);
  Matrix mCopy(rows_,cols_);

  //Copy matrix so as not to modify it later on
  for(int i=0;i<rows_;i++)
    {
      for(int j=0;j<cols_;j++)
	{
	  mCopy(i,j) = m[i][j];
	}
    }
  
  //search for maximum
  for(int i = 0; i<rows_; i++)
    {
      maxEl = std::abs(mCopy(i,i));
      maxRow = i;
      for (int k = i+1; k<rows_; k++)
	{
	  if(std::abs(mCopy(k,i)) > maxEl)
	    {
	      maxEl = std::abs(mCopy(k,i));
	      maxRow = k;
	    }
	}
      //Swap max row wit current row coum by colum
      for(int k=i; k<rows_+1; k++)
	{
	  tmp = mCopy(maxRow,k);
	  mCopy(maxRow,k) = mCopy(i,k);
	  mCopy(i,k) = tmp;
	}
      //make all rows below this one 0 in the current column
      for(int k=i+1; k<rows_; k++)
	{
	  double c = -mCopy(k,i)/mCopy(i,i);
	  for(int j=i; j<rows_+1;j++)
	    {
	      if(i==j)
		{
		  mCopy(k,j) =0;
		}
	      else
		{
		  mCopy(k,j) += c*mCopy(i,j);
		}
	    }
	}
    }
  for(int i = rows_-1; i>=0;i--)
    {
      rtn[i] = mCopy(i,rows_)/mCopy(i,i);
      for(int k=i-1;k>=0;k--)
	{
	  mCopy(k,rows_) -= mCopy(k,i)*rtn[i];
	}
    }

  return rtn;
}


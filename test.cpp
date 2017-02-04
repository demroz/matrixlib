#include<iostream>
#include<cmath>
#include<fstream>
#include<string>

#include "Matrix.hpp"

using namespace std;
int main()
{
  // Create Matrix A
  Matrix A(3,3);
  int k = 0;
  for (int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
	{
	  k++;
	  A(i,j) = k;
	}
    }
  cout << "Matrix A " << endl << A << endl;

  //Create Matrix B
  Matrix B(3,3);
  k=0;
  for (int i = 0; i<3; i++)
    {
      for(int j = 0; j<3; j++)
	{
	  k++;
	  B(i,j) = k;
	}
    }
  cout << "Matrix B " << endl <<  B << endl;
  // Matrx Matrix operations

  //Addition
  
  cout << "A + B" << endl<< A+B << endl;

  //Subtraction
  cout <<"A-B" << endl << A-B << endl;

  // Multiplication

  cout <<"A*B"<<endl<< A*B << endl;

  // Matrix/Scalar
  // Addition
  cout << "A+2" << endl <<A+2. << endl;
  //  Subtraction
  cout << "A-2" << endl << A-2. << endl;
  // multiplication
  cout << "A*2" << endl<< A*2. << endl;
  // division
  cout << "A/2" << endl << A/2.0 << endl;

  // Vector definiton
  Vector b(3);
  b[0] = 1.;
  b[1] = 2.;
  b[2] = 3.;
  cout << "Vector b " << endl << b << endl;

  //Matrix/vector multiplication
  cout << "A*b" << endl << A*b << endl;

  ///

  // Solve

  //Define Augmented Matrix

  Matrix U(3,4);
  k=0;
  for(int i =0;i<3;i++)
    {
      for(int j=0;j<4;j++)
	{
	  k++;
	  U(i,j) = k;
	}
    }
  // define solution vector
  Vector Sol(3);
  Sol = U.Solve();

  cout << "solution to matrix U " << endl << U << endl;
  cout <<"is"<< endl<<Sol <<endl;

  // Define Identity Matrix

  Matrix Identity(3,3);
  Identity=I(3,3);

  cout <<"3x3 identity matrix"<<endl<< Identity << endl;
}

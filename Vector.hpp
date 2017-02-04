#ifndef VECTOR_HPP
#define VECTOR_HPP

#include <vector>
#include <cmath>
#include <iostream>



class Vector
{
public:
  Vector(int n) : x_(n) {}

  int dim() const {return x_.size();}

  const double& operator[](int i) const {return x_[i];}

  double& operator[](int i) {return x_[i];}

  Vector operator*(const double& a) const ;

  double operator*(const Vector& v) const ;

  Vector operator+(const Vector& v) const ;

  Vector operator-(const Vector& v) const ;

  Vector operator-() const ;

  double norm2() const {return ::sqrt((*this)*(*this));}
private:
  std::vector<double> x_;
};

inline Vector operator*(const double& a, const Vector& v) {return v*a;}

inline std::ostream& operator<<(std::ostream& os, const Vector& v)
{
  os << "[";
  for (int i=0; i<v.dim(); i++)
    {
      if (i!=0) os << ", ";
      os << v[i];
    }
  os << "]";
  return os;
}

#endif

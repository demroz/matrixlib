#include "Vector.hpp"

Vector Vector::operator*(const double& a) const
{
  Vector rtn(dim());
  for (int i=0; i<dim(); i++) rtn[i] = x_[i]*a;

  return rtn;
}

double Vector::operator*(const Vector& v) const
{
  double rtn = 0.0;
  
  for (int i=0; i<dim(); i++) rtn += x_[i]*v.x_[i];

  return rtn;
}

Vector Vector::operator+(const Vector& v) const
{
  Vector rtn(dim());

  for (int i=0; i<dim(); i++) rtn[i] = x_[i] + v[i];

  return rtn;
}

Vector Vector::operator-(const Vector& v) const
{
  Vector rtn(dim());

  for (int i=0; i<dim(); i++) rtn[i] = x_[i] - v[i];

  return rtn;
}

Vector Vector::operator-() const
{
  Vector rtn(dim());

  for (int i=0; i<dim(); i++) rtn[i] = -x_[i];

  return rtn;
}


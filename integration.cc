#include<iostream>
#include<cmath>
#include<fstream>
#include<string>

#include "Matrix.hpp"

double f(double y, double t)
{
  return -5000*y+cos(t);
}

double exact_sol(double t)
{
  return (-5000*exp(-5000*t)+sin(t)+5000*cos(t))/(25000001);
}

double backward_euler(double y0, double tstart, double tend, double h)
{
  double y,yexact,err;
  y = y0;

  std::ofstream outfile;
  outfile.open("./backward_euler.dat");

  
  for(double t = tstart; t< tend; t=t+h)
    {
      //backward euler algorithm
      yexact = exact_sol(t);
      err = pow(pow((y-yexact),2),0.5);
      outfile << t << " " << y << " "<<yexact<<" " << err << std::endl;

      y = (y + h*cos(t))/(5000*h+1);

      // yexact
      //yexact = exact_sol(t+h);
      //err = pow(pow((y-yexact),2),0.5);
      //outfile << t << " " << y << " " << err << std::endl;
    }
  outfile.close();
  return err;
}


double GaussLegendre(double y0, double tstart, double tend, double h)
{
  double y, yexact, err;
  //Initialize A, bTranspose, k, and c matrices
  Matrix A(2,2),k(2,3);
  Vector bT(2),c(2);
  // Matrix A
  A(0,0) = 1./4.;
  A(0,1) = 1./4.-sqrt(3.)/6.;
  A(1,0) = 1./4.+sqrt(3.)/6.;
  A(1,1) = 1./4.;
  //matrix bT
  bT[0] = 0.5;
  bT[1] = 0.5;
  //Matrix c
  c[0] = 1./2.-sqrt(3.)/6.;
  c[1] = 1./2.+sqrt(3.)/6.;
  //Matrix k

  ///////////////////////////
  ///  initialize soln vec///
  ///////////////////////////
  Vector solK(2);

  std::ofstream outfile;
  outfile.open("./gauss_legendre.dat");
  k(0,0) = 1.+5000*h*A(0,0);
  k(0,1) = 5000.*h*A(0,1);
  k(1,0) = 5000.*h*A(1,0);
  k(1,1) = 1.+5000.*h*A(1,1);
  for(double t = tstart; t<tend; t=t+h)
    {

      k(0,2) = -5000.*y+cos(t+c[0]*h);
      k(1,2) = -5000.*y+cos(t+c[1]*h);
      //Solve system of equations
      solK = k.Gauss();
      yexact = exact_sol(t);
      err = pow(pow((y-yexact),2),0.5);
      outfile << t << " " << y << " " << yexact << " "<< err << std::endl;
      y = y+h*(bT[0]*solK[0]+bT[1]*solK[1]);

    }

  return err;
  
}

double Radau(double y0, double tstart, double tend, double h)
{
  double y, yexact, err;
  //Initialize A, bTranspose, k, and c matrices
  Matrix A(3,3),k(3,4);
  Vector bT(3),c(3),solK(3);
  // Matrix A
  A(0,0) = 11./45.-7.*sqrt(6.)/360.;
  A(0,1) = 37./225.-169.*sqrt(6.)/1800.;
  A(0,2) = -2./225.+sqrt(6.)/75.;
  
  A(1,0) = 37./225.+169.*sqrt(6.)/1800.;
  A(1,1) = 11./45.+7.*sqrt(6.)/360.;
  A(1,2) = -2./225.-sqrt(6.)/75.;
  
  A(2,0) = 4./9.-sqrt(6.)/36.;
  A(2,1) = 4./9.+sqrt(6.)/36.;
  A(2,2) = 1./9.;
  //matrix bT
  bT[0] = 4./9.-sqrt(6.)/36.;
  bT[1] = 4./9.+sqrt(6.)/36.;
  bT[2] = 1./9.;
  //Matrix c
  c[0] = 2./5.-sqrt(6.)/10.;
  c[1] = 2./5.+sqrt(6.)/10.;
  c[2] = 1.;
  // matrix k
  k(0,0) = 1.+5000.*h*A(0,0);
  k(0,1) = 5000.*h*A(0,1);
  k(0,2) = 5000.*h*A(0,2);
  
  k(1,0) = 5000.*h*A(1,0);
  k(1,1) = 1.+5000.*h*A(1,1);
  k(1,2) = 5000.*h*A(1,2);

  k(2,0) = 5000.*h*A(2,0);
  k(2,1) = 5000.*h*A(2,1);
  k(2,2) = 1.+5000.*h*A(2,2);

  /*  std::cout<<std::endl;

  for(int i=0;i<3;i++)
    {
      for(int j=0;j<3;j++)
	{
	  std::cout<<k(i,j)<<" ";
	}
      std::cout<<std::endl;
      }*/ //debugging purposes
  
  std::ofstream outfile;
  outfile.open("./radau.dat");

  for(double t = tstart; t<tend; t=t+h)
    {
      k(0,3) = -5000.*y+cos(t+c[0]*h);
      k(1,3) = -5000.*y+cos(t+c[1]*h);
      k(2,3) = -5000.*y+cos(t+c[2]*h);

      solK = k.Gauss();

      //std::cout << solK << std::endl;
      y = y+h*(bT[0]*solK[0]+bT[1]*solK[1]+bT[2]*solK[2]);

      
      yexact = exact_sol(t+h);
      err = pow(pow((y-yexact),2),0.5);
      outfile << t << " " << y << " " << yexact << " "<< err << std::endl;
    }
  return err;
}

double trapez(double y0, double tstart, double tend, double h)
{
  double y, yexact, err;
  //Initialize A, bTranspose, k, and c matrices
  Matrix A(2,2),k(2,3);
  Vector bT(2),c(2),solK(2);
  // Matrix A
  A(0,0) = 0.0;
  A(0,1) = 0.0;
  A(1,0) = 0.5;
  A(1,1) = 0.5;
  // Matrix bT
  bT[0] = 0.5;
  bT[1] = 0.5;
  //Matrix c
  c[0] = 0;
  c[1] = 1;
  //k matrix

  ///////////////////////////
  ///  initialize soln vec///
  ///////////////////////////
  
  std::ofstream outfile;
  outfile.open("./trapezoid.dat");
  k(0,0) = 1.+5000*h*A(0,0);
  k(0,1) = 5000.*h*A(0,1);
  k(1,0) = 5000.*h*A(1,0);
  k(1,1) = 1.+5000.*h*A(1,1);

  for(double t = tstart; t<tend; t=t+h)
    {
      
      k(0,2) = -5000.*y+cos(t+c[0]*h);
      k(1,2) = -5000.*y+cos(t+c[1]*h);
      //Solve system of equations
      solK = k.Gauss();
      yexact = exact_sol(t);
      err = pow(pow((y-yexact),2),0.5);

      outfile << t << " " << y << " " << yexact << " "<< err << std::endl;
      y = y+h*(bT[0]*solK[0]+bT[1]*solK[1]);

    }

  return err;
}

main()
{
  double y0, tstart,tend,h;
  y0 = 0;
  tstart = 0;
  tend = 1000;
  h = 10;

  std::ofstream outfile;
  outfile.open("./err_v_step.dat");
    for(int i = 1; i<3; i++)
      {      
	double euler_err, gauss_err,radau_err,trapez_err;
      
	euler_err = backward_euler(y0,tstart,tend,h);
	gauss_err = GaussLegendre(y0,tstart,tend,h);
	radau_err = Radau(y0,tstart,tend,h);
	trapez_err = trapez(y0,tstart,tend,h);
	
	outfile << h << " " << euler_err << " " << gauss_err << " " << radau_err << " " << trapez_err << std::endl;

	h = h/10;

	
      }
}

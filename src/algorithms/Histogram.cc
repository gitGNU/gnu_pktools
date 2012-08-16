/**********************************************************************
Histogram.cc: class for statistical operations on vectors
Copyright (C) 2008-2012 Pieter Kempeneers

This file is part of pktools

pktools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

pktools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with pktools.  If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/
#include "Histogram.h"


Histogram::Histogram(){
}

void Histogram::signature(double m1, double m2, double& k, double& alpha, double& beta, double e)
{
  double y=m1*m1/m2;
  beta=F_1(y,0.1,10.0,e);
  double fb=F(beta);
  double g=exp(lgamma(1.0/beta));
  alpha=m1*g/exp(lgamma(2.0/beta));
  k=beta/(2*alpha*g);
//   cout << "y, alpha, beta: " << y << ", " << alpha << ", " << beta << endl;
}

double Histogram::F(double x)
{
  double g2=exp(lgamma(2.0/x));
  return(g2*g2/exp(lgamma(3.0/x))/exp(lgamma(1.0/x)));
}

//x1 is under estimate, x2 is over estimate, e is error
double Histogram::F_1(double y, double x1, double x2, double e)
{
  double f1=F(x1);
  double f2=F(x2);
  assert(f1!=f2);
  double x=x1+(x2-x1)*(y-f1)/(f2-f1);
  double f=F(x);
  while(f-y>=e||y-f>=e){
    if(f<y)
      x1=x;
    else 
      x2=x;
    if(x1==x2)
      return x1;
    assert(f1!=f2);
    x=x1+(x2-x1)*(y-f1)/(f2-f1);
    f=F(x);
  }
  return x;
}

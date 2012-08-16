/**********************************************************************
Vector2d.h: 2-dimensional vector class (inherits from stl vector class)
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
#ifndef _VECTOR2D_H_
#define _VECTOR2D_H_

#include <vector>
#include <list>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_matrix.h>
#include "PosValue.h"
#include "algorithms/Histogram.h"

using namespace std;
template<class T> class Vector2d: public vector<vector <T> >
{
public:
  Vector2d();
  Vector2d(const Vector2d<T>& v1);//copy constructor
  ~Vector2d();
  Vector2d(int nrow);
  Vector2d(int nrow, int ncol);
  Vector2d(int nrow, int ncol, const T& value);
  Vector2d(const gsl_matrix* gsl_m);
  void resize(int nrow)
  {
    vector< vector<T> >::resize(nrow);
  };
  void resize(int nrow, int ncol);
  int nRows() const {return this->size();};
  int nCols() const {if(this->size()) return this->begin()->size(); else return 0;};
  void selectCol(int col, vector<T> &output) const;
  void selectCol(int col, T* output) const;
  vector<T> selectCol(int col);
  void selectCols(const list<int> &cols, Vector2d<T> &output) const;
  void selectCols(const list<int> &cols);
  void sort(Vector2d<T>& output);  
  void scale(const vector<double> &scaleVector, const vector<double> &offsetVector, Vector2d<T>& scaledOutput);
  void scale(const T lbound, const T ubound, vector<double> &scaleVector, vector<double> &offsetVector, Vector2d<T>& scaledOutput);
  Vector2d<T> operator=(const Vector2d<T>& v1);
//   ostream& operator<<(ostream& os, const Vector2d<T>& v);
//   template<class T> ostream& operator<<(ostream& os, const Vector2d<T>& v);
  Vector2d<T> sum(const Vector2d<T>& v1, const Vector2d<T>& v2) const;
  T max(int& x, int& y, double maxValue) const;

  T sum() const;
};
  
template<class T> Vector2d<T>::Vector2d() 
  : vector< vector<T> >()
{
}

template<class T> Vector2d<T>::~Vector2d() 
{
}

//copy constructor
template<class T> Vector2d<T>::Vector2d(const Vector2d<T>& v1){
  this->resize(v1.size());
  for(int irow=0;irow<v1.size();++irow)
    this->at(irow)=v1[irow];
}

template<class T> Vector2d<T> Vector2d<T>::operator=(const Vector2d<T>& v1){
  //check for assignment to self (of the form v=v)
  if(this==&v1)
     return *this;
  else{
    this->resize(v1.size());
    for(int irow=0;irow<v1.size();++irow)
      this->at(irow)=v1[irow];
    return *this;
  }
}

template<class T> Vector2d<T>::Vector2d(int nrow) 
  : vector< vector<T> >(nrow)
{
}

template<class T> Vector2d<T>::Vector2d(int nrow, int ncol) 
//   : vector< vector<T> >(nrow)
{
  this->resize(nrow);
  for(int irow=0;irow<nrow;++irow){
    (this->operator[](irow)).resize(ncol);
//     (*this)[irow].resize(ncol);
  }
}

template<class T> Vector2d<T>::Vector2d(int nrow, int ncol, const T& value) 
{
  this->resize(nrow);
  for(int irow=0;irow<nrow;++irow){
    (this->operator[](irow)).resize(ncol);
    for(int icol=0;icol<ncol;++icol)
      (this->operator[](irow))[icol]=value;
  }
}

template<class T> Vector2d<T>::Vector2d(const gsl_matrix* gsl_m)
{
  this->resize(gsl_m->size1);
  for(int irow=0;irow<this->size();++irow){
    (this->operator[](irow)).resize(gsl_m->size2);
    for(int icol=0;icol<this->operator[](irow).size();++icol)
      (this->operator[](irow))[icol]=gsl_matrix_get(gsl_m,irow,icol);
  }
}


template<class T> void Vector2d<T>::resize(int nrow, int ncol)
{
    this->vector< vector<T> >::resize(nrow);
    for(int irow=0;irow<nrow;++irow){
      (this->operator[](irow)).resize(ncol);
    }
}

template<class T> void Vector2d<T>::selectCols(const list<int> &cols, Vector2d<T> &output) const
{
  output.resize(this->size());
  list<int>::const_iterator it;
  for(int irow=0;irow<this->size();++irow){
    output[irow].resize(cols.size());
    it=cols.begin();
    for(int icol=0;icol<cols.size();++icol)
      output[irow][icol]=(*this)[irow][*(it++)];
  }
}

template<class T> void Vector2d<T>::selectCol(int col, vector<T> &output) const
{
  assert(col>=0);
  assert(col<(*this)[0].size());
  output.resize(this->size());
  for(int irow=0;irow<this->size();++irow){
    output[irow]=(*this)[irow][col];
  }
}

template<class T> vector<T> Vector2d<T>::selectCol(int col)
{
  assert(col>=0);
  assert(col<(*this)[0].size());
  vector<T> output(this->size());
  for(int irow=0;irow<this->size();++irow)
    output[irow]=(*this)[irow][col];
  return(output);
}

template<class T> void Vector2d<T>::selectCol(int col, T* output) const
{
  assert(col>=0);
  assert(col<(*this)[0].size());
  for(int irow=0;irow<this->size();++irow){
    output[irow]=(*this)[irow][col];
  }
}

template<class T> void Vector2d<T>::selectCols(const list<int> &cols)
{
  for(int irow=0;irow<this->size();++irow)
    for(int icol=((*this)[irow]).size()-1;icol>=0;--icol)
      if(find(cols.begin(),cols.end(),icol)==cols.end())
	(*this)[irow].erase(((*this)[irow]).begin()+icol);
}

// template<class T> ostream& operator<<(ostream& os, const Vector2d<T>& v)
// {
//   for(int irow=0;irow<v.size();++irow){
//     for(int icol=0;icol<v[irow].size();++icol){
//       os << v[irow][icol] << "\t";
//     }
//     os << endl;
//   }
//   return os;
// }

template<class T> void Vector2d<T>::sort(Vector2d<T>& output)
{
  //sort according to first sample (ex. wavelength)
  int nsample=this->size();//including first sample (ex. wavelength)
  int nband=(*this)[0].size();  
  vector<PosValue> sortW(nband);
  for(int ilevel=0;ilevel<nband;++ilevel){
    PosValue pv;
    pv.position=ilevel;
    pv.value=(*this)[0][ilevel];
    sortW[ilevel]=pv;
  }
  std::sort(sortW.begin(),sortW.end(),Increase_PosValue());
  output.resize(nsample);  
  for(int isample=0;isample<nsample;++isample){
    output[isample].resize(nband);
    for(int iband=0;iband<nband;++iband)
      output[isample][iband]=(*this)[isample][sortW[iband].position];
  }
}

template<class T> void Vector2d<T>::scale(const vector<double> &scaleVector,const vector<double> &offsetVector, Vector2d<T>& scaledOutput)
{
  int nsample=this->size();//including first sample (ex. wavelength)
  int nband=(*this)[0].size();
  assert(scaleVector.size()==nband);
  assert(offsetVector.size()==nband);
  vector<T> pixel(nband);
  scaledOutput.resize(nsample,nband);
  for(int isample=0;isample<nsample;++isample)
    for(int iband=0;iband<nband;++iband)
      scaledOutput[isample][iband]=((*this)[isample][iband])*scaleVector[iband]+offsetVector[iband];
}

template<class T> void Vector2d<T>::scale(const T lbound, const T ubound, vector<double> &scaleVector, vector<double> &offsetVector, Vector2d<T>& scaledOutput)
{
  //scale to lbound and ubound
  int nsample=this->size();//including first sample (ex. wavelength)
  int nband=(*this)[0].size();
  scaleVector.resize(nband);
  offsetVector.resize(nband);
  vector<T> pixel(nsample);
  T theMin;
  T theMax;
  Histogram hist;
  scaledOutput.resize(nsample,nband);
  for(int iband=0;iband<nband;++iband){
    pixel=selectCol(iband);
    hist.minmax(pixel, pixel.begin(), pixel.end(), theMin, theMax);
    scaleVector[iband]=static_cast<double>(ubound-lbound)/(theMax-theMin);
    offsetVector[iband]=static_cast<double>(-theMin*scaleVector[iband])-lbound;
    for(int isample=0;isample<pixel.size();++isample)
      scaledOutput[isample][iband]=((*this)[isample][iband])*scaleVector[iband]+offsetVector[iband];
  }
}

template<class T> Vector2d<T> Vector2d<T>::sum(const Vector2d<T>& v1, const Vector2d<T>& v2) const{
  Vector2d<T> vsum(v1.size());
  assert(v1.size()==v2.size());
  for(int irow=0;irow<v1.size();++irow){
    assert(v1[irow].size()==v2[irow].size());
    vsum[irow].resize(v1[irow].size());
    for(int icol=0;icol<v1.size();++icol)
      vsum[irow][icol]=v1[irow][icol]+v2[irow][icol];
  }
  return vsum;
}

template<class T> T Vector2d<T>::sum() const{
  double theSum=0;
  for(int irow=0;irow<this->size();++irow){
    for(int icol=0;icol<this->operator[](irow).size();++icol)
      theSum+=(this->operator[](irow))[icol];
  }
  return theSum;
}

template<class T> T Vector2d<T>::max(int& x, int& y, double maxValue) const{
  //todo: what if this->operator[](0)[0] >=maxValue?
  // double theMax=(this->operator[](0))[0];
  double theMax=0;
  for(int irow=0;irow<this->size();++irow){
    for(int icol=0;icol<(this->operator[](irow)).size();++icol){
      double currentValue=(this->operator[](irow))[icol];
      if(currentValue<maxValue&&currentValue>theMax){
        assert(theMax<maxValue);
        y=irow;
        x=icol;
        theMax=currentValue;
      }
    }
  }
  assert(theMax<maxValue);
  return theMax;
}

#endif /* _VECTOR2D_H_ */

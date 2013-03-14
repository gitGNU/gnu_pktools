/**********************************************************************
StatFactory.h: class for statistical operations on vectors
Copyright (C) 2008-2013 Pieter Kempeneers

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
#ifndef _STATFACTORY_H_
#define _STATFACTORY_H_

#include <iostream>
#include <vector>
#include <map>
#include <math.h>
#include <assert.h>
#include <string>
#include <sstream>
#include <fstream>
#include <algorithm>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

using namespace std;
// using namespace NEWMAT;

namespace statfactory
{

class StatFactory{

public:
  enum INTERPOLATION_TYPE {undefined=0,linear=1,polynomial=2,cspline=3,cspline_periodic=4,akima=5,akima_periodic=6};
  StatFactory(void){};
  virtual ~StatFactory(void){};
  INTERPOLATION_TYPE getInterpolationType(const std::string interpolationType){
    std::map<std::string, INTERPOLATION_TYPE> m_interpMap;
    initMap(m_interpMap);
    // gsl_interp_accel *acc=gsl_interp_accel_alloc();
    // gsl_spline *spline=gsl_spline_alloc(m_interpMap["interpolationType"],theSize);
    return m_interpMap[interpolationType];
  };
  static void allocAcc(gsl_interp_accel *&acc){
    acc = gsl_interp_accel_alloc ();
  };

  static void getSpline(const std::string type, int size, gsl_spline *& spline){
    std::map<std::string, INTERPOLATION_TYPE> m_interpMap;
    initMap(m_interpMap);
    switch(m_interpMap[type]){
    case(polynomial):
      spline=gsl_spline_alloc(gsl_interp_polynomial,size);
      break;
    case(cspline):
      spline=gsl_spline_alloc(gsl_interp_cspline,size);
      break;
    case(cspline_periodic):
      spline=gsl_spline_alloc(gsl_interp_cspline_periodic,size);
      break;
    case(akima):
      spline=gsl_spline_alloc(gsl_interp_akima,size);
      break;
    case(akima_periodic):
      spline=gsl_spline_alloc(gsl_interp_akima_periodic,size);
      break;
    case(linear):
    default:
      spline=gsl_spline_alloc(gsl_interp_linear,size);
    break;
    }
    assert(spline);
  };
  static void initSpline(gsl_spline *spline, const double *x, const double *y, int size){
    gsl_spline_init (spline, x, y, size);
  };
  static double evalSpline(gsl_spline *spline, double x, gsl_interp_accel *acc){
    return gsl_spline_eval (spline, x, acc);
  };


  template<class T> T max(const vector<T>& v) const;
  template<class T> T min(const vector<T>& v) const;
//   template<class T> typename vector<T>::const_iterator max(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const;
  template<class T> typename vector<T>::const_iterator max(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const;
  template<class T> typename vector<T>::iterator max(const vector<T>& v, typename vector<T>::iterator begin, typename vector<T>::iterator end) const;
  template<class T> typename vector<T>::const_iterator absmax(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const;
  template<class T> typename vector<T>::const_iterator min(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const;
  template<class T> typename vector<T>::iterator min(const vector<T>& v, typename vector<T>::iterator begin, typename vector<T>::iterator end) const;
  template<class T> typename vector<T>::const_iterator absmin(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const;
  template<class T> void minmax(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, T& theMin, T& theMax) const;  
  template<class T> T sum(const vector<T>& v) const;
  template<class T> double mean(const vector<T>& v) const;
  template<class T> T median(const vector<T>& v) const;
  template<class T> double var(const vector<T>& v) const;
  template<class T> double moment(const vector<T>& v, int n) const;
  template<class T> double cmoment(const vector<T>& v, int n) const;
  template<class T> double skewness(const vector<T>& v) const;
  template<class T> double kurtosis(const vector<T>& v) const;
  template<class T> void meanVar(const vector<T>& v, double& m1, double& v1) const;
  template<class T1, class T2> void  scale2byte(const vector<T1>& input, vector<T2>& output, unsigned char lbound=0, unsigned char ubound=255) const;
  template<class T> void distribution(const vector<T>& input, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end,  vector<int>& output, int nbin, T &minimum=0.0, T &maximum=0.0, const string &filename="");
  template<class T> void cumulative (const vector<T>& input, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, vector<int>& output, int nbin, T &minimum, T &maximum);
  template<class T> void  percentiles (const vector<T>& input, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, vector<T>& output, int nbin=10, T &minimum=0.0, T &maximum=0.0, const string &filename="");
  template<class T> void signature(const vector<T>& input, double& k, double& alpha, double& beta, double e);
  void signature(double m1, double m2, double& k, double& alpha, double& beta, double e);
  template<class T> void normalize(const vector<T>& input, vector<double>& output);
  template<class T> void normalize_pct(vector<T>& input);
  template<class T> double rmse(const vector<T>& x, const vector<T>& y) const;
  template<class T> double correlation(const vector<T>& x, const vector<T>& y, int delay=0) const;
  template<class T> double cross_correlation(const vector<T>& x, const vector<T>& y, int maxdelay, vector<T>& z) const;
  template<class T> double linear_regression(const vector<T>& x, const vector<T>& y, double &c0, double &c1) const;
  template<class T> void interpolateUp(const vector<double>& wavelengthIn, const vector<T>& input, const vector<double>& wavelengthOut, const std::string& type, vector<T>& output, bool verbose=false);
  template<class T> void interpolateUp(const vector<double>& wavelengthIn, const vector< vector<T> >& input, const vector<double>& wavelengthOut, const std::string& type, vector< vector<T> >& output, bool verbose=false);
  // template<class T> void interpolateUp(const vector< vector<T> >& input, vector< vector<T> >& output, double start, double end, double step, const gsl_interp_type* type);
  // template<class T> void interpolateUp(const vector< vector<T> >& input, const vector<double>& wavelengthIn, vector< vector<T> >& output, vector<double>& wavelengthOut, double start, double end, double step, const gsl_interp_type* type);
  template<class T> void interpolateUp(const vector<T>& input, vector<T>& output, int nbin);
  template<class T> void interpolateUp(double* input, int dim, vector<T>& output, int nbin);
  template<class T> void interpolateDown(const vector<T>& input, vector<T>& output, int nbin);
  template<class T> void interpolateDown(double* input, int dim, vector<T>& output, int nbin);

private:
  static void initMap(std::map<std::string, INTERPOLATION_TYPE>& m_interpMap){
    //initialize selMap
    m_interpMap["linear"]=linear;
    m_interpMap["polynomial"]=polynomial;
    m_interpMap["cspline"]=cspline;
    m_interpMap["cspline_periodic"]=cspline_periodic;
    m_interpMap["akima"]=akima;
    m_interpMap["akima_periodic"]=akima_periodic;
  }

};


template<class T> typename vector<T>::iterator StatFactory::max(const vector<T>& v, typename vector<T>::iterator begin, typename vector<T>::iterator end) const
{
  typename vector<T>::iterator tmpIt=begin;
  for (typename vector<T>::iterator it = begin; it!=end; ++it){
    if(*tmpIt<*it)
      tmpIt=it;
  }
  return tmpIt;
}

template<class T> T StatFactory::max(const vector<T>& v) const
{
  T maxValue=*(v.begin());
  for (typename vector<T>::const_iterator it = v.begin(); it!=v.end(); ++it){
    if(maxValue<*it)
      maxValue=*it;
  }
  return maxValue;
}

template<class T> T StatFactory::min(const vector<T>& v) const
{
  T minValue=*(v.begin());
  for (typename vector<T>::const_iterator it = v.begin(); it!=v.end(); ++it){
    if(minValue>*it)
      minValue=*it;
  }
  return minValue;
}

template<class T> typename vector<T>::const_iterator StatFactory::absmax(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const
{
  typename vector<T>::const_iterator tmpIt=begin;
  for (typename vector<T>::const_iterator it = begin; it!=end; ++it){
    if(abs(*tmpIt)<abs(*it))
      tmpIt=it;
  }
  return tmpIt;
}

template<class T> typename vector<T>::const_iterator StatFactory::min(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const
{
  typename vector<T>::const_iterator tmpIt=begin;
  for (typename vector<T>::const_iterator it = begin; it!=end; ++it){
    if(*tmpIt>*it)
      tmpIt=it;
  }
  return tmpIt;
}

template<class T> typename vector<T>::const_iterator StatFactory::absmin(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end) const
{
  typename vector<T>::const_iterator tmpIt=begin;
  for (typename vector<T>::const_iterator it = begin; it!=end; ++it){
    if(abs(*tmpIt)>abs(*it))
      tmpIt=it;
  }
}

template<class T> void StatFactory::minmax(const vector<T>& v, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, T& theMin, T& theMax) const
{
  theMin=*begin;
  theMax=*begin;
  for (typename vector<T>::const_iterator it = begin; it!=end; ++it){
    if(theMin>*it)
      theMin=*it;
    if(theMax<*it)
      theMax=*it;
  }
}

template<class T> T StatFactory::sum(const vector<T>& v) const
{
  typename vector<T>::const_iterator it;
  T tmpSum=0;
  for (it = v.begin(); it!= v.end(); ++it)
    tmpSum+=*it;
  return tmpSum;
}

template<class T> double StatFactory::mean(const vector<T>& v) const
{
  assert(v.size());
  return static_cast<double>(sum(v))/v.size();
}

template<class T> T StatFactory::median(const vector<T>& v) const
{
  vector<T> tmpV=v;
  sort(tmpV.begin(),tmpV.end());
  if(tmpV.size()%2)
    return tmpV[tmpV.size()/2];
  else
    return 0.5*(tmpV[tmpV.size()/2-1]+tmpV[tmpV.size()/2]);
}

template<class T> double StatFactory::var(const vector<T>& v) const
{
  typename vector<T>::const_iterator it;
  double v1=0;
  double m1=mean(v);
  double n=v.size();
  assert(n>1);
  for (it = v.begin(); it!= v.end(); ++it)
    v1+=(*it-m1)*(*it-m1);
  v1/=(n-1);
//   if(v1<0){
//     for (it = v.begin(); it!= v.end(); ++it)
//       cout << *it << " ";
//     cout << endl;
//   }
  assert(v1>=0);
  return v1;
}

template<class T> double StatFactory::moment(const vector<T>& v, int n) const
{
  assert(v.size());
  typename vector<T>::const_iterator it;
  double m=0;
//   double m1=mean(v);
  for(it = v.begin(); it!= v.end(); ++it){
    m+=pow((*it),n);
  }
  return m/v.size();
}

  //central moment
template<class T> double StatFactory::cmoment(const vector<T>& v, int n) const
{
  assert(v.size());
  typename vector<T>::const_iterator it;
  double m=0;
  double m1=mean(v);
  for(it = v.begin(); it!= v.end(); ++it){
    m+=pow((*it-m1),n);
  }
  return m/v.size();
}

template<class T> double StatFactory::skewness(const vector<T>& v) const
{
  return cmoment(v,3)/pow(var(v),1.5);
}

template<class T> double StatFactory::kurtosis(const vector<T>& v) const
{
  double m2=cmoment(v,2);
  double m4=cmoment(v,4);
  return m4/m2/m2-3.0;
}

template<class T> void StatFactory::meanVar(const vector<T>& v, double& m1, double& v1) const
{
  typename vector<T>::const_iterator it;
  v1=0;
  m1=mean(v);
  double n=v.size();
  assert(n>1);
  for (it = v.begin(); it!= v.end(); ++it)
    v1+=(*(it)-m1)*(*(it)-m1);
  v1/=(n-1);
  assert(v1>=0);
}

template<class T1, class T2> void StatFactory::scale2byte(const vector<T1>& input, vector<T2>& output, unsigned char lbound,  unsigned char ubound) const
{
  output.resize(input.size());
  T1 minimum=min(input);
  T1 maximum=max(input);
  assert(maximum>minimum);
  double scale=(ubound-lbound)/(maximum-minimum);
  for (int i=0;i<input.size();++i)
    output[i]=scale*(input[i]-(minimum))+lbound;
}

template<class T> void  StatFactory::distribution (const vector<T>& input, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, vector<int>& output, int nbin, T &minimum, T &maximum, const string &filename)
{
  if(maximum<=minimum)
    minmax(input,begin,end,minimum,maximum);
  // if(!minimum)
  //   minimum=*(min(input,begin,end));
  // if(!maximum)
  //   maximum=*(max(input,begin,end));
  if(maximum<=minimum){
    ostringstream s;
    s<<"Error: could not calculate distribution (min>=max) in " << filename;
    throw(s.str());
  }
  assert(nbin>1);
  assert(input.size());
  if(output.size()!=nbin){
    output.resize(nbin);
    for(int i=0;i<nbin;output[i++]=0);
  }
  typename vector<T>::const_iterator it;
  for(it=begin;it!=end;++it){
    if(*it==maximum)
      ++output[nbin-1];
    else if(*it>=minimum && *it<maximum)
      ++output[static_cast<int>(static_cast<double>((*it)-minimum)/(maximum-minimum)*nbin)];
  }
  if(!filename.empty()){
    ofstream outputfile;
    outputfile.open(filename.c_str());
    if(!outputfile){
      ostringstream s;
      s<<"Error opening distribution file , " << filename;
      throw(s.str());
    }
    for(int bin=0;bin<nbin;++bin)
      outputfile << (maximum-minimum)*bin/(nbin-1)+minimum << " " << static_cast<double>(output[bin])/input.size() << endl;
    outputfile.close();
  }
}

template<class T> void  StatFactory::percentiles (const vector<T>& input, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, vector<T>& output, int nbin, T &minimum, T &maximum, const string &filename)
{
  if(maximum<=minimum)
    minmax(input,begin,end,minimum,maximum);
  // if(!minimum)
  //   minimum=*(min(input,begin,end));
  // if(!maximum)
  //   maximum=*(max(input,begin,end));
  assert(maximum>minimum);
  assert(nbin>1);
  assert(input.size());
  output.resize(nbin);
  std::vector<T> inputSort;
  inputSort.assign(begin,end);
  typename std::vector<T>::iterator vit=inputSort.begin();
  while(vit!=inputSort.end()){
    if(*vit<minimum||*vit>maximum)
      inputSort.erase(vit);
    else
      ++vit;
  }
  std::sort(inputSort.begin(),inputSort.end());
  vit=inputSort.begin();
  std::vector<T> inputBin;
  for(int ibin=0;ibin<nbin;++ibin){
    inputBin.clear();
    while(inputBin.size()<inputSort.size()/nbin&&vit!=inputSort.end()){
      inputBin.push_back(*vit);
      ++vit;
    }
    if(inputBin.size())
      output[ibin]=median(inputBin);
  }
  if(!filename.empty()){
    ofstream outputfile;
    outputfile.open(filename.c_str());
    if(!outputfile){
      ostringstream s;
      s<<"error opening distribution file , " << filename;
      throw(s.str());
    }
    for(int ibin=0;ibin<nbin;++ibin)
      outputfile << ibin*100.0/nbin << " " << static_cast<double>(output[ibin])/input.size() << endl;
    outputfile.close();
  }
}

// template<class T> void  StatFactory::cumulative (const vector<T>& input, typename vector<T>::const_iterator begin, typename vector<T>::const_iterator end, vector<int>& output, int nbin, T &minimum, T &maximum)
// {
//   assert(nbin>1);
//   assert(input.size());
//   distribution(input,output,nbin,minimum,maximum);
//   for(vector<int>::iterator it=begin+1;it!=end;++it)
//     *it+=*(it-1);
//   if(!filename.empty()){
//     ofstream outputfile;
//     outputfile.open(filename.c_str());
//     if(!outputfile){
//       ostringstream s;
//       s<<"error opening cumulative file , " << filename;
//       throw(s.str());
//     }
//     for(int bin=0;bin<nbin;++bin)
//       outputfile << (maximum-minimum)*bin/(nbin-1)+minimum << " " << static_cast<double>(output[bin])/input.size() << endl;
//     outputfile.close();
//   }
// }

template<class T> void StatFactory::signature(const vector<T>& input, double&k, double& alpha, double& beta, double e)
{
  double m1=moment(input,1);
  double m2=moment(input,2);
  signature(m1,m2,k,alpha,beta,e);
}

template<class T> void StatFactory::normalize(const vector<T>& input, vector<double>& output){
  double total=sum(input);
  if(total){
    output.resize(input.size());
    for(int index=0;index<input.size();++index)
      output[index]=input[index]/total;
  }
  else
    output=input;
}

template<class T> void StatFactory::normalize_pct(vector<T>& input){
  double total=sum(input);
  if(total){
    typename vector<T>::iterator it;
    for(it=input.begin();it!=input.end();++it)
      *it=100.0*(*it)/total;
  }
}
 
template<class T> double StatFactory::rmse(const vector<T>& x, const vector<T>& y) const{
  assert(x.size()==y.size());
  assert(x.size());
  double mse=0;
  for(int isample=0;isample<x.size();++isample){
    double e=x[isample]-y[isample];
    mse+=e*e/x.size();
  }
  return sqrt(mse);
}

template<class T> double StatFactory::correlation(const vector<T>& x, const vector<T>& y, int delay) const{
  double meanX=0;
  double meanY=0;
  double varX=0;
  double varY=0;
  double sXY=0;
  meanVar(x,meanX,varX);
  meanVar(y,meanY,varY);
  double denom = sqrt(varX*varY);
  if(denom){
    //Calculate the correlation series
    sXY = 0;
    for (int i=0;i<x.size();++i) {
      int j = i + delay;
      if (j < 0 || j >= y.size())
        continue;
      else{
        assert(i>=0&&i<x.size());
        assert(j>=0&&j<y.size());
        sXY += (x[i] - meanX) * (y[j] - meanY);
      }
    }
    double minSize=(x.size()<y.size())?x.size():y.size();
    return(sXY / denom / (minSize-1));
  }
  else
    return 0;
}

template<class T> double StatFactory::cross_correlation(const vector<T>& x, const vector<T>& y, int maxdelay, vector<T>& z) const{
  z.clear();
  double sumCorrelation=0;
  for (int delay=-maxdelay;delay<maxdelay;delay++) {
    z.push_back(correlation(x,y,delay));
    sumCorrelation+=z.back();
  }
  return sumCorrelation;
}

  template<class T> double StatFactory::linear_regression(const vector<T>& x, const vector<T>& y, double &c0, double &c1) const{
  assert(x.size()==y.size());
  assert(x.size());
  double cov00;
  double cov01;
  double  cov11;
  double sumsq;
  gsl_fit_linear(&(x[0]),1,&(y[0]),1,x.size(),&c0,&c1,&cov00,&cov01,&cov11,&sumsq);
  return (1-sumsq/var(y)/(y.size()-1));
}

//alternatively: use GNU scientific library:
// gsl_stats_correlation (const double data1[], const size_t stride1, const double data2[], const size_t stride2, const size_t n)

template<class T> void StatFactory::interpolateUp(const vector<double>& wavelengthIn, const vector<T>& input, const vector<double>& wavelengthOut, const std::string& type, vector<T>& output, bool verbose){
  assert(wavelengthIn.size());
  assert(input.size()==wavelengthIn.size());
  assert(wavelengthOut.size());
  int nband=wavelengthIn.size();
  output.clear();
  gsl_interp_accel *acc;
  allocAcc(acc);
  gsl_spline *spline;
  getSpline(type,nband,spline);
  assert(spline);
  assert(&(wavelengthIn[0]));
  assert(&(input[0]));
  initSpline(spline,&(wavelengthIn[0]),&(input[0]),nband);
  for(int index=0;index<wavelengthOut.size();++index){
    if(type=="linear"){
      if(wavelengthOut[index]<wavelengthIn.back()){
        output.push_back(*(input.begin()));
        continue;
      }
      else if(wavelengthOut[index]>wavelengthIn.back()){
        output.push_back(input.back());
        continue;
      }
    }
    double dout=evalSpline(spline,wavelengthOut[index],acc);
    output.push_back(dout);
  }
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
}

// template<class T> void StatFactory::interpolateUp(const vector<double>& wavelengthIn, const vector< vector<T> >& input, const vector<double>& wavelengthOut, const std::string& type, vector< vector<T> >& output, bool verbose){
//   assert(wavelengthIn.size());
//   assert(wavelengthOut.size());
//   int nsample=input.size();  
//   int nband=wavelengthIn.size();
//   output.clear();
//   output.resize(nsample);
//   gsl_interp_accel *acc;
//   allocAcc(acc);
//   gsl_spline *spline;
//   getSpline(type,nband,spline);
//   for(int isample=0;isample<nsample;++isample){
//     assert(input[isample].size()==wavelengthIn.size());
//     initSpline(spline,&(wavelengthIn[0]),&(input[isample][0]),nband);      
//     for(int index=0;index<wavelengthOut.size();++index){
//       if(type=="linear"){
//         if(wavelengthOut[index]<wavelengthIn.back())
//           output[isample].push_back(*(input.begin()));
//         else if(wavelengthOut[index]>wavelengthIn.back())
//           output[isample].push_back(input.back());
//       }
//       else{
//         double dout=evalSpline(spline,wavelengthOut[index],acc);
//         output[isample].push_back(dout);
//       }
//     }
//   }
//   gsl_spline_free(spline);
//   gsl_interp_accel_free(acc);
// }

template<class T> void StatFactory::interpolateUp(const vector<T>& input, vector<T>& output, int nbin)
{
  assert(input.size());
  assert(nbin);
  output.clear();
  int dim=input.size();
  for(int i=0;i<dim;++i){
    double deltaX=0;
    double left=input[i];
    if(i<dim-1){
      double right=(i<dim-1)? input[i+1]:input[i];
      deltaX=(right-left)/static_cast<double>(nbin);
      for(int x=0;x<nbin;++x){
        output.push_back(left+x*deltaX);
      }
    }
    else
      output.push_back(input.back());
  }
}

template<class T> void StatFactory::interpolateUp(double* input, int dim, vector<T>& output, int nbin)
{
  assert(nbin);
  output.clear();
  for(int i=0;i<dim;++i){
    double deltaX=0;
    double left=input[i];
    if(i<dim-1){
      double right=(i<dim-1)? input[i+1]:input[i];
      deltaX=(right-left)/static_cast<double>(nbin);
      for(int x=0;x<nbin;++x){
        output.push_back(left+x*deltaX);
      }
    }
    else
      output.push_back(input[dim-1]);
  }
}

template<class T> void StatFactory::interpolateDown(const vector<T>& input, vector<T>& output, int nbin)
{
  assert(input.size());
  assert(nbin);
  output.clear();
  int dim=input.size();
  int x=0;
  output.push_back(input[0]);
  for(int i=1;i<dim;++i){
    if(i%nbin)
      continue;
    else{
      x=(i-1)/nbin+1;
      output.push_back(input[i]);
    }
  }
}

template<class T> void StatFactory::interpolateDown(double* input, int dim, vector<T>& output, int nbin)
{
  assert(nbin);
  output.clear();
  int x=0;
  output.push_back(input[0]);
  for(int i=1;i<dim;++i){
    if(i%nbin)
      continue;
    else{
      x=(i-1)/nbin+1;
      output.push_back(input[i]);
    }
  }
}
}

#endif /* _STATFACTORY_H_ */

// void Histogram::signature(double m1, double m2, double& k, double& alpha, double& beta, double e)
// {
//   double y=m1*m1/m2;
//   beta=F_1(y,0.1,10.0,e);
//   double fb=F(beta);
//   double g=exp(lgamma(1.0/beta));
//   alpha=m1*g/exp(lgamma(2.0/beta));
//   k=beta/(2*alpha*g);
// //   cout << "y, alpha, beta: " << y << ", " << alpha << ", " << beta << endl;
// }

// double Histogram::F(double x)
// {
//   double g2=exp(lgamma(2.0/x));
//   return(g2*g2/exp(lgamma(3.0/x))/exp(lgamma(1.0/x)));
// }

// //x1 is under estimate, x2 is over estimate, e is error
// double Histogram::F_1(double y, double x1, double x2, double e)
// {
//   double f1=F(x1);
//   double f2=F(x2);
//   assert(f1!=f2);
//   double x=x1+(x2-x1)*(y-f1)/(f2-f1);
//   double f=F(x);
//   while(f-y>=e||y-f>=e){
//     if(f<y)
//       x1=x;
//     else 
//       x2=x;
//     if(x1==x2)
//       return x1;
//     assert(f1!=f2);
//     x=x1+(x2-x1)*(y-f1)/(f2-f1);
//     f=F(x);
//   }
//   return x;
// }

/**********************************************************************
Filter.h: class for filtering
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
#ifndef _MYFILTER_H_
#define _MYFILTER_H_

#include <vector>
#include <iostream>
#include <gslwrap/vector_double.h>
#include <gslwrap/matrix_double.h>
#include <gslwrap/matrix_vector_operators.h>
#include "StatFactory.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;
namespace filter
{
  
  enum FILTER_TYPE { median=0, var=1 , min=2, max=3, sum=4, mean=5, minmax=6, dilate=7, erode=8, close=9, open=10, homog=11, sobelx=12, sobely=13, sobelxy=14, sobelyx=-14, smooth=15, density=16, majority=17, mixed=18, smoothnodata=19, threshold=20, ismin=21, ismax=22, heterog=23, order=24, stdev=25};

class Filter
{
public:
  Filter(void);
  Filter(const vector<double> &taps);
  virtual ~Filter(){};
  FILTER_TYPE getFilterType(const std::string filterType){
    std::map<std::string, FILTER_TYPE> m_filterMap;
    initMap(m_filterMap);
    return m_filterMap[filterType];
  };
  void setTaps(const vector<double> &taps);
  void pushClass(short theClass=1){m_class.push_back(theClass);};
  void pushMask(short theMask=0){m_mask.push_back(theMask);};
  template<class T> void doit(const vector<T>& input, vector<T>& output, int down=1, int offset=0);
  template<class T> void doit(T* input, int inputSize, vector<T>& output, int down=1, int offset=0);
  template<class T> void morphology(const vector<T>& input, vector<T>& output, const std::string& method, int dim, short down=1, int offset=0, bool verbose=0);
  void morphology(const ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dim, short down=1, int offset=0);
  void doit(const ImgReaderGdal& input, ImgWriterGdal& output, short down=1, int offset=0);

  template<class T> double applySrf(const vector<double> &wavelengthIn, const vector<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, T& output, double delta=1.0, bool normalize=false, bool verbose=false);
  template<class T> double applySrf(const vector<double> &wavelengthIn, const Vector2d<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, vector<T>& output, double delta=1.0, bool normalize=false, int down=1, bool verbose=false);

  // void applySrf(const vector<double> &wavelengthIn, const ImgReaderGdal& input, const vector< Vector2d<double> > &srf, const std::string& interpolationType, ImgWriterGdal& output, bool verbose=false);
  template<class T> void applyFwhm(const vector<double> &wavelengthIn, const vector<T>& input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, vector<T>& output, bool verbose=false);
  template<class T> void applyFwhm(const vector<double> &wavelengthIn, const Vector2d<T>& input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, Vector2d<T>& output, int down=1, bool verbose=false);
  // void applyFwhm(const vector<double> &wavelengthIn, const ImgReaderGdal& input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, ImgWriterGdal& output, bool verbose=false);
// int fir(double* input, int nbandIn, vector<double>& output, int startBand, const string& wavelength, const string& fwhm, bool verbose);
// int fir(const vector<double>&input, vector<double>& output, int startBand, double fwhm, int ntaps, int down, int offset, bool verbose);
// int fir(double* input, int nbandIn, vector<double>& output, int startBand, double fwhm, int ntaps, int down, int offset, bool verbose);

private:
  static void initMap(std::map<std::string, FILTER_TYPE>& m_filterMap){
    //initialize selMap
    m_filterMap["stdev"]=filter::stdev;
    m_filterMap["var"]=filter::var;
    m_filterMap["min"]=filter::min;
    m_filterMap["max"]=filter::max;
    m_filterMap["sum"]=filter::sum;
    m_filterMap["mean"]=filter::mean;
    m_filterMap["minmax"]=filter::minmax;
    m_filterMap["dilate"]=filter::dilate;
    m_filterMap["erode"]=filter::erode;
    m_filterMap["close"]=filter::close;
    m_filterMap["open"]=filter::open;
    m_filterMap["homog"]=filter::homog;
    m_filterMap["sobelx"]=filter::sobelx;
    m_filterMap["sobely"]=filter::sobely;
    m_filterMap["sobelxy"]=filter::sobelxy;
    m_filterMap["sobelyx"]=filter::sobelyx;
    m_filterMap["smooth"]=filter::smooth;
    m_filterMap["density"]=filter::density;
    m_filterMap["majority"]=filter::majority;
    m_filterMap["mixed"]=filter::mixed;
    m_filterMap["smoothnodata"]=filter::smoothnodata;
    m_filterMap["threshold"]=filter::threshold;
    m_filterMap["ismin"]=filter::ismin;
    m_filterMap["ismax"]=filter::ismax;
    m_filterMap["heterog"]=filter::heterog;
    m_filterMap["order"]=filter::order;
    m_filterMap["median"]=filter::median;
  }
  vector<double> m_taps;
  vector<short> m_class;
  vector<short> m_mask;
};

//input[band], output
//returns wavelength for which srf is maximum
  template<class T> double Filter::applySrf(const vector<double> &wavelengthIn, const vector<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, T& output, double delta, bool normalize, bool verbose)
{  
  assert(srf.size()==2);//[0]: wavelength, [1]: response function
  int nband=srf[0].size(); 
  double start=floor(wavelengthIn[0]);
  double end=ceil(wavelengthIn.back());
  if(verbose)
    std::cout << "wavelengths in [" << start << "," << end << "]" << std::endl << std::flush;

  statfactory::StatFactory stat;

  gsl_interp_accel *acc;
  stat.allocAcc(acc);
  gsl_spline *spline;
  stat.getSpline(interpolationType,nband,spline);
  stat.initSpline(spline,&(srf[0][0]),&(srf[1][0]),nband);
  // gsl_interp_accel *acc=gsl_interp_accel_alloc();
  // gsl_spline *spline=gsl_spline_alloc(gsl_interp_linear,nband);
  // gsl_spline_init(spline,&(srf[0][0]),&(srf[1][0]),nband);
  if(verbose)
    std::cout << "calculating norm of srf" << std::endl << std::flush;
  double norm=0;
  norm=gsl_spline_eval_integ(spline,srf[0].front(),srf[0].back(),acc);
  if(verbose)
    std::cout << "norm of srf: " << norm << std::endl << std::flush;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);  
  //interpolate input and srf to delta
  
  vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);

  if(verbose)
    std::cout << "interpolate wavelengths to " << wavelength_fine.size() << " entries " << std::endl;
  vector<double> srf_fine;//spectral response function, interpolated for wavelength_fine

  stat.interpolateUp(srf[0],srf[1],wavelength_fine,interpolationType,srf_fine,verbose);
  assert(srf_fine.size()==wavelength_fine.size());

  gsl_interp_accel *accOut;
  stat.allocAcc(accOut);
  gsl_spline *splineOut;
  stat.getSpline(interpolationType,wavelength_fine.size(),splineOut);
  assert(splineOut);

  assert(wavelengthIn.size()==input.size());
  vector<double> input_fine;
  vector<double> product(wavelength_fine.size());
  stat.interpolateUp(wavelengthIn,input,wavelength_fine,interpolationType,input_fine,verbose);

  if(verbose)
    std::cout << "input_fine.size(): " << input_fine.size() << std::endl;
  for(int iband=0;iband<input_fine.size();++iband)
    product[iband]=input_fine[iband]*srf_fine[iband];

  assert(input_fine.size()==srf_fine.size());
  assert(input_fine.size()==wavelength_fine.size());
  stat.initSpline(splineOut,&(wavelength_fine[0]),&(product[0]),wavelength_fine.size());
  if(normalize)
    output=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
  else
    output=gsl_spline_eval_integ(splineOut,start,end,accOut);
  gsl_spline_free(splineOut);
  gsl_interp_accel_free(accOut);

  double maxResponse=0;
  int maxIndex=0;
  for(int index=0;index<srf[1].size();++index){
    if(maxResponse<srf[1][index]){
      maxResponse=srf[1][index];
      maxIndex=index;
    }
  }
  return(srf[0][maxIndex]);
}

//input[band][sample], output[sample]
//returns wavelength for which srf is maximum
  template<class T> double Filter::applySrf(const vector<double> &wavelengthIn, const Vector2d<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, vector<T>& output, double delta, bool normalize, int down, bool verbose)
{  
  assert(srf.size()==2);//[0]: wavelength, [1]: response function
  int nband=srf[0].size(); 
  unsigned int nsample=input[0].size();
  output.resize((nsample+down-1)/down);
  double start=floor(wavelengthIn[0]);
  double end=ceil(wavelengthIn.back());
  if(verbose)
    std::cout << "wavelengths in [" << start << "," << end << "]" << std::endl << std::flush;

  statfactory::StatFactory stat;

  gsl_interp_accel *acc;
  stat.allocAcc(acc);
  gsl_spline *spline;
  stat.getSpline(interpolationType,nband,spline);
  stat.initSpline(spline,&(srf[0][0]),&(srf[1][0]),nband);
  // gsl_interp_accel *acc=gsl_interp_accel_alloc();
  // gsl_spline *spline=gsl_spline_alloc(gsl_interp_linear,nband);
  // gsl_spline_init(spline,&(srf[0][0]),&(srf[1][0]),nband);
  if(verbose)
    std::cout << "calculating norm of srf" << std::endl << std::flush;
  double norm=0;
  norm=gsl_spline_eval_integ(spline,srf[0].front(),srf[0].back(),acc);
  if(verbose)
    std::cout << "norm of srf: " << norm << std::endl << std::flush;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);  
  //interpolate input and srf to delta
  
  vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);

  if(verbose)
    std::cout << "interpolate wavelengths to " << wavelength_fine.size() << " entries " << std::endl;
  vector<double> srf_fine;//spectral response function, interpolated for wavelength_fine

  stat.interpolateUp(srf[0],srf[1],wavelength_fine,interpolationType,srf_fine,verbose);
  assert(srf_fine.size()==wavelength_fine.size());

  gsl_interp_accel *accOut;
  stat.allocAcc(accOut);
  gsl_spline *splineOut;
  stat.getSpline(interpolationType,wavelength_fine.size(),splineOut);
  assert(splineOut);

  for(int isample=0;isample<nsample;++isample){
    if((isample+1+down/2)%down)
      continue;
    vector<T> inputValues;
    input.selectCol(isample,inputValues);
    assert(wavelengthIn.size()==inputValues.size());
    vector<double> input_fine;
    vector<double> product(wavelength_fine.size());
    stat.interpolateUp(wavelengthIn,inputValues,wavelength_fine,interpolationType,input_fine,verbose);

    for(int iband=0;iband<input_fine.size();++iband)
      product[iband]=input_fine[iband]*srf_fine[iband];

    assert(input_fine.size()==srf_fine.size());
    assert(input_fine.size()==wavelength_fine.size());
    stat.initSpline(splineOut,&(wavelength_fine[0]),&(product[0]),wavelength_fine.size());
    //hiero
    if(normalize)
      output[isample/down]=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
    else
      output[isample/down]=gsl_spline_eval_integ(splineOut,start,end,accOut);
  }
  gsl_spline_free(splineOut);
  gsl_interp_accel_free(accOut);

  double maxResponse=0;
  int maxIndex=0;
  for(int index=0;index<srf[1].size();++index){
    if(maxResponse<srf[1][index]){
      maxResponse=srf[1][index];
      maxIndex=index;
    }
  }
  return(srf[0][maxIndex]);
}

template<class T> void Filter::applyFwhm(const vector<double> &wavelengthIn, const vector<T>& input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, vector<T>& output, bool verbose){
  double delta=1;//1 nm resolution
  vector<double> stddev(fwhm.size());
  for(int index=0;index<fwhm.size();++index)
    stddev[index]=fwhm[index]/2.0/sqrt(2*log(2.0));//http://mathworld.wolfram.com/FullWidthatHalfMaximum.html
  assert(wavelengthOut.size()==fwhm.size());
  assert(wavelengthIn.size()==input.size());
  assert(wavelengthIn[0]<wavelengthOut[0]);
  assert(wavelengthIn.back()>wavelengthOut.back());
  statfactory::StatFactory stat;
  vector<double> input_fine;
  vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);
  if(verbose){
    for(int index=0;index<wavelength_fine.size();++index)
      std::cout << " " << wavelength_fine[index];
    std::cout << std::endl;
    std::cout << "interpolate input wavelength to " << delta << " nm resolution (size=" << wavelength_fine.size() << ")" << std::endl;
  }
  stat.interpolateUp(wavelengthIn,input,wavelength_fine,interpolationType,input_fine,verbose);
  int nbandIn=wavelength_fine.size();
    
  int nbandOut=wavelengthOut.size();
  output.resize(nbandOut);
  gsl::matrix tf(nbandIn,nbandOut);
  for(int indexOut=0;indexOut<nbandOut;++indexOut){
    double norm=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn){
      tf(indexIn,indexOut)=
        exp((wavelengthOut[indexOut]-wavelength_fine[indexIn])
            *(wavelength_fine[indexIn]-wavelengthOut[indexOut])
            /2.0/stddev[indexOut]
            /stddev[indexOut]);
      tf(indexIn,indexOut)/=sqrt(2.0*M_PI);
      tf(indexIn,indexOut)/=stddev[indexOut];
      norm+=tf(indexIn,indexOut);
    }
    output[indexOut]=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn)
      output[indexOut]+=input_fine[indexIn]*tf(indexIn,indexOut)/norm;
  }
}


  //input[inBand][sample], output[outBand][sample]
  template<class T> void Filter::applyFwhm(const vector<double> &wavelengthIn, const Vector2d<T>& input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, Vector2d<T>& output, int down, bool verbose){
  double delta=1;//1 nm resolution
  vector<double> stddev(fwhm.size());
  for(int index=0;index<fwhm.size();++index)
    stddev[index]=fwhm[index]/2.0/sqrt(2*log(2.0));//http://mathworld.wolfram.com/FullWidthatHalfMaximum.html
  statfactory::StatFactory stat;
  vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);
  assert(wavelengthOut.size()==fwhm.size());
  assert(wavelengthIn[0]<wavelengthOut[0]);
  assert(wavelengthIn.back()>wavelengthOut.back());
  if(verbose){
    for(int index=0;index<wavelength_fine.size();++index)
      std::cout << " " << wavelength_fine[index];
    std::cout << std::endl;
    std::cout << "interpolate input wavelength to " << delta << " nm resolution (size=" << wavelength_fine.size() << ")" << std::endl;
  }
  int nbandIn=wavelength_fine.size();
  int nbandOut=wavelengthOut.size();
  output.resize(nbandOut,(input[0].size()+down-1)/down);

  gsl::matrix tf(nbandIn,nbandOut);
  vector<double> norm(nbandOut);
  for(int indexOut=0;indexOut<nbandOut;++indexOut){
    norm[indexOut]=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn){
      tf(indexIn,indexOut)=
        exp((wavelengthOut[indexOut]-wavelength_fine[indexIn])
            *(wavelength_fine[indexIn]-wavelengthOut[indexOut])
            /2.0/stddev[indexOut]
            /stddev[indexOut]);
      tf(indexIn,indexOut)/=sqrt(2.0*M_PI);
      tf(indexIn,indexOut)/=stddev[indexOut];
      norm[indexOut]+=tf(indexIn,indexOut);
    }
  }

  for(int isample=0;isample<input[0].size();++isample){
    if((isample+1+down/2)%down)
      continue;
    vector<T> inputValues;
    input.selectCol(isample,inputValues);
    assert(wavelengthIn.size()==inputValues.size());
    for(int indexOut=0;indexOut<nbandOut;++indexOut){
      vector<double> input_fine;
      stat.interpolateUp(wavelengthIn,inputValues,wavelength_fine,interpolationType,input_fine,verbose);
      output[indexOut][(isample+down-1)/down]=0;
      for(int indexIn=0;indexIn<nbandIn;++indexIn){
        output[indexOut][(isample+down-1)/down]+=input_fine[indexIn]*tf(indexIn,indexOut)/norm[indexOut];
      }
    }
  }
}

template<class T> void Filter::doit(const vector<T>& input, vector<T>& output, int down, int offset)
{
  output.resize((input.size()-offset+down-1)/down);
  int i=0;
  //start: extend input with mirrored version of itself
  for(i=offset;i<m_taps.size()/2;++i){
    if((i-offset)%down)
      continue;
    output[(i-offset+down-1)/down]=m_taps[m_taps.size()/2]*input[i];
    for(int t=1;t<=m_taps.size()/2;++t)
      output[(i-offset+down-1)/down]+=(m_taps[m_taps.size()/2+t]+m_taps[m_taps.size()/2-t])*input[i+t];
  }
  //main
  for(i=offset+m_taps.size()/2;i<input.size()-m_taps.size()/2;++i){
    if((i-offset)%down)
      continue;
    T leaveOut=(*(m_taps.begin()))*input[i-m_taps.size()/2];
    T include=(m_taps.back())*input[i+m_taps.size()/2];
    output[(i-offset+down-1)/down]=0;
    for(int t=0;t<m_taps.size();++t)
      output[(i-offset+down-1)/down]+=input[i-m_taps.size()/2+t]*m_taps[t];
  }
  //end: extend input with mirrored version of itself
  for(i=input.size()-m_taps.size()/2;i<input.size();++i){
    if((i-offset)%down)
      continue;
    output[(i-offset+down-1)/down]=m_taps[m_taps.size()/2]*input[i];
    for(int t=1;t<=m_taps.size()/2;++t)
      output[(i-offset+down-1)/down]+=(m_taps[m_taps.size()/2+t]+m_taps[m_taps.size()/2-t])*input[i-t];
  }
}

template<class T> void Filter::morphology(const vector<T>& input, vector<T>& output, const std::string& method, int dim, short down, int offset, bool verbose)
{
  assert(dim);
  output.resize((input.size()-offset+down-1)/down);
  int i=0;
  statfactory::StatFactory stat;
  vector<T> statBuffer;
  short binValue=0;
  //start: extend input with mirrored version of itself
  for(i=offset;i<dim/2;++i){
    binValue=0;
    for(int iclass=0;iclass<m_class.size();++iclass){
      if(input[i]==m_class[iclass]){
        binValue=m_class[0];
        break;
      }
    }
    if(m_class.size())
      statBuffer.push_back(binValue);
    else
      statBuffer.push_back(input[i]);
    for(int t=1;t<=dim/2;++t){
      binValue=0;
      for(int iclass=0;iclass<m_class.size();++iclass){
        if(input[i+t]==m_class[iclass]){
          binValue=m_class[0];
          break;
        }
      }
      if(m_class.size()){
        statBuffer.push_back(binValue);
        statBuffer.push_back(binValue);
      }
      else{
        statBuffer.push_back(input[i+t]);
        statBuffer.push_back(input[i+t]);
      }
    }
    assert(statBuffer.size()==dim);
    if((i-offset)%down){
      statBuffer.clear();
      continue;
    }
    switch(getFilterType(method)){
    case(filter::dilate):
      output[(i-offset+down-1)/down]=stat.max(statBuffer);
      break;
    case(filter::erode):
      output[(i-offset+down-1)/down]=stat.min(statBuffer);
      break;
    default:
      string errorString="method not supported";
      throw(errorString);
      break;
    }
    if(verbose){
      cout << "buffer: ";
      for(int ibuf=0;ibuf<statBuffer.size();++ibuf)
        cout << statBuffer[ibuf] << " ";
      cout << "->" << output[(i-offset+down-1)/down] << endl;
    }
  }
  //main
  statBuffer.clear();
  for(i=offset+dim/2;i<input.size()-dim/2;++i){
    binValue=0;
    for(int t=0;t<dim;++t){
      for(int iclass=0;iclass<m_class.size();++iclass){
        if(input[i-dim/2+t]==m_class[iclass]){
          binValue=m_class[0];
          break;
        }
      }
      if(m_class.size())
        statBuffer.push_back(binValue);
      else
        statBuffer.push_back(input[i-dim/2+t]);
    }
    assert(statBuffer.size()==dim);
    if((i-offset)%down){
      statBuffer.clear();
      continue;
    }
    switch(getFilterType(method)){
    case(filter::dilate):
      output[(i-offset+down-1)/down]=stat.max(statBuffer);
      break;
    case(filter::erode):
      output[(i-offset+down-1)/down]=stat.min(statBuffer);
      break;
    default:
      string errorString="method not supported";
      throw(errorString);
      break;
    }
    if(verbose){
      cout << "buffer: ";
      for(int ibuf=0;ibuf<statBuffer.size();++ibuf)
        cout << statBuffer[ibuf] << " ";
      cout << "->" << output[(i-offset+down-1)/down] << endl;
    }
    statBuffer.clear();
  }
  //end: extend input with mirrored version of itself
  for(i=input.size()-dim/2;i<input.size();++i){
      binValue=0;
      for(int iclass=0;iclass<m_class.size();++iclass){
        if(input[i]==m_class[iclass]){
          binValue=m_class[0];
          break;
        }
      }
      if(m_class.size())
        statBuffer.push_back(binValue);
      else
        statBuffer.push_back(input[i]);
      for(int t=1;t<=dim/2;++t){
        binValue=0;
        for(int iclass=0;iclass<m_class.size();++iclass){
          if(input[i-t]==m_class[iclass]){
            binValue=m_class[0];
            break;
          }
        }
        if(m_class.size()){
          statBuffer.push_back(binValue);
          statBuffer.push_back(binValue);
        }
        else{
          statBuffer.push_back(input[i-t]);
          statBuffer.push_back(input[i-t]);
        }
      }
    if((i-offset)%down){
      statBuffer.clear();
      continue;
    }
    switch(getFilterType(method)){
    case(filter::dilate):
      output[(i-offset+down-1)/down]=stat.max(statBuffer);
      break;
    case(filter::erode):
      output[(i-offset+down-1)/down]=stat.min(statBuffer);
      break;
    default:
      string errorString="method not supported";
      throw(errorString);
      break;
    }
    if(verbose){
      cout << "buffer: ";
      for(int ibuf=0;ibuf<statBuffer.size();++ibuf)
        cout << statBuffer[ibuf] << " ";
      cout << "->" << output[(i-offset+down-1)/down] << endl;
    }
  }
}

template<class T> void Filter::doit(T* input, int inputSize, vector<T>& output, int down, int offset)
{
  output.resize((inputSize-offset+down-1)/down);
  int i=0;
  //start: extend input with mirrored version of itself
  for(i=offset;i<m_taps.size()/2;++i){
    if((i-offset)%down)
      continue;
    output[(i-offset+down-1)/down]=m_taps[m_taps.size()/2]*input[i];
    for(int t=1;t<=m_taps.size()/2;++t)
      output[(i-offset+down-1)/down]+=(m_taps[m_taps.size()/2+t]+m_taps[m_taps.size()/2-t])*input[i+t];
  }
  //main
  for(i=offset+m_taps.size()/2;i<inputSize-m_taps.size()/2;++i){
    if((i-offset)%down)
      continue;
    T leaveOut=(*(m_taps.begin()))*input[i-m_taps.size()/2];
    T include=(m_taps.back())*input[i+m_taps.size()/2];
    output[(i-offset+down-1)/down]=0;
    for(int t=0;t<m_taps.size();++t)
      output[(i-offset+down-1)/down]+=input[i-m_taps.size()/2+t]*m_taps[t];
  }
  //end: extend input with mirrored version of itself
  for(i=inputSize-m_taps.size()/2;i<inputSize;++i){
    if((i-offset)%down)
      continue;
    output[(i-offset+down-1)/down]=m_taps[m_taps.size()/2]*input[i];
    for(int t=1;t<=m_taps.size()/2;++t)
      output[(i-offset+down-1)/down]+=(m_taps[m_taps.size()/2+t]+m_taps[m_taps.size()/2-t])*input[i-t];
  }
}
}

#endif /* _MYFILTER_H_ */

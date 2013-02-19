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
  
class Filter
{
public:
  enum Type { MEDIAN=0, VAR=1 , MIN=2, MAX=3, SUM=4, MEAN=5, MINMAX=6, DILATE=7, ERODE=8, CLOSE=9, OPEN=10, HOMOG=11, SOBELX=12, SOBELY=13, SOBELXY=14, SMOOTH=15, DENSITY=16, MAJORITY=17, MIXED=18};
  Filter(void);
  Filter(const vector<double> &taps);
  virtual ~Filter(){};
  void setTaps(const vector<double> &taps);
  void pushClass(short theClass=1){m_class.push_back(theClass);};
  void pushMask(short theMask=0){m_mask.push_back(theMask);};
  template<class T> void doit(const vector<T>& input, vector<T>& output, int down=1, int offset=0);
  template<class T> void doit(T* input, int inputSize, vector<T>& output, int down=1, int offset=0);
  template<class T> void morphology(const vector<T>& input, vector<T>& output, int method, int dim, short down=1, int offset=0, bool verbose=0);
  void morphology(const ImgReaderGdal& input, ImgWriterGdal& output, int method, int dim, short down=1, int offset=0);
  void doit(const ImgReaderGdal& input, ImgWriterGdal& output, short down=1, int offset=0);

  template<class T> void applySrf(const Vector2d<T>& input, const Vector2d<double>& srf, Vector2d<T>& output, double delta=1, bool normalize=false, double centreWavelength=0, bool verbose=false);
  template<class T> void applyFwhm(const vector<double>& input, const vector<double> &wavelengthIn, vector<double>& output, const vector<double> &wavelengthOut, const vector<double> &fwhm, bool verbose=false);
  void applyFwhm(const ImgReaderGdal& input, ImgWriterGdal& output, const vector<double> &wavelengthIn, const vector<double> &wavelengthOut, const vector<double> &fwhm, bool verbose=false);
// int fir(double* input, int nbandIn, vector<double>& output, int startBand, const string& wavelength, const string& fwhm, bool verbose);
// int fir(const vector<double>&input, vector<double>& output, int startBand, double fwhm, int ntaps, int down, int offset, bool verbose);
// int fir(double* input, int nbandIn, vector<double>& output, int startBand, double fwhm, int ntaps, int down, int offset, bool verbose);

private:
  vector<double> m_taps;
  vector<short> m_class;
  vector<short> m_mask;
};

  template<class T> void Filter::applySrf(const Vector2d<T>& input, const Vector2d<double>& srf, Vector2d<T>& output, double delta, bool normalize, double centreWavelength, bool verbose)
{  
  output.resize(input.size());
  assert(srf.size()==2);//[0]: wavelength, [1]: response function
  assert(input.size()>1);//[0]: wavelength, [1],[2],...: value
  double start=floor(input[0][0]);
  double end=ceil(input[0].back());
  Vector2d<double> product(input.size());  
  assert(input.size());
  assert(input[0].size()>1);
  int nband=srf[0].size();  
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  gsl_spline *spline=gsl_spline_alloc(gsl_interp_linear,nband);
  gsl_spline_init(spline,&(srf[0][0]),&(srf[1][0]),nband);
//   double norm=gsl_spline_eval_integ(spline,start,end,acc);
  double norm=gsl_spline_eval_integ(spline,srf[0].front(),srf[0].back(),acc);
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);  
  //interpolate input and srf to delta
  statfactory::StatFactory stat;
  Vector2d<double> input_d;
  Vector2d<double> srf_d;
  stat.interpolateUp(input,input_d,start,end,delta);
  stat.interpolateUp(srf,srf_d,start,end,delta);
  nband=input_d[0].size();
  if(verbose)
    cout << "number of interpolated bands: " << nband << endl;
  for(int isample=0;isample<input_d.size();++isample){
    product[isample].resize(nband);
    for(int iband=0;iband<nband;++iband){
      if(!isample)
	product[isample][iband]=input_d[isample][iband];
      else{
// 	if(verbose&&isample==1)
// 	  cout << srf_d[0][iband] << " " << srf_d[1][iband] << endl;
	product[isample][iband]=input_d[isample][iband]*srf_d[1][iband];
      }
    }
  }
  output[0].resize(1);
  if(centreWavelength)
    output[0][0]=centreWavelength;
  else{
    double maxResponse=0;
    int maxIndex=0;
    for(int index=0;index<srf[1].size();++index){
      if(maxResponse<srf[1][index]){
	maxResponse=srf[1][index];
	maxIndex=index;
      }
    }
    output[0][0]=srf[0][maxIndex];
  }
  for(int isample=1;isample<input_d.size();++isample){
    output[isample].resize(1);    
    gsl_interp_accel *acc=gsl_interp_accel_alloc();
    gsl_spline *spline=gsl_spline_alloc(gsl_interp_linear,nband);
    gsl_spline_init(spline,&(product[0][0]),&(product[isample][0]),nband);
    if(normalize)
      output[isample][0]=gsl_spline_eval_integ(spline,start,end,acc)/norm;
    else
      output[isample][0]=gsl_spline_eval_integ(spline,start,end,acc);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
  }
}

  template<class T> void Filter::applyFwhm(const vector<double>& input, const vector<double> &wavelengthIn, vector<double>& output, const vector<double> &wavelengthOut, const vector<double> &fwhm, bool verbose){
  double delta=1;//1 nm resolution
  vector<double> stddev(fwhm.size());
  for(int index=0;index<fwhm.size();++index)
    stddev[index]=fwhm[index]/2.0/sqrt(2*log(2.0));//http://mathworld.wolfram.com/FullWidthatHalfMaximum.html
  assert(wavelengthOut.size()==fwhm.size());
  assert(wavelengthIn.size()==input.size());
  //todo densify input to 1 nm
  assert(wavelengthIn[0]<wavelengthOut[0]);
  assert(wavelengthIn.back()<wavelengthOut.back());
  statfactory::StatFactory stat;
  Vector2d<double> input_course(1);
  input_course[0]=input;
  Vector2d<double> input_fine;
  vector<double> wavelength_fine;
  stat.interpolateUp(input_course,wavelengthIn,input_fine,wavelength_fine,wavelengthIn[0],wavelengthIn.back(),delta);
  int nbandIn=wavelength_fine.size();
  int nbandOut=wavelengthOut.size();
  gsl::matrix tf(nbandIn,nbandOut);
  for(int indexOut=0;indexOut<nbandOut;++indexOut){
    for(int indexIn=0;indexIn<wavelengthIn.size();++indexIn){
      tf(indexIn,indexOut)=
	exp((wavelengthOut[indexOut]-wavelengthIn[indexIn])
	    *(wavelengthIn[indexIn]-wavelengthOut[indexOut])
	    /2.0/stddev[indexOut]
	    /stddev[indexOut]);
      tf(indexIn,indexOut)/=sqrt(2.0*M_PI);
      tf(indexIn,indexOut)/=stddev[indexOut];
    }
    double norm=exp(tf.LU_lndet());//(tf.Column(indexOut+1)).NormFrobenius();
    if(norm)
      tf.get_col(indexOut+1)/=norm;
  }
  //create filtered vector
  output.resize(nbandOut);
  for(int indexOut=0;indexOut<nbandOut;++indexOut){
    output[indexOut]=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn)
      output[indexOut]+=input_fine[0][indexIn]*tf(indexIn,indexOut);
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

template<class T> void Filter::morphology(const vector<T>& input, vector<T>& output, int method, int dim, short down, int offset, bool verbose)
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
    switch(method){
    case(DILATE):
      output[(i-offset+down-1)/down]=stat.max(statBuffer);
      break;
    case(ERODE):
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
    switch(method){
    case(DILATE):
      output[(i-offset+down-1)/down]=stat.max(statBuffer);
      break;
    case(ERODE):
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
    switch(method){
    case(DILATE):
      output[(i-offset+down-1)/down]=stat.max(statBuffer);
      break;
    case(ERODE):
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

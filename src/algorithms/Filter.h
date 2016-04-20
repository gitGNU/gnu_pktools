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
extern "C" {
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
}
#include "StatFactory.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

namespace filter
{
  
  enum FILTER_TYPE { median=0, var=1 , min=2, max=3, sum=4, mean=5, minmax=6, dilate=7, erode=8, close=9, open=10, homog=11, sobelx=12, sobely=13, sobelxy=14, sobelyx=-14, smooth=15, density=16, mode=17, mixed=18, smoothnodata=19, threshold=20, ismin=21, ismax=22, heterog=23, order=24, stdev=25, dwt=26, dwti=27, dwt_cut=28, dwt_cut_from=29, savgolay=30, percentile=31, nvalid=32};

   enum PADDING { symmetric=0, replicate=1, circular=2, zero=3};

class Filter
{
public:
  Filter(void);
  Filter(const std::vector<double> &taps);
  virtual ~Filter(){};

  void setPadding(const std::string& padString){
    m_padding=padString;
  };

  static const gsl_wavelet_type* getWaveletType(const std::string waveletType){
    if(waveletType=="daubechies") return(gsl_wavelet_daubechies);
    if(waveletType=="daubechies_centered") return(gsl_wavelet_daubechies_centered);
    if(waveletType=="haar") return(gsl_wavelet_haar);
    if(waveletType=="haar_centered") return(gsl_wavelet_haar_centered);
    if(waveletType=="bspline") return(gsl_wavelet_bspline);
    if(waveletType=="bspline_centered") return(gsl_wavelet_bspline_centered);
  }
  static FILTER_TYPE getFilterType(const std::string filterType){
    std::map<std::string, FILTER_TYPE> m_filterMap;
    initFilterMap(m_filterMap);
    return m_filterMap[filterType];
  };

  void setTaps(const std::vector<double> &taps, bool normalize=true);
  void pushClass(short theClass=1){m_class.push_back(theClass);};
  void pushMask(short theMask=0){m_mask.push_back(theMask);};
  unsigned int pushNoDataValue(double noDataValue);
  unsigned int setNoDataValues(std::vector<double> vnodata);
  void pushThreshold(double theThreshold){m_threshold.push_back(theThreshold);};
  void setThresholds(const std::vector<double>& theThresholds){m_threshold=theThresholds;};
  template<class T> void filter(const std::vector<T>& input, std::vector<T>& output);
  template<class T> void filter(const std::vector<T>& input, std::vector<T>& output, const std::string& method, int dim);
  template<class T> void smooth(const std::vector<T>& input, std::vector<T>& output, short dim);
  template<class T> void smoothNoData(const std::vector<T>& input, const std::string& interpolationType, std::vector<T>& output);
  template<class T> void filter(T* input, int inputSize, std::vector<T>& output);
  template<class T> void smooth(T* input, int inputSize, std::vector<T>& output, short dim);
  //template<class T> void morphology(const std::vector<T>& input, std::vector<T>& output, const std::string& method, int dim, bool verbose=false);
  void morphology(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dim, short verbose=0);
  void filter(ImgReaderGdal& input, ImgWriterGdal& output);
  void stat(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method);
  void stats(ImgReaderGdal& input, ImgWriterGdal& output, const std::vector<std::string >& methods);
  void filter(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dim);
  void getSavGolayCoefficients(std::vector<double> &c, int np, int nl, int nr, int ld, int m);
  void ludcmp(std::vector<double> &a, std::vector<int> &indx, double &d);
  void lubksb(std::vector<double> &a, std::vector<int> &indx, std::vector<double> &b);
  /* void savgolay(const ImgReaderGdal& input, ImgWriterGdal& output, int np, int nl, int nr, int m); */
  void smooth(ImgReaderGdal& input, ImgWriterGdal& output, short dim);
  void smoothNoData(ImgReaderGdal& input, const std::string& interpolationType, ImgWriterGdal& output);
  double getCentreWavelength(const std::vector<double> &wavelengthIn, const Vector2d<double>& srf, const std::string& interpolationType, double delta=1.0, bool verbose=false);
  template<class T> double applySrf(const std::vector<double> &wavelengthIn, const std::vector<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, T& output, double delta=1.0, bool normalize=false, bool verbose=false);
  template<class T> double applySrf(const std::vector<double> &wavelengthIn, const Vector2d<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, std::vector<T>& output, double delta=1.0, bool normalize=false, int down=1, bool transposeInput=false, bool verbose=false);

  template<class T> void applyFwhm(const std::vector<double> &wavelengthIn, const std::vector<T>& input, const std::vector<double> &wavelengthOut, const std::vector<double> &fwhm, const std::string& interpolationType, std::vector<T>& output, bool verbose=false);
  template<class T> void applyFwhm(const std::vector<double> &wavelengthIn, const Vector2d<T>& input, const std::vector<double> &wavelengthOut, const std::vector<double> &fwhm, const std::string& interpolationType, Vector2d<T>& output, int down=1, bool verbose=false);
  void dwtForward(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& wavelet_type, int family);
  void dwtInverse(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& wavelet_type, int family);
  void dwtCut(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& wavelet_type, int family, double cut);
  void dwtForward(std::vector<double>& data, const std::string& wavelet_type, int family);
  void dwtInverse(std::vector<double>& data, const std::string& wavelet_type, int family);
  void dwtCut(std::vector<double>& data, const std::string& wavelet_type, int family, double cut);
  void dwtCutFrom(ImgReaderGdal& input, ImgWriterGdal& output, const std::string& wavelet_type, int family, int band);

private:

  static void initFilterMap(std::map<std::string, FILTER_TYPE>& m_filterMap){
    //initialize Map
    m_filterMap["dwt"]=filter::dwt;
    m_filterMap["dwti"]=filter::dwti;
    m_filterMap["dwt_cut"]=filter::dwt_cut;
    m_filterMap["dwt_cut_from"]=filter::dwt_cut_from;
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
    m_filterMap["mode"]=filter::mode;
    m_filterMap["mixed"]=filter::mixed;
    m_filterMap["smoothnodata"]=filter::smoothnodata;
    m_filterMap["threshold"]=filter::threshold;
    m_filterMap["ismin"]=filter::ismin;
    m_filterMap["ismax"]=filter::ismax;
    m_filterMap["heterog"]=filter::heterog;
    m_filterMap["order"]=filter::order;
    m_filterMap["nvalid"]=filter::nvalid;
    m_filterMap["median"]=filter::median;
    m_filterMap["savgolay"]=filter::savgolay;
    m_filterMap["percentile"]=filter::percentile;
  }


  static PADDING getPadding(const std::string& padString){
    std::map<std::string, PADDING> padMap;
    padMap["zero"]=filter::zero;
    padMap["symmetric"]=filter::symmetric;
    padMap["replicate"]=filter::replicate;
    padMap["circular"]=filter::circular;
    return(padMap[padString]);
  };

  std::vector<double> m_taps;
  std::vector<short> m_class;
  std::vector<short> m_mask;
  std::string m_padding;
  std::vector<double> m_noDataValues;
  std::vector<double> m_threshold;
};


//input[band], output
//returns wavelength for which srf is maximum
  template<class T> double Filter::applySrf(const std::vector<double> &wavelengthIn, const std::vector<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, T& output, double delta, bool normalize, bool verbose)
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
  if(verbose)
    std::cout << "calculating norm of srf" << std::endl << std::flush;
  double norm=0;
  norm=gsl_spline_eval_integ(spline,srf[0].front(),srf[0].back(),acc);
  if(verbose)
    std::cout << "norm of srf: " << norm << std::endl << std::flush;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);  
  //interpolate input and srf to delta
  
  std::vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);

  if(verbose)
    std::cout << "interpolate wavelengths to " << wavelength_fine.size() << " entries " << std::endl;
  std::vector<double> srf_fine;//spectral response function, interpolated for wavelength_fine

  stat.interpolateUp(srf[0],srf[1],wavelength_fine,interpolationType,srf_fine,verbose);
  assert(srf_fine.size()==wavelength_fine.size());

  gsl_interp_accel *accOut;
  stat.allocAcc(accOut);
  gsl_spline *splineOut;
  stat.getSpline(interpolationType,wavelength_fine.size(),splineOut);
  assert(splineOut);

  assert(wavelengthIn.size()==input.size());
  std::vector<double> input_fine;
  std::vector<double> product(wavelength_fine.size());
  std::vector<double> wavelengthOut(wavelength_fine.size());
  stat.interpolateUp(wavelengthIn,input,wavelength_fine,interpolationType,input_fine,verbose);

  if(verbose)
    std::cout << "input_fine.size(): " << input_fine.size() << std::endl;
  for(int iband=0;iband<input_fine.size();++iband){
    product[iband]=input_fine[iband]*srf_fine[iband];
    wavelengthOut[iband]=wavelength_fine[iband]*srf_fine[iband];
  }

  assert(input_fine.size()==srf_fine.size());
  assert(input_fine.size()==wavelength_fine.size());
  stat.initSpline(splineOut,&(wavelength_fine[0]),&(product[0]),wavelength_fine.size());
  if(normalize)
    output=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
  else
    output=gsl_spline_eval_integ(splineOut,start,end,accOut);

  stat.initSpline(splineOut,&(wavelength_fine[0]),&(wavelengthOut[0]),wavelength_fine.size());
  double centreWavelength=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
  
  gsl_spline_free(splineOut);
  gsl_interp_accel_free(accOut);

  return(centreWavelength);
}

//input[band][sample], output[sample] (if !transposeInput)
//returns wavelength for which srf is maximum
  template<class T> double Filter::applySrf(const std::vector<double> &wavelengthIn, const Vector2d<T>& input, const Vector2d<double>& srf, const std::string& interpolationType, std::vector<T>& output, double delta, bool normalize, int down, bool transposeInput, bool verbose)
{  
  assert(srf.size()==2);//[0]: wavelength, [1]: response function
  int nband=srf[0].size(); 
  unsigned int nsample=(transposeInput)? input.size():input[0].size();
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
  if(verbose)
    std::cout << "calculating norm of srf" << std::endl << std::flush;
  double norm=0;
  norm=gsl_spline_eval_integ(spline,srf[0].front(),srf[0].back(),acc);
  if(verbose)
    std::cout << "norm of srf: " << norm << std::endl << std::flush;
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);  
  //interpolate input and srf to delta
  
  std::vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);

  if(verbose)
    std::cout << "interpolate wavelengths to " << wavelength_fine.size() << " entries " << std::endl;
  std::vector<double> srf_fine;//spectral response function, interpolated for wavelength_fine

  stat.interpolateUp(srf[0],srf[1],wavelength_fine,interpolationType,srf_fine,verbose);
  assert(srf_fine.size()==wavelength_fine.size());

  gsl_interp_accel *accOut;
  stat.allocAcc(accOut);
  gsl_spline *splineOut;
  stat.getSpline(interpolationType,wavelength_fine.size(),splineOut);
  assert(splineOut);

  std::vector<double> wavelengthOut;
  double centreWavelength=0;
  for(int isample=0;isample<nsample;++isample){
    if((isample+1+down/2)%down)
      continue;
    std::vector<T> inputValues;
    if(transposeInput)
      inputValues=input[isample];
    else
      input.selectCol(isample,inputValues);
    assert(wavelengthIn.size()==inputValues.size());
    std::vector<double> input_fine;
    std::vector<double> product(wavelength_fine.size());
    stat.interpolateUp(wavelengthIn,inputValues,wavelength_fine,interpolationType,input_fine,verbose);

    for(int iband=0;iband<input_fine.size();++iband){
      product[iband]=input_fine[iband]*srf_fine[iband];
      if(wavelengthOut.size()<input_fine.size())
        wavelengthOut.push_back(wavelength_fine[iband]*srf_fine[iband]);
    }

    assert(input_fine.size()==srf_fine.size());
    assert(input_fine.size()==wavelength_fine.size());
    stat.initSpline(splineOut,&(wavelength_fine[0]),&(product[0]),wavelength_fine.size());
    if(normalize)
      output[isample/down]=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
    else
      output[isample/down]=gsl_spline_eval_integ(splineOut,start,end,accOut);

    stat.initSpline(splineOut,&(wavelength_fine[0]),&(wavelengthOut[0]),wavelength_fine.size());
    if(centreWavelength>0);
    else
      centreWavelength=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
  }
  gsl_spline_free(splineOut);
  gsl_interp_accel_free(accOut);

  return(centreWavelength);
}

template<class T> void Filter::applyFwhm(const std::vector<double> &wavelengthIn, const std::vector<T>& input, const std::vector<double> &wavelengthOut, const std::vector<double> &fwhm, const std::string& interpolationType, std::vector<T>& output, bool verbose){
  double delta=1;//1 nm resolution
  std::vector<double> stddev(fwhm.size());
  for(int index=0;index<fwhm.size();++index)
    stddev[index]=fwhm[index]/2.0/sqrt(2*log(2.0));//http://mathworld.wolfram.com/FullWidthatHalfMaximum.html
  assert(wavelengthOut.size()==fwhm.size());
  assert(wavelengthIn.size()==input.size());
  assert(wavelengthIn[0]<=wavelengthOut[0]);
  assert(wavelengthIn.back()>=wavelengthOut.back());
  statfactory::StatFactory stat;
  std::vector<double> input_fine;
  std::vector<double> wavelength_fine;
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
  Vector2d<double> tf(nbandIn,nbandOut);
  for(int indexOut=0;indexOut<nbandOut;++indexOut){
    double norm=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn){
      // tf(indexIn,indexOut)=
      tf[indexIn][indexOut]=
        exp((wavelengthOut[indexOut]-wavelength_fine[indexIn])
            *(wavelength_fine[indexIn]-wavelengthOut[indexOut])
            /2.0/stddev[indexOut]
            /stddev[indexOut]);
      tf[indexIn][indexOut]/=sqrt(2.0*M_PI);
      tf[indexIn][indexOut]/=stddev[indexOut];
      norm+=tf[indexIn][indexOut];
    }
    output[indexOut]=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn)
      output[indexOut]+=input_fine[indexIn]*tf[indexIn][indexOut]/norm;
  }
}


  //input[inBand][sample], output[outBand][sample]
  template<class T> void Filter::applyFwhm(const std::vector<double> &wavelengthIn, const Vector2d<T>& input, const std::vector<double> &wavelengthOut, const std::vector<double> &fwhm, const std::string& interpolationType, Vector2d<T>& output, int down, bool verbose){
  double delta=1;//1 nm resolution
  std::vector<double> stddev(fwhm.size());
  for(int index=0;index<fwhm.size();++index)
    stddev[index]=fwhm[index]/2.0/sqrt(2*log(2.0));//http://mathworld.wolfram.com/FullWidthatHalfMaximum.html
  statfactory::StatFactory stat;
  std::vector<double> wavelength_fine;
  for(double win=floor(wavelengthIn[0]);win<=ceil(wavelengthIn.back());win+=delta)
    wavelength_fine.push_back(win);
  assert(wavelengthOut.size()==fwhm.size());
  assert(wavelengthIn[0]<=wavelengthOut[0]);
  assert(wavelengthIn.back()>=wavelengthOut.back());
  if(verbose){
    for(int index=0;index<wavelength_fine.size();++index)
      std::cout << " " << wavelength_fine[index];
    std::cout << std::endl;
    std::cout << "interpolate input wavelength to " << delta << " nm resolution (size=" << wavelength_fine.size() << ")" << std::endl;
  }
  int nbandIn=wavelength_fine.size();
  int nbandOut=wavelengthOut.size();
  output.resize(nbandOut,(input[0].size()+down-1)/down);

  Vector2d<double> tf(nbandIn,nbandOut);
  std::vector<double> norm(nbandOut);
  for(int indexOut=0;indexOut<nbandOut;++indexOut){
    norm[indexOut]=0;
    for(int indexIn=0;indexIn<nbandIn;++indexIn){
      tf[indexIn][indexOut]=
        exp((wavelengthOut[indexOut]-wavelength_fine[indexIn])
            *(wavelength_fine[indexIn]-wavelengthOut[indexOut])
            /2.0/stddev[indexOut]
            /stddev[indexOut]);
      tf[indexIn][indexOut]/=sqrt(2.0*M_PI);
      tf[indexIn][indexOut]/=stddev[indexOut];
      norm[indexOut]+=tf[indexIn][indexOut];
    }
  }

  for(int isample=0;isample<input[0].size();++isample){
    if((isample+1+down/2)%down)
      continue;
    std::vector<T> inputValues;
    input.selectCol(isample,inputValues);
    assert(wavelengthIn.size()==inputValues.size());
    for(int indexOut=0;indexOut<nbandOut;++indexOut){
      std::vector<double> input_fine;
      stat.interpolateUp(wavelengthIn,inputValues,wavelength_fine,interpolationType,input_fine,verbose);
      output[indexOut][(isample+down-1)/down]=0;
      for(int indexIn=0;indexIn<nbandIn;++indexIn){
        output[indexOut][(isample+down-1)/down]+=input_fine[indexIn]*tf[indexIn][indexOut]/norm[indexOut];
      }
    }
  }
}

  template<class T> void Filter::smooth(const std::vector<T>& input, std::vector<T>& output, short dim)
{
  assert(dim>0);
  m_taps.resize(dim);
  for(int itap=0;itap<dim;++itap)
    m_taps[itap]=1.0/dim;
  filter(input,output);
 }

  template<class T> void Filter::smoothNoData(const std::vector<T>& input, const std::string& interpolationType, std::vector<T>& output)
{
  statfactory::StatFactory stat;
  stat.setNoDataValues(m_noDataValues);
  std::vector<double> abscis(input.size());
  for(int i=0;i<abscis.size();++i)
    abscis[i]=i;
  stat.interpolateNoData(abscis,input,interpolationType,output);
 }

template<class T> void Filter::filter(const std::vector<T>& input, std::vector<T>& output)
{
  assert(input.size()>=m_taps.size());
  output.resize(input.size());
  int i=0;
  //start: extend input by padding
  for(i=0;i<m_taps.size()/2;++i){
    //todo:introduce nodata?
    output[i]=m_taps[m_taps.size()/2]*input[i];
    for(int t=1;t<=m_taps.size()/2;++t){
      output[i]+=m_taps[m_taps.size()/2+t]*input[i+t];
      if(i>=t)
	output[i]+=m_taps[m_taps.size()/2-t]*input[i-t];
      else{
	switch(getPadding(m_padding)){
	case(replicate):
	  output[i]+=m_taps[m_taps.size()/2-t]*input[0];
	  break;
	case(circular):
	  output[i]+=m_taps[m_taps.size()/2-t]*input[input.size()+i-t];
	  break;
	case(zero):
	  output[i]+=m_taps[m_taps.size()/2-t]*0;
	  break;
	case(symmetric):
	default:
	  output[i]+=m_taps[m_taps.size()/2-t]*input[t-i];
	  break;
	}
      }
    }
  }
  //main
  for(i=m_taps.size()/2;i<input.size()-m_taps.size()/2;++i){
    //todo:introduce nodata
    T leaveOut=(*(m_taps.begin()))*input[i-m_taps.size()/2];
    T include=(m_taps.back())*input[i+m_taps.size()/2];
    output[i]=0;
    for(int t=0;t<m_taps.size();++t)
      output[i]+=input[i-m_taps.size()/2+t]*m_taps[t];
  }
  //end: extend input by padding
  for(i=input.size()-m_taps.size()/2;i<input.size();++i){
    //todo:introduce nodata?
    output[i]=m_taps[m_taps.size()/2]*input[i];
    //todo:introduce nodata?
    for(int t=1;t<=m_taps.size()/2;++t){
      output[i]+=m_taps[m_taps.size()/2-t]*input[i-t];
      if(i+t<input.size())
	output[i]+=m_taps[m_taps.size()/2+t]*input[i+t];
      else{
	switch(getPadding(m_padding)){
	case(replicate):
	  output[i]+=m_taps[m_taps.size()/2+t]*input.back();
	  break;
	case(circular):
	  output[i]+=m_taps[m_taps.size()/2+t]*input[t-1];
	  break;
	case(zero):
	  output[i]+=m_taps[m_taps.size()/2+t]*0;
	  break;
	case(symmetric):
	default:
	  output[i]+=m_taps[m_taps.size()/2+t]*input[i-t];
	  break;
	}
      }
    //output[i]+=(m_taps[m_taps.size()/2+t]+m_taps[m_taps.size()/2-t])*input[i-t];
    }
  }
}

//todo: filling statBuffer can be optimized (no need to clear and fill entire buffer, just push back new value...)
 template<class T> void Filter::filter(const std::vector<T>& input, std::vector<T>& output, const std::string& method, int dim)
{
  bool verbose=false;
  assert(dim);
  output.resize(input.size());
  int i=0;
  statfactory::StatFactory stat;
  stat.setNoDataValues(m_noDataValues);
  std::vector<T> statBuffer;
  short binValue=0;
  //start: extend input by padding
  for(i=0;i<dim/2;++i){
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
      T theValue=input[i+t];
      for(int iclass=0;iclass<m_class.size();++iclass){
        if(theValue==m_class[iclass]){
          binValue=m_class[0];
          break;
        }
      }
      if(m_class.size())
	statBuffer.push_back(binValue);
      else
	statBuffer.push_back(theValue);

      if(i>=t){
	theValue=input[i-t];
      }
      else{
	switch(getPadding(m_padding)){
	case(replicate):
	  theValue=input[0];
	  break;
	case(circular):
	  theValue=input[input.size()+i-t];
	  break;
	case(zero):
	  theValue=0;
	  break;
	case(symmetric):
	default:
	  theValue=input[t-i];
	  break;
	}
      }
      for(int iclass=0;iclass<m_class.size();++iclass){
        if(theValue==m_class[iclass]){
          binValue=m_class[0];
          break;
        }
      }
      if(m_class.size())
	statBuffer.push_back(binValue);
      else
	statBuffer.push_back(theValue);
    }

    switch(getFilterType(method)){
    case(filter::nvalid):
      output[i]=stat.nvalid(statBuffer);
      break;
    case(filter::median):
      output[i]=stat.median(statBuffer);
      break;
    case(filter::min):
    case(filter::erode):
      output[i]=stat.mymin(statBuffer);
      break;
    case(filter::max):
    case(filter::dilate):
      output[i]=stat.mymax(statBuffer);
      break;
    case(filter::sum):
      output[i]=sqrt(stat.sum(statBuffer));
      break;
    case(filter::var):
      output[i]=stat.var(statBuffer);
      break;
    case(filter::stdev):
      output[i]=sqrt(stat.var(statBuffer));
      break;
    case(filter::mean):
      output[i]=stat.mean(statBuffer);
      break;
    case(filter::percentile):
      assert(m_threshold.size());
      output[i]=stat.percentile(statBuffer,statBuffer.begin(),statBuffer.end(),m_threshold[0]);
      break;
    default:{
      std::ostringstream ess;
      ess << "method " << method << " (" << getFilterType(method) << ") not supported";
      throw(ess.str());
      break;
    }
    }
  }
  //main
  statBuffer.clear();
  for(i=dim/2;i<input.size()-dim/2;++i){
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
    switch(getFilterType(method)){
    case(filter::nvalid):
      output[i]=stat.nvalid(statBuffer);
      break;
    case(filter::median):
      output[i]=stat.median(statBuffer);
      break;
    case(filter::min):
    case(filter::erode):
      output[i]=stat.mymin(statBuffer);
      break;
    case(filter::max):
    case(filter::dilate):
      output[i]=stat.mymax(statBuffer);
      break;
    case(filter::sum):
      output[i]=sqrt(stat.sum(statBuffer));
      break;
    case(filter::var):
      output[i]=stat.var(statBuffer);
      break;
    case(filter::mean):
      output[i]=stat.mean(statBuffer);
      break;
    case(filter::percentile):
      assert(m_threshold.size());
      output[i]=stat.percentile(statBuffer,statBuffer.begin(),statBuffer.end(),m_threshold[0]);
      break;
    default:
      std::string errorString="method not supported";
      throw(errorString);
      break;
    }
    statBuffer.clear();
  }
  //end: extend input by padding
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
      T theValue=input[i-t];
      for(int iclass=0;iclass<m_class.size();++iclass){
	if(theValue==m_class[iclass]){
	  binValue=m_class[0];
	  break;
	}
      }
      if(m_class.size())
	statBuffer.push_back(binValue);
      else
	statBuffer.push_back(theValue);
      if(i+t<input.size())
	theValue=input[i+t];
      else{
	switch(getPadding(m_padding)){
	case(replicate):
	  theValue=input.back();
	  break;
	case(circular):
	  theValue=input[t-1];
	  break;
	case(zero):
	  theValue=0;
	  break;
	case(symmetric):
	default:
	  theValue=input[i-t];
	  break;
	}
      }
      for(int iclass=0;iclass<m_class.size();++iclass){
	if(theValue==m_class[iclass]){
	  binValue=m_class[0];
	  break;
	}
      }
      if(m_class.size())
	statBuffer.push_back(binValue);
      else
	statBuffer.push_back(theValue);
    }
    switch(getFilterType(method)){
    case(filter::nvalid):
      output[i]=stat.nvalid(statBuffer);
      break;
    case(filter::median):
      output[i]=stat.median(statBuffer);
      break;
    case(filter::min):
    case(filter::erode):
      output[i]=stat.mymin(statBuffer);
      break;
    case(filter::max):
    case(filter::dilate):
      output[i]=stat.mymax(statBuffer);
      break;
    case(filter::sum):
      output[i]=sqrt(stat.sum(statBuffer));
      break;
    case(filter::var):
      output[i]=stat.var(statBuffer);
      break;
    case(filter::mean):
      output[i]=stat.mean(statBuffer);
      break;
    case(filter::percentile):
      assert(m_threshold.size());
      output[i]=stat.percentile(statBuffer,statBuffer.begin(),statBuffer.end(),m_threshold[0]);
      break;
    default:
      std::string errorString="method not supported";
      throw(errorString);
      break;
    }
  }
 }

 template<class T> void Filter::smooth(T* input, int inputSize, std::vector<T>& output, short dim)
{
  assert(dim>0);
  m_taps.resize(dim);
  for(int itap=0;itap<dim;++itap)
    m_taps[itap]=1.0/dim;
  filter(input,output);
 }

template<class T> void Filter::filter(T* input, int inputSize, std::vector<T>& output)
{
  assert(inputSize>=m_taps.size());
  output.resize(inputSize);
  int i=0;

  //start: extend input by padding
  for(i=0;i<m_taps.size()/2;++i){
    //todo:introduce nodata
    output[i]=m_taps[m_taps.size()/2]*input[i];

    for(int t=1;t<=m_taps.size()/2;++t){
      output[i]+=m_taps[m_taps.size()/2+t]*input[i+t];
      if(i>=t)
	output[i]+=m_taps[m_taps.size()/2-t]*input[i-t];
      else{
	switch(getPadding(m_padding)){
	case(replicate):
	  output[i]+=m_taps[m_taps.size()/2-t]*input[0];
	  break;
	case(circular):
	  output[i]+=m_taps[m_taps.size()/2-t]*input[input.size()+i-t];
	  break;
	case(zero):
	  output[i]+=m_taps[m_taps.size()/2-t]*0;
	  break;
	case(symmetric):
	default:
	  output[i]+=m_taps[m_taps.size()/2-t]*input[t-i];
	  break;
	}
      }
    }
  }
  //main
  for(i=m_taps.size()/2;i<input.size()-m_taps.size()/2;++i){
    //todo:introduce nodata
    T leaveOut=(*(m_taps.begin()))*input[i-m_taps.size()/2];
    T include=(m_taps.back())*input[i+m_taps.size()/2];
    output[i]=0;
    for(int t=0;t<m_taps.size();++t)
      output[i]+=input[i-m_taps.size()/2+t]*m_taps[t];
  }
  //end: extend input by padding
  for(i=input.size()-m_taps.size()/2;i<input.size();++i){
    //todo:introduce nodata
    output[i]=m_taps[m_taps.size()/2]*input[i];
    //todo:introduce nodata
    for(int t=1;t<=m_taps.size()/2;++t){
      output[i]+=m_taps[m_taps.size()/2-t]*input[i-t];
      if(i+t<input.size())
	output[i]+=m_taps[m_taps.size()/2+t]*input[i+t];
      else{
	switch(getPadding(m_padding)){
	case(replicate):
	  output[i]+=m_taps[m_taps.size()/2+t]*input.back();
	  break;
	case(circular):
	  output[i]+=m_taps[m_taps.size()/2+t]*input[t-1];
	  break;
	case(zero):
	  output[i]+=m_taps[m_taps.size()/2+t]*0;
	  break;
	case(symmetric):
	default:
	  output[i]+=m_taps[m_taps.size()/2+t]*input[i-t];
	  break;
	}
      }
    }
  }
}

}

#endif /* _MYFILTER_H_ */

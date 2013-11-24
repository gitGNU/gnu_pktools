/**********************************************************************
Filter.cc: class for filtering
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
#include "Filter.h"
#include <assert.h>
#include <math.h>
#include <iostream>

filter::Filter::Filter(void)
{
}


filter::Filter::Filter(const vector<double> &taps)
  : m_taps(taps)
{
  assert(m_taps.size()%2);
}

void filter::Filter::setTaps(const vector<double> &taps)
{
  m_taps=taps;
  assert(m_taps.size()%2);
}

void filter::Filter::dwtForward(std::vector<double>& data, const std::string& wavelet_type, int family){
  //make sure data size if power of 2
  while(data.size()&(data.size()-1))
    data.push_back(data.back());
  int nsize=data.size();
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  assert(nsize);
  w=gsl_wavelet_alloc(getWaveletType(wavelet_type),family);
  work=gsl_wavelet_workspace_alloc(nsize);
  gsl_wavelet_transform_forward(w,&(data[0]),1,nsize,work);
}

void filter::Filter::dwtInverse(std::vector<double>& data, const std::string& wavelet_type, int family){
  //make sure data size if power of 2
  while(data.size()&(data.size()-1))
    data.push_back(data.back());
  int nsize=data.size();
  assert(nsize);
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  w=gsl_wavelet_alloc(getWaveletType(wavelet_type),family);
  work=gsl_wavelet_workspace_alloc(nsize);
  gsl_wavelet_transform_inverse(w,&(data[0]),1,nsize,work);
}

void filter::Filter::morphology(const ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dim, short down, int offset)
{
  Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
  Vector2d<double> lineOutput(input.nrOfBand(),input.nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int y=0;y<input.nrOfRow();++y){
    for(int iband=0;iband<input.nrOfBand();++iband)
      input.readData(lineInput[iband],GDT_Float64,y,iband);
    vector<double> pixelInput(input.nrOfBand());
    vector<double> pixelOutput(input.nrOfBand());
    for(int x=0;x<input.nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      morphology(pixelInput,pixelOutput,method,dim,down,offset);
      for(int iband=0;iband<input.nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(int iband=0;iband<input.nrOfBand();++iband){
      try{
        output.writeData(lineOutput[iband],GDT_Float64,y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::doit(const ImgReaderGdal& input, ImgWriterGdal& output, short down, int offset)
{
  Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
  Vector2d<double> lineOutput(input.nrOfBand(),input.nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int y=0;y<input.nrOfRow();++y){
    for(int iband=0;iband<input.nrOfBand();++iband)
      input.readData(lineInput[iband],GDT_Float64,y,iband);
    vector<double> pixelInput(input.nrOfBand());
    vector<double> pixelOutput(input.nrOfBand());
    for(int x=0;x<input.nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      doit(pixelInput,pixelOutput,down,offset);
      for(int iband=0;iband<input.nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(int iband=0;iband<input.nrOfBand();++iband){
      try{
        output.writeData(lineOutput[iband],GDT_Float64,y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

double filter::Filter::getCentreWavelength(const std::vector<double> &wavelengthIn, const Vector2d<double>& srf, const std::string& interpolationType, double delta, bool verbose)
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

  std::vector<double> wavelengthOut(wavelength_fine.size());

  for(int iband=0;iband<wavelengthOut.size();++iband)
    wavelengthOut[iband]=wavelength_fine[iband]*srf_fine[iband];

  stat.initSpline(splineOut,&(wavelength_fine[0]),&(wavelengthOut[0]),wavelength_fine.size());
  double centreWavelength=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
  
  gsl_spline_free(splineOut);
  gsl_interp_accel_free(accOut);

  return(centreWavelength);
}

// void filter::Filter::applyFwhm(const vector<double> &wavelengthIn, const ImgReaderGdal& input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, ImgWriterGdal& output, bool verbose){
//   Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
//   Vector2d<double> lineOutput(wavelengthOut.size(),input.nrOfCol());
//   const char* pszMessage;
//   void* pProgressArg=NULL;
//   GDALProgressFunc pfnProgress=GDALTermProgress;
//   double progress=0;
//   pfnProgress(progress,pszMessage,pProgressArg);
//   for(int y=0;y<input.nrOfRow();++y){
//     for(int iband=0;iband<input.nrOfBand();++iband)
//       input.readData(lineInput[iband],GDT_Float64,y,iband);
//     applyFwhm<double>(wavelengthIn,lineInput,wavelengthOut,fwhm, interpolationType, lineOutput, verbose);
//     for(int iband=0;iband<output.nrOfBand();++iband){
//       try{
//         output.writeData(lineOutput[iband],GDT_Float64,y,iband);
//       }
//       catch(string errorstring){
//         cerr << errorstring << "in band " << iband << ", line " << y << endl;
//       }
//     }
//     progress=(1.0+y)/output.nrOfRow();
//     pfnProgress(progress,pszMessage,pProgressArg);
//   }
// }

// void filter::Filter::applySrf(const vector<double> &wavelengthIn, const ImgReaderGdal& input, const vector< Vector2d<double> > &srf, const std::string& interpolationType, ImgWriterGdal& output, bool verbose){
//   assert(output.nrOfBand()==srf.size());
//   double centreWavelength=0;
//   Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
//   const char* pszMessage;
//   void* pProgressArg=NULL;
//   GDALProgressFunc pfnProgress=GDALTermProgress;
//   double progress=0;
//   pfnProgress(progress,pszMessage,pProgressArg);
//   for(int y=0;y<input.nrOfRow();++y){
//     for(int iband=0;iband<input.nrOfBand();++iband)
//       input.readData(lineInput[iband],GDT_Float64,y,iband);
//     for(int isrf=0;isrf<srf.size();++isrf){
//       vector<double> lineOutput(input.nrOfCol());
//       centreWavelength=applySrf<double>(wavelengthIn,lineInput,srf[isrf], interpolationType, lineOutput, verbose);
//       for(int iband=0;iband<output.nrOfBand();++iband){
//         try{
//           output.writeData(lineOutput,GDT_Float64,y,isrf);
//         }
//         catch(string errorstring){
//           cerr << errorstring << "in band " << iband << ", line " << y << endl;
//         }
//       }
//     }
//     progress=(1.0+y)/output.nrOfRow();
//     pfnProgress(progress,pszMessage,pProgressArg);
//   }
// }

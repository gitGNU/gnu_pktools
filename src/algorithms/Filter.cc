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

using namespace std;

filter::Filter::Filter(void)
  : m_padding("symmetric")
{
}


filter::Filter::Filter(const vector<double> &taps)
  : m_padding("symmetric")
{
  setTaps(taps);
}

void filter::Filter::setTaps(const vector<double> &taps, bool normalize)
{
  m_taps.resize(taps.size());
  double norm=0;
  for(int itap=0;itap<taps.size();++itap)
    norm+=taps[itap];
  if(norm){
    for(int itap=0;itap<taps.size();++itap)
      m_taps[itap]=taps[itap]/norm;
  }
  else
    m_taps=taps;
  assert(m_taps.size()%2);
}

unsigned int filter::Filter::pushNoDataValue(double noDataValue){
  if(find(m_noDataValues.begin(),m_noDataValues.end(),noDataValue)==m_noDataValues.end())
    m_noDataValues.push_back(noDataValue);
  return m_noDataValues.size();
};

unsigned int filter::Filter::setNoDataValues(std::vector<double> vnodata){
  m_noDataValues=vnodata;
  return m_noDataValues.size();
};

void filter::Filter::dwtForward(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& wavelet_type, int family){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      dwtForward(pixelInput,wavelet_type,family);
      for(int iband=0;iband<input->nrOfBand();++iband)
        lineOutput[iband][x]=pixelInput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::dwtInverse(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& wavelet_type, int family){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      dwtInverse(pixelInput,wavelet_type,family);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband)
        lineOutput[iband][x]=pixelInput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::dwtCut(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& wavelet_type, int family, double cut){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      dwtCut(pixelInput,wavelet_type,family,cut);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband)
        lineOutput[iband][x]=pixelInput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::dwtCutFrom(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& wavelet_type, int family, int band){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      dwtForward(pixelInput,wavelet_type,family);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband){
	if(iband>=band)
	  pixelInput[iband]=0;
      }
      dwtInverse(pixelInput,wavelet_type,family);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband)
	lineOutput[iband][x]=pixelInput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

//todo: support different padding strategies
void filter::Filter::dwtForward(std::vector<double>& data, const std::string& wavelet_type, int family){
  int origsize=data.size();
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
  data.erase(data.begin()+origsize,data.end());
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
}

//todo: support different padding strategies
void filter::Filter::dwtInverse(std::vector<double>& data, const std::string& wavelet_type, int family){
  int origsize=data.size();
  //make sure data size if power of 2
  while(data.size()&(data.size()-1))
    data.push_back(data.back());
  int nsize=data.size();
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  assert(nsize);
  w=gsl_wavelet_alloc(getWaveletType(wavelet_type),family);
  work=gsl_wavelet_workspace_alloc(nsize);
  gsl_wavelet_transform_inverse(w,&(data[0]),1,nsize,work);
  data.erase(data.begin()+origsize,data.end());
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
}

//todo: support different padding strategies
void filter::Filter::dwtCut(std::vector<double>& data, const std::string& wavelet_type, int family, double cut){
  int origsize=data.size();
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
  std::vector<double> abscoeff(data.size());
  size_t* p=new size_t[data.size()];
  for(int index=0;index<data.size();++index){
    abscoeff[index]=fabs(data[index]);
  }
  int nc=(100-cut)/100.0*nsize;
  gsl_sort_index(p,&(abscoeff[0]),1,nsize);
  for(int i=0;(i+nc)<nsize;i++)
    data[p[i]]=0;
  gsl_wavelet_transform_inverse(w,&(data[0]),1,nsize,work);
  data.erase(data.begin()+origsize,data.end());
  delete[] p;
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
}

void filter::Filter::morphology(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& method, int dim, short verbose)
{
  // bool bverbose=(verbose>1)? true:false;
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  if(verbose)
    std::cout << "Number of bands in input: " << lineInput.size() << std::endl;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<input->nrOfRow();++y){
    try{
      for(unsigned int iband=0;iband<input->nrOfBand();++iband){
        input->readData(lineInput[iband],y,iband);
      }
    }
    catch(string errorString){
      std::cerr << "Error: could not read data from input" << endl;
      exit(1);
    }
    catch(...){
      std::cerr << "Error: " << std::endl;
    }
    vector<double> pixelInput(input->nrOfBand());
    vector<double> pixelOutput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      filter(pixelInput,pixelOutput,method,dim);
      // morphology(pixelInput,pixelOutput,method,dim,bverbose);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::smoothNoData(std::shared_ptr<ImgRaster> input, const std::string& interpolationType, std::shared_ptr<ImgRaster> output)
{
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    vector<double> pixelOutput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      smoothNoData(pixelInput,interpolationType,pixelOutput);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::smooth(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, short dim)
{
  assert(dim>0);
  m_taps.resize(dim);
  for(int itap=0;itap<dim;++itap)
    m_taps[itap]=1.0/dim;
  filter(input,output);
}

void filter::Filter::filter(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output)
{
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    vector<double> pixelOutput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      filter(pixelInput,pixelOutput);
      for(unsigned int iband=0;iband<input->nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::stat(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& method)
{
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  assert(output->nrOfCol()==input->nrOfCol());
  vector<double> lineOutput(output->nrOfCol());
  statfactory::StatFactory stat;
  stat.setNoDataValues(m_noDataValues);
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      switch(getFilterType(method)){
      case(filter::median):
	lineOutput[x]=stat.median(pixelInput);
	break;
      case(filter::min):
	lineOutput[x]=stat.mymin(pixelInput);
	break;
      case(filter::max):
	lineOutput[x]=stat.mymax(pixelInput);
	break;
      case(filter::sum):
	lineOutput[x]=stat.sum(pixelInput);
	break;
      case(filter::var):
	lineOutput[x]=stat.var(pixelInput);
	break;
      case(filter::stdev):
	lineOutput[x]=sqrt(stat.var(pixelInput));
	break;
      case(filter::mean):
	lineOutput[x]=stat.mean(pixelInput);
	break;
      case(filter::percentile):
	assert(m_threshold.size());
	lineOutput[x]=stat.percentile(pixelInput,pixelInput.begin(),pixelInput.end(),m_threshold[0]);
	break;
      default:
	std::string errorString="method not supported";
	throw(errorString);
	break;
      }
    }
    try{
      output->writeData(lineOutput,y);
    }
    catch(string errorstring){
      cerr << errorstring << "in line " << y << endl;
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::stats(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const vector<std::string>& methods)
{
  assert(output->nrOfBand()==methods.size());
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  assert(output->nrOfCol()==input->nrOfCol());
  Vector2d<double> lineOutput(methods.size(),output->nrOfCol());
  statfactory::StatFactory stat;
  stat.setNoDataValues(m_noDataValues);
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      int ithreshold=0;//threshold to use for percentiles
      for(int imethod=0;imethod<methods.size();++imethod){
	switch(getFilterType(methods[imethod])){
	case(filter::nvalid):
	  lineOutput[imethod][x]=stat.nvalid(pixelInput);
	  break;
	case(filter::median):
	  lineOutput[imethod][x]=stat.median(pixelInput);
	  break;
	case(filter::min):
	  lineOutput[imethod][x]=stat.mymin(pixelInput);
	  break;
	case(filter::max):
	  lineOutput[imethod][x]=stat.mymax(pixelInput);
	  break;
	case(filter::sum):
	  lineOutput[imethod][x]=stat.sum(pixelInput);
	  break;
	case(filter::var):
	  lineOutput[imethod][x]=stat.var(pixelInput);
	  break;
	case(filter::stdev):
	  lineOutput[imethod][x]=sqrt(stat.var(pixelInput));
	  break;
	case(filter::mean):
	  lineOutput[imethod][x]=stat.mean(pixelInput);
	  break;
	case(filter::percentile):{
	  assert(m_threshold.size());
	  double threshold=(ithreshold<m_threshold.size())? m_threshold[ithreshold] : m_threshold[0];
	  lineOutput[imethod][x]=stat.percentile(pixelInput,pixelInput.begin(),pixelInput.end(),threshold);
	  ++ithreshold;
	  break;
	}
	default:
	  std::string errorString="method not supported";
	  throw(errorString);
	  break;
	}
      }
    }
    for(int imethod=0;imethod<methods.size();++imethod){
      try{
	output->writeData(lineOutput[imethod],y,imethod);
      }
      catch(string errorstring){
	cerr << errorstring << "in line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::filter(std::shared_ptr<ImgRaster> input, std::shared_ptr<ImgRaster> output, const std::string& method, int dim)
{
  Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
  Vector2d<double> lineOutput(input->nrOfBand(),input->nrOfCol());;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int y=0;y<input->nrOfRow();++y){
    for(unsigned int iband=0;iband<input->nrOfBand();++iband)
      input->readData(lineInput[iband],y,iband);
    vector<double> pixelInput(input->nrOfBand());
    vector<double> pixelOutput;
    for(unsigned int x=0;x<input->nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      filter(pixelInput,pixelOutput,method,dim);
      for(unsigned int iband=0;iband<pixelOutput.size();++iband){
        lineOutput[iband][x]=pixelOutput[iband];
	// if(pixelInput[iband]!=0)
	//   assert(pixelOutput[iband]!=0);
      }
    }
    for(unsigned int iband=0;iband<input->nrOfBand();++iband){
      try{
        output->writeData(lineOutput[iband],y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::getSavGolayCoefficients(vector<double> &tapz, int np, int nl, int nr, int ld, int m) {
  int j, k, imj, ipj, kk, mm;
  double d, fac, sum;

  // c.resize(nl+1+nr);
  vector<double> tmpc(np);
  if(np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m) {
    cerr << "bad args in savgol" << endl;
    return;
  }
  vector<int> indx(m + 1, 0);
  vector<double> a((m + 1) * (m + 1), 0.0);
  vector<double> b(m + 1, 0.0);

  for(ipj = 0; ipj <= (m << 1); ++ipj) {
    sum = (ipj ? 0.0 : 1.0);
    for(k = 1; k <= nr; ++k)
      sum += (int)pow((double)k, (double)ipj);
    for(k = 1; k <= nl; ++k)
      sum += (int)pow((double) - k, (double)ipj);
    mm = (ipj < 2 * m - ipj ? ipj : 2 * m - ipj);
    for(imj = -mm; imj <= mm; imj += 2)
      a[(ipj + imj) / 2 * (m + 1) + (ipj - imj) / 2] = sum;
  }
  ludcmp(a, indx, d);

  for(j = 0; j < m + 1; ++j)
    b[j] = 0.0;
  b[ld] = 1.0;

  lubksb(a, indx, b);
  // for(kk = 0; kk < np; ++kk)
  //   c[kk] = 0.0;
  for(k = -nl; k <= nr; ++k) {
  // for(k = -nl; k < nr; ++k) {
    sum = b[0];
    fac = 1.0;
    for(mm = 1; mm <= m; ++mm)
      sum += b[mm] * (fac *= k);
    // store in wrap=around order
    kk = (np - k) % np;
    //re-order c as I need for taps
    // kk=k+nl;
    tmpc[kk] = sum;
  }
  tapz.resize(nl+1+nr);
  //  for(k=0;k<nl+1+nr)
  tapz[tapz.size()/2]=tmpc[0];
  //past data points
  for(k=1;k<=tapz.size()/2;++k)
    tapz[tapz.size()/2-k]=tmpc[k];
  //future data points
  for(k=1;k<=tapz.size()/2;++k)
    tapz[tapz.size()/2+k]=tmpc[np-k];
}

void filter::Filter::ludcmp(vector<double> &a, vector<int> &indx, double &d) {
  const double TINY = 1.0e-20;
  int i, imax = -1, j, k;
  double big, dum, sum, temp;

  int n = indx.size();
  vector<double> vv(n, 0.0);

  d = 1.0;
  for(i = 0; i < n; ++i) {
    big = 0.0;
    for(j = 0; j < n; ++j)
      if((temp = fabs(a[i * n + j])) > big) big = temp;

    if(big == 0.0) {
      cerr << "Singular matrix in routine ludcmp" << endl;
      return;
    }
    vv[i] = 1. / big;
  }

  for(j = 0; j < n; ++j) {
    for(i = 0; i < j; ++i) {
      sum = a[i * n + j];
      for(k = 0; k < i; ++k)
	sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
    }
    big = 0.0;
    for(i = j; i < n; ++i) {
      sum = a[i * n + j];
      for(k = 0; k < j; ++k)
	sum -= a[i * n + k] * a[k * n + j];
      a[i * n + j] = sum;
      if((dum = vv[i] * fabs(sum)) >= big) {
	big = dum;
	imax = i;
      }
    }

    if(j != imax) {
      for(k = 0; k < n; ++k) {
	dum = a[imax * n + k];
	a[imax * n + k] = a[j * n + k];
	a[j * n + k] = dum;
      }
      d = -d;
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if(a[j * n + j] == 0.0) a[j * n + j] = TINY;
    if(j != n - 1) {
      dum = 1. / a[j * n + j];
      for(i = j + 1; i < n; ++i)
	a[i * n + j] *= dum;
    }
  }
}

void filter::Filter::lubksb(vector<double> &a, vector<int> &indx, vector<double> &b) {
  int i, ii = 0, ip, j;
  double sum;
  int n = indx.size();

  for(i = 0; i < n; ++i) {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if(ii != 0)
      for(j = ii - 1; j < i; ++j)
	sum -= a[i * n + j] * b[j];
    else if(sum != 0.0)
      ii = i + 1;
    b[i] = sum;
  }
  for(i = n - 1; i >= 0; --i) {
    sum = b[i];
    for(j = i + 1; j < n; ++j)
      sum -= a[i * n + j] * b[j];
    b[i] = sum / a[i * n + i];
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

  for(unsigned int iband=0;iband<wavelengthOut.size();++iband)
    wavelengthOut[iband]=wavelength_fine[iband]*srf_fine[iband];

  stat.initSpline(splineOut,&(wavelength_fine[0]),&(wavelengthOut[0]),wavelength_fine.size());
  double centreWavelength=gsl_spline_eval_integ(splineOut,start,end,accOut)/norm;
  
  gsl_spline_free(splineOut);
  gsl_interp_accel_free(accOut);

  return(centreWavelength);
}

// void filter::Filter::applyFwhm(const vector<double> &wavelengthIn, const std::shared_ptr<ImgRaster> input, const vector<double> &wavelengthOut, const vector<double> &fwhm, const std::string& interpolationType, std::shared_ptr<ImgRaster> output, bool verbose){
//   Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
//   Vector2d<double> lineOutput(wavelengthOut.size(),input->nrOfCol());
//   const char* pszMessage;
//   void* pProgressArg=NULL;
//   GDALProgressFunc pfnProgress=GDALTermProgress;
//   double progress=0;
//   pfnProgress(progress,pszMessage,pProgressArg);
//   for(unsigned int y=0;y<input->nrOfRow();++y){
//     for(unsigned int iband=0;iband<input->nrOfBand();++iband)
//       input->readData(lineInput[iband],GDT_Float64,y,iband);
//     applyFwhm<double>(wavelengthIn,lineInput,wavelengthOut,fwhm, interpolationType, lineOutput, verbose);
//     for(unsigned int iband=0;iband<output->nrOfBand();++iband){
//       try{
//         output->writeData(lineOutput[iband],GDT_Float64,y,iband);
//       }
//       catch(string errorstring){
//         cerr << errorstring << "in band " << iband << ", line " << y << endl;
//       }
//     }
//     progress=(1.0+y)/output->nrOfRow();
//     pfnProgress(progress,pszMessage,pProgressArg);
//   }
// }

// void filter::Filter::applySrf(const vector<double> &wavelengthIn, const std::shared_ptr<ImgRaster> input, const vector< Vector2d<double> > &srf, const std::string& interpolationType, std::shared_ptr<ImgRaster> output, bool verbose){
//   assert(output->nrOfBand()==srf.size());
//   double centreWavelength=0;
//   Vector2d<double> lineInput(input->nrOfBand(),input->nrOfCol());
//   const char* pszMessage;
//   void* pProgressArg=NULL;
//   GDALProgressFunc pfnProgress=GDALTermProgress;
//   double progress=0;
//   pfnProgress(progress,pszMessage,pProgressArg);
//   for(unsigned int y=0;y<input->nrOfRow();++y){
//     for(unsigned int iband=0;iband<input->nrOfBand();++iband)
//       input->readData(lineInput[iband],GDT_Float64,y,iband);
//     for(int isrf=0;isrf<srf.size();++isrf){
//       vector<double> lineOutput(input->nrOfCol());
//       centreWavelength=applySrf<double>(wavelengthIn,lineInput,srf[isrf], interpolationType, lineOutput, verbose);
//       for(unsigned int iband=0;iband<output->nrOfBand();++iband){
//         try{
//           output->writeData(lineOutput,GDT_Float64,y,isrf);
//         }
//         catch(string errorstring){
//           cerr << errorstring << "in band " << iband << ", line " << y << endl;
//         }
//       }
//     }
//     progress=(1.0+y)/output->nrOfRow();
//     pfnProgress(progress,pszMessage,pProgressArg);
//   }
// }

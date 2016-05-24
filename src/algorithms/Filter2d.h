/**********************************************************************
Filter2d.h: class for filtering images
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
#ifndef _MYFILTER2D_H_
#define _MYFILTER2D_H_

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

#ifndef DEG2RAD
#define DEG2RAD(DEG) (DEG/180.0*PI)
#endif

#ifndef RAD2DEG
#define RAD2DEG(RAD) (RAD/PI*180)
#endif

#ifdef WIN32
#include <process.h>
#define getpid _getpid
#endif

#include <assert.h>
#include <math.h>
#include <limits>
#include <vector>
#include <string>
#include <map>
extern "C" {
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>
#include <gsl/gsl_wavelet2d.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}
#include "base/Vector2d.h"
#include "Filter.h"
#include "imageclasses/ImgRasterGdal.h"
#include "algorithms/StatFactory.h"

namespace filter2d
{
  enum FILTER_TYPE { median=100, var=101 , min=102, max=103, sum=104, mean=105, minmax=106, dilate=107, erode=108, close=109, open=110, homog=111, sobelx=112, sobely=113, sobelxy=114, sobelyx=115, smooth=116, density=117, mode=118, mixed=119, threshold=120, ismin=121, ismax=122, heterog=123, order=124, stdev=125, mrf=126, dwt=127, dwti=128, dwt_cut=129, scramble=130, shift=131, linearfeature=132, smoothnodata=133, countid=134, dwt_cut_from=135, savgolay=136, percentile=137, proportion=138, nvalid=139, sauvola=140};
  
class Filter2d
{
public:
  enum RESAMPLE { NEAR = 0, BILINEAR = 1, BICUBIC = 2 };//bicubic not supported yet...
  Filter2d(void);
  Filter2d(const Vector2d<double> &taps);
  virtual ~Filter2d(){};
  static FILTER_TYPE getFilterType(const std::string filterType){
    std::map<std::string, FILTER_TYPE> m_filterMap;
    initMap(m_filterMap);
    return m_filterMap[filterType];
  };
  static const RESAMPLE getResampleType(const std::string resampleType){
    if(resampleType=="near") return(NEAR);
    else if(resampleType=="bilinear") return(BILINEAR);
    else{
      std::string errorString="resampling type not supported: ";
      errorString+=resampleType;
      errorString+=" use near or bilinear";
      throw(errorString);
    }
  };

  void setTaps(const Vector2d<double> &taps);
  /* void setNoValue(double noValue=0){m_noValue=noValue;}; */
  void pushClass(short theClass=1){m_class.push_back(theClass);};
  int pushNoDataValue(double noDataValue=0);//{m_mask.push_back(theMask);};
  void pushThreshold(double theThreshold){m_threshold.push_back(theThreshold);};
  void setThresholds(const std::vector<double>& theThresholds){m_threshold=theThresholds;};
  void setClasses(const std::vector<short>& theClasses){m_class=theClasses;};
  void filter(ImgRasterGdal& input, ImgRasterGdal& output, bool absolute=false, bool normalize=false, bool noData=false);
  void smooth(ImgRasterGdal& input, ImgRasterGdal& output,int dim);
  void smooth(ImgRasterGdal& input, ImgRasterGdal& output,int dimX, int dimY);
  void smoothNoData(ImgRasterGdal& input, ImgRasterGdal& output,int dim);
  void smoothNoData(ImgRasterGdal& input, ImgRasterGdal& output,int dimX, int dimY);
  template<class T1, class T2> void filter(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector);
  template<class T1, class T2> void smooth(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector,int dim);
  template<class T1, class T2> void smooth(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector,int dimX, int dimY);
  void dwtForward(ImgRasterGdal& input, ImgRasterGdal& output, const std::string& wavelet_type, int family);
  void dwtInverse(ImgRasterGdal& input, ImgRasterGdal& output, const std::string& wavelet_type, int family);
  void dwtCut(ImgRasterGdal& input, ImgRasterGdal& output, const std::string& wavelet_type, int family, double cut, bool verbose=false);
  template<class T> void dwtForward(Vector2d<T>& data, const std::string& wavelet_type, int family);
  template<class T> void dwtInverse(Vector2d<T>& data, const std::string& wavelet_type, int family);
  template<class T> void dwtCut(Vector2d<T>& data, const std::string& wavelet_type, int family, double cut);
  void majorVoting(ImgRasterGdal& input, ImgRasterGdal& output, int dim=0,const std::vector<int> &prior=std::vector<int>());
  /* void homogeneousSpatial(const std::string& inputFilename, const std::string& outputFilename, int dim, bool disc=false, int noValue=0); */
  void doit(ImgRasterGdal& input, ImgRasterGdal& output, const std::string& method, int dim, short down=1, bool disc=false);
  void doit(ImgRasterGdal& input, ImgRasterGdal& output, const std::string& method, int dimX, int dimY, short down=1, bool disc=false);
  void mrf(ImgRasterGdal& input, ImgRasterGdal& output, int dimX, int dimY, double beta, bool eightConnectivity=true, short down=1, bool verbose=false);
  void mrf(ImgRasterGdal& input, ImgRasterGdal& output, int dimX, int dimY, Vector2d<double> beta, bool eightConnectivity=true, short down=1, bool verbose=false);
  template<class T1, class T2> void doit(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector, const std::string& method, int dimX, int dimY, short down=1, bool disc=false);
  void median(ImgRasterGdal& input, ImgRasterGdal& output, int dim, bool disc=false);
  void var(ImgRasterGdal& input, ImgRasterGdal& output, int dim, bool disc=false);
  void morphology(ImgRasterGdal& input, ImgRasterGdal& output, const std::string& method, int dimX, int dimY, const std::vector<double> &angle=std::vector<double>(), bool disc=false);
  template<class T> unsigned long int morphology(const Vector2d<T>& input, Vector2d<T>& output, const std::string& method, int dimX, int dimY, bool disc=false, double hThreshold=0);
  template<class T> unsigned long int dsm2dtm_nwse(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim=3);
  template<class T> unsigned long int dsm2dtm_nesw(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim=3);
  template<class T> unsigned long int dsm2dtm_senw(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim=3);
  template<class T> unsigned long int dsm2dtm_swne(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim=3);
  template<class T> void shadowDsm(const Vector2d<T>& input, Vector2d<T>& output, double sza, double saa, double pixelSize, short shadowFlag=1);
  void shadowDsm(ImgRasterGdal& input, ImgRasterGdal& output, double sza, double saa, double pixelSize, short shadowFlag=1);
  //  void dwt_texture(ImgRasterGdal& input, ImgRasterGdal& output, int dim, int scale, int down=1, int iband=0, bool verbose=false);
  void shift(ImgRasterGdal& input, ImgRasterGdal& output, double offsetX=0, double offsetY=0, double randomSigma=0, RESAMPLE resample=BILINEAR, bool verbose=false);
  template<class T> void shift(const Vector2d<T>& input, Vector2d<T>& output, double offsetX=0, double offsetY=0, double randomSigma=0, RESAMPLE resample=NEAR, bool verbose=false);
  void linearFeature(const Vector2d<float>& input, std::vector< Vector2d<float> >& output, float angle=361, float angleStep=1, float maxDistance=0, float eps=0, bool l1=true, bool a1=true, bool l2=true, bool a2=true, bool verbose=false);
  void linearFeature(ImgRasterGdal& input, ImgRasterGdal& output, float angle=361, float angleStep=1, float maxDistance=0, float eps=0, bool l1=true, bool a1=true, bool l2=true, bool a2=true, int band=0, bool verbose=false);
  
private:
  static void initMap(std::map<std::string, FILTER_TYPE>& m_filterMap){
    //initialize selMap
    m_filterMap["median"]=filter2d::median;
    m_filterMap["nvalid"]=filter2d::nvalid;
    m_filterMap["var"]=filter2d::var;
    m_filterMap["min"]=filter2d::min;
    m_filterMap["max"]=filter2d::max;
    m_filterMap["sum"]=filter2d::sum;
    m_filterMap["mean"]=filter2d::mean;
    m_filterMap["minmax"]=filter2d::minmax;
    m_filterMap["dilate"]=filter2d::dilate;
    m_filterMap["erode"]=filter2d::erode;
    m_filterMap["close"]=filter2d::close;
    m_filterMap["open"]=filter2d::open;
    m_filterMap["homog"]=filter2d::homog;
    m_filterMap["sobelx"]=filter2d::sobelx;
    m_filterMap["sobely"]=filter2d::sobely;
    m_filterMap["sobelxy"]=filter2d::sobelxy;
    m_filterMap["sobelyx"]=filter2d::sobelyx;
    m_filterMap["smooth"]=filter2d::smooth;
    m_filterMap["density"]=filter2d::density;
    m_filterMap["mode"]=filter2d::mode;
    m_filterMap["mixed"]=filter2d::mixed;
    m_filterMap["smoothnodata"]=filter2d::smoothnodata;
    m_filterMap["threshold"]=filter2d::threshold;
    m_filterMap["ismin"]=filter2d::ismin;
    m_filterMap["ismax"]=filter2d::ismax;
    m_filterMap["heterog"]=filter2d::heterog;
    m_filterMap["sauvola"]=filter2d::sauvola;
    m_filterMap["order"]=filter2d::order;
    m_filterMap["stdev"]=filter2d::stdev;
    m_filterMap["mrf"]=filter2d::mrf;
    m_filterMap["dwt"]=filter2d::dwt;
    m_filterMap["dwti"]=filter2d::dwti;
    m_filterMap["dwt_cut"]=filter2d::dwt_cut;
    m_filterMap["dwt_cut_from"]=filter2d::dwt_cut_from;
    m_filterMap["scramble"]=filter2d::scramble;
    m_filterMap["shift"]=filter2d::shift;
    m_filterMap["linearfeature"]=filter2d::linearfeature;
    m_filterMap["countid"]=filter2d::countid;
    m_filterMap["savgolay"]=filter2d::savgolay;
    m_filterMap["percentile"]=filter2d::percentile;
    m_filterMap["proportion"]=filter2d::proportion;
  }

  Vector2d<double> m_taps;
  /* double m_noValue; */
  std::vector<short> m_class;
  /* std::vector<short> m_mask; */
  std::vector<double> m_noDataValues;
  std::vector<double> m_threshold;
};


 template<class T1, class T2> void Filter2d::smooth(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector,int dim)
  {
    smooth(inputVector,outputVector,dim,dim);
  }

 template<class T1, class T2> void Filter2d::smooth(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector,int dimX, int dimY)
  {
    m_taps.resize(dimY);
    for(int j=0;j<dimY;++j){
      m_taps[j].resize(dimX);
      for(int i=0;i<dimX;++i)
	m_taps[j][i]=1.0/dimX/dimY;
    }
    filter(inputVector,outputVector);
  }
  
  template<class T1, class T2> void Filter2d::filter(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector)
  {
    outputVector.resize(inputVector.size());
    int dimX=m_taps[0].size();//horizontal!!!
    int dimY=m_taps.size();//vertical!!!
    Vector2d<T1> inBuffer(dimY);
    std::vector<T2> outBuffer(inputVector[0].size());
    //initialize last half of inBuffer
    int indexI=0;
    int indexJ=0;
    //initialize last half of inBuffer
    for(int j=-(dimY-1)/2;j<=dimY/2;++j){
      inBuffer[indexJ]=inputVector[abs(j)];
      ++indexJ;
    }

    for(int y=0;y<inputVector.size();++y){
      if(y){//inBuffer already initialized for y=0
        //erase first line from inBuffer
        inBuffer.erase(inBuffer.begin());
        //read extra line and push back to inBuffer if not out of bounds
        if(y+dimY/2<inputVector.size()){
	  //allocate buffer
          inBuffer.push_back(inputVector[y+dimY/2]);
        }
        else{
          int over=y+dimY/2-inputVector.nRows();
          int index=(inBuffer.size()-1)-over;
          assert(index>=0);
          assert(index<inBuffer.size());
          inBuffer.push_back(inBuffer[index]);
        }
      }
      for(int x=0;x<inputVector.nCols();++x){
        outBuffer[x]=0;
	for(int j=-(dimY-1)/2;j<=dimY/2;++j){
	  for(int i=-(dimX-1)/2;i<=dimX/2;++i){
	    indexI=x+i;
	    indexJ=(dimY-1)/2+j;
	    //check if out of bounds
	    if(x<(dimX-1)/2)
	      indexI=x+abs(i);
	    else if(x>=inputVector.nCols()-(dimX-1)/2)
	      indexI=x-abs(i);
	    if(y<(dimY-1)/2)
	      indexJ=(dimY-1)/2+abs(j);
	    else if(y>=inputVector.nRows()-(dimY-1)/2)
	      indexJ=(dimY-1)/2-abs(j);
            outBuffer[x]+=(m_taps[(dimY-1)/2+j][(dimX-1)/2+i]*inBuffer[indexJ][indexI]);
          }
        }
      }
      //copy outBuffer to outputVector
      outputVector[y]=outBuffer;
    }
  }

template<class T1, class T2> void Filter2d::doit(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector, const std::string& method, int dimX, int dimY, short down, bool disc)
{
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  statfactory::StatFactory stat;
  double noDataValue=0;
  if(m_noDataValues.size()){
    stat.setNoDataValues(m_noDataValues);
    noDataValue=m_noDataValues[0];
  }
  assert(dimX);
  assert(dimY);

  outputVector.resize((inputVector.size()+down-1)/down);
  Vector2d<T1> inBuffer(dimY);
  std::vector<T2> outBuffer((inputVector[0].size()+down-1)/down);
  
  int indexI=0;
  int indexJ=0;
  //initialize last half of inBuffer
  for(int j=-(dimY-1)/2;j<=dimY/2;++j){
    inBuffer[indexJ]=inputVector[abs(j)];
    ++indexJ;
  }
  for(int y=0;y<inputVector.size();++y){
    
    if(y){//inBuffer already initialized for y=0
      //erase first line from inBuffer
      inBuffer.erase(inBuffer.begin());
      //read extra line and push back to inBuffer if not out of bounds
      if(y+dimY/2<inputVector.size())
	inBuffer.push_back(inputVector[y+dimY/2]);
      else{
        int over=y+dimY/2-inputVector.size();
        int index=(inBuffer.size()-1)-over;
        assert(index>=0);
        assert(index<inBuffer.size());
        inBuffer.push_back(inBuffer[index]);
      }
    }
    if((y+1+down/2)%down)
      continue;
    for(int x=0;x<inputVector[0].size();++x){
      if((x+1+down/2)%down)
        continue;
      outBuffer[x/down]=0;

      std::vector<T1> windowBuffer;
      std::map<int,int> occurrence;
      int centre=dimX*(dimY-1)/2+(dimX-1)/2;
      bool centreMasked=false;
      for(int j=-(dimY-1)/2;j<=dimY/2;++j){
        for(int i=-(dimX-1)/2;i<=dimX/2;++i){
	  indexI=x+i;
	  //check if out of bounds
          if(indexI<0)
            indexI=-indexI;
          else if(indexI>=inputVector[0].size())
            indexI=inputVector[0].size()-i;
          if(y+j<0)
            indexJ=-j;
          else if(y+j>=inputVector.size())
            indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
          else
            indexJ=(dimY-1)/2+j;
          windowBuffer.push_back(inBuffer[indexJ][indexI]);
          if(!stat.isNoData(inBuffer[indexJ][indexI])){
            std::vector<short>::const_iterator vit=m_class.begin();
            //todo: test if this works (only add occurrence if within defined classes)!
            if(!m_class.size())
              ++occurrence[inBuffer[indexJ][indexI]];
            else{
              while(vit!=m_class.end()){
                if(inBuffer[indexJ][indexI]==*(vit++))
                  ++occurrence[inBuffer[indexJ][indexI]];
              }
            }
          }
        }
      }
      switch(getFilterType(method)){
      case(filter2d::nvalid):
	outBuffer[x/down]=stat.nvalid(windowBuffer);
        break;
      case(filter2d::median):
        outBuffer[x/down]=stat.median(windowBuffer);
        break;
      case(filter2d::var):{
        outBuffer[x/down]=stat.var(windowBuffer);
        break;
      }
      case(filter2d::stdev):{
        T2 varValue=stat.var(windowBuffer);
        if(stat.isNoData(varValue))
          outBuffer[x/down]=noDataValue;
        else          
          outBuffer[x/down]=sqrt(varValue);
        break;
      }
      case(filter2d::mean):{
        if(windowBuffer.empty())
          outBuffer[x/down]=noDataValue;
        else
          outBuffer[x/down]=stat.mean(windowBuffer);
        break;
      }
      case(filter2d::min):{
        outBuffer[x/down]=stat.mymin(windowBuffer);
        break;
      }
      case(filter2d::ismin):{
        T1 minValue=stat.mymin(windowBuffer);
        if(stat.isNoData(minValue))
          outBuffer[x/down]=noDataValue;
        else
          outBuffer[x/down]=(windowBuffer[centre]==minValue)? 1:0;
        break;
      }
      case(filter2d::minmax):{
        T1 min=0;
        T1 max=0;
        stat.minmax(windowBuffer,windowBuffer.begin(),windowBuffer.end(),min,max);
        if(min!=max)
          outBuffer[x/down]=0;
        else
          outBuffer[x/down]=windowBuffer[centre];//centre pixels
        break;
      }
      case(filter2d::max):{
        outBuffer[x/down]=stat.mymax(windowBuffer);
        break;
      }
      case(filter2d::ismax):{
        T1 maxValue=stat.mymax(windowBuffer);
        if(stat.isNoData(maxValue))
          outBuffer[x/down]=noDataValue;
        else
          outBuffer[x/down]=(windowBuffer[centre]==maxValue)? 1:0;
        break;
      }
      case(filter2d::order):{
        stat.eraseNoData(windowBuffer);
        if(windowBuffer.empty())
          outBuffer[x/down]=noDataValue;
        else{
          double lbound=0;
          double ubound=dimX*dimY;
          double theMin=stat.mymin(windowBuffer);
          double theMax=stat.mymax(windowBuffer);
          double scale=(ubound-lbound)/(theMax-theMin);
          outBuffer[x/down]=static_cast<short>(scale*(windowBuffer[centre]-theMin)+lbound);
        }
        break;
      }
      case(filter2d::sum):{
        outBuffer[x/down]=stat.sum(windowBuffer);
        break;
      }
      case(filter2d::percentile):{
	assert(m_threshold.size());
        outBuffer[x/down]=stat.percentile(windowBuffer,windowBuffer.begin(),windowBuffer.end(),m_threshold[0]);
        break;
      }
      case(filter2d::proportion):{
        stat.eraseNoData(windowBuffer);
	T2 sum=stat.sum(windowBuffer);
	if(sum)
	  outBuffer[x/down]=windowBuffer[centre]/sum;
	else
	  outBuffer[x/down]=noDataValue;
        break;
      }
      case(filter2d::homog):{
        T1 centreValue=inBuffer[(dimY-1)/2][x];
        bool isHomog=true;
        stat.eraseNoData(windowBuffer);
        typename std::vector<T1>::const_iterator wit;
        for(wit=windowBuffer.begin();wit!=windowBuffer.end();++wit){
          if(*wit==centreValue)
            continue;
          else{
            isHomog=false;
            break;
          }
        }
        if(isHomog)
          outBuffer[x/down]=1;
        else
          outBuffer[x/down]=noDataValue;
        break;
      }
      case(filter2d::heterog):{
        T1 centreValue=inBuffer[(dimY-1)/2][x];
        bool isHeterog=true;
        stat.eraseNoData(windowBuffer);
        typename std::vector<T1>::const_iterator wit;
        for(wit=windowBuffer.begin();wit!=windowBuffer.end();++wit){
          if(*wit!=centreValue)
            continue;
          else{
            isHeterog=false;
            break;
          }
        }
        if(isHeterog)
          outBuffer[x/down]=1;
        else
          outBuffer[x/down]=noDataValue;
        break;
      }
      case(filter2d::sauvola):{
        try{
          double theMean=0;
          double theStdev=0;
          bool invalid=false;
          T1 centreValue=inBuffer[(dimY-1)/2][x];
          if(windowBuffer.empty()||stat.isNoData(centreValue)){
            invalid=true;
            throw(invalid);
          }
          stat.meanVar(windowBuffer,theMean,theStdev);
          theStdev=sqrt(theStdev);
          double kValue=0.5;
          double rValue=128;
          if(m_threshold.size()==2){
            kValue=m_threshold[0];
            rValue=m_threshold[1];
          }
          //from http://fiji.sc/Auto_Local_Threshold
          //pixel = ( pixel > mean * ( 1 + k * ( standard_deviation / r - 1 ) ) ) ? object : background
          double theThreshold=theMean*(1+kValue*(theStdev/rValue - 1));
          //isdata value hardcoded as 1 for now
          outBuffer[x/down]=(centreValue>theThreshold) ? 1 : noDataValue;
        }
        catch(bool invalid){
          outBuffer[x/down]=noDataValue;
        }
        break;
      }
      case(filter2d::density):{
        int nvalid=stat.nvalid(windowBuffer);
        if(nvalid){
          std::vector<short>::const_iterator vit=m_class.begin();
          while(vit!=m_class.end())
            outBuffer[x/down]+=100.0*occurrence[*(vit++)]/nvalid;
        }
        else
          outBuffer[x/down]=noDataValue;
        break;
      }
      case(filter2d::countid):{
        if(occurrence.size())
	  outBuffer[x/down]=occurrence.size();
	else
	  outBuffer[x/down]=noDataValue;
	break;
      }
      case(filter2d::mode):{
        if(occurrence.size()){
          std::map<int,int>::const_iterator maxit=occurrence.begin();
          for(std::map<int,int>::const_iterator mit=occurrence.begin();mit!=occurrence.end();++mit){
            if(mit->second>maxit->second)
              maxit=mit;
          }
          if(occurrence[inBuffer[(dimY-1)/2][x]]<maxit->second)//
            outBuffer[x/down]=maxit->first;
          else//favorize original value in case of ties
            outBuffer[x/down]=inBuffer[(dimY-1)/2][x];
        }
        else
          outBuffer[x/down]=noDataValue;
        break;
      }
      case(filter2d::threshold):{
        assert(m_class.size()==m_threshold.size());
        int nvalid=stat.nvalid(windowBuffer);
        if(nvalid>0){
          outBuffer[x/down]=inBuffer[(dimY-1)/2][x];//initialize with original value (in case thresholds not met)
          for(int iclass=0;iclass<m_class.size();++iclass){
            if(100.0*(occurrence[m_class[iclass]])/nvalid>m_threshold[iclass])
              outBuffer[x/down]=m_class[iclass];
          }
        }
        else
          outBuffer[x/down]=noDataValue;
        break;
      }
      case(filter2d::scramble):{//could be done more efficiently window by window with random shuffling entire buffer and assigning entire buffer at once to output image...
	if(windowBuffer.size()){
	  int randomIndex=std::rand()%windowBuffer.size();
	  if(randomIndex>=windowBuffer.size())
	    outBuffer[x/down]=windowBuffer.back();
	  else if(randomIndex<0)
	    outBuffer[x/down]=windowBuffer[0];
	  else
	    outBuffer[x/down]=windowBuffer[randomIndex];
	}
	else
	  outBuffer[x/down]=noDataValue;
	break;
      }
      case(filter2d::mixed):{
        enum MixType { BF=11, CF=12, MF=13, NF=20, W=30 };
        double nBF=occurrence[BF];
        double nCF=occurrence[CF];
        double nMF=occurrence[MF];
        double nNF=occurrence[NF];
        double nW=occurrence[W];
        if(windowBuffer.size()){
          if((nBF+nCF+nMF)&&(nBF+nCF+nMF>=nNF+nW)){//forest
            if(nBF/(nBF+nCF)>=0.75)
              outBuffer[x/down]=BF;
            else if(nCF/(nBF+nCF)>=0.75)
              outBuffer[x/down]=CF;
            else
              outBuffer[x/down]=MF;
          }
          else{//non-forest
            if(nW&&(nW>=nNF))
              outBuffer[x/down]=W;
            else
              outBuffer[x/down]=NF;
          }
        }
        else
          outBuffer[x/down]=inBuffer[indexJ][indexI];
        break;
      }
      default:
        break;
      }
    }
    progress=(1.0+y/down);
    progress/=outputVector.size();
    pfnProgress(progress,pszMessage,pProgressArg);
    //copy outBuffer to outputVector
    outputVector[y/down]=outBuffer;
  }
}

// class Compare_mapValue{
// public:
//   int operator() (const map<int,int>::value_type& v1, const map<int, int>::value_type& v2) const{
//     return (v1.second)>(v2.second);
//   }
// };

template<class T> void Filter2d::shift(const Vector2d<T>& input, Vector2d<T>& output, double offsetX, double offsetY, double randomSigma, RESAMPLE resample, bool verbose){
  output.resize(input.nRows(),input.nCols());
  const gsl_rng_type *rangenType;
  gsl_rng *rangen;
  gsl_rng_env_setup();
  rangenType=gsl_rng_default;
  rangen=gsl_rng_alloc(rangenType);
  long seed=time(NULL)*getpid();
  gsl_rng_set(rangen,seed);
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int j=0;j<input.nRows();++j){
    for(int i=0;i<input.nCols();++i){
      T theValue=0;
      double randomX=0;
      double randomY=0;
      if(randomSigma>0){
        randomX=gsl_ran_gaussian(rangen,randomSigma);
        randomY=gsl_ran_gaussian(rangen,randomSigma);
      }
      double readCol=i+offsetX+randomX;
      double readRow=j+offsetY+randomY;
      if(readRow<0)
        readRow=0;
      if(readRow>input.nRows()-1)
        readRow=input.nRows()-1;
      if(readCol<0)
        readCol=0;
      if(readCol>input.nCols()-1)
        readCol=input.nCols()-1;
      switch(resample){
      case(BILINEAR):{
        double lowerRow=readRow-0.5;
        double upperRow=readRow+0.5;
        lowerRow=static_cast<int>(lowerRow);
        upperRow=static_cast<int>(upperRow);
        double lowerCol=readCol-0.5;
        double upperCol=readCol+0.5;
        lowerCol=static_cast<int>(lowerCol);
        upperCol=static_cast<int>(upperCol);
        assert(lowerRow>=0);
        assert(lowerRow<input.nRows());
        assert(lowerCol>=0);
        assert(lowerCol<input.nCols());
        assert(upperRow>=0);
        assert(upperRow<input.nRows());
        assert(upperCol>=0);
        if(upperCol>=input.nCols()){
          std::cout << "upperCol: " << upperCol << std::endl;
          std::cout << "readCol: " << readCol << std::endl;
          std::cout << "readCol+0.5: " << readCol+0.5 << std::endl;
          std::cout << "static_cast<int>(readCol+0.5): " << static_cast<int>(readCol+0.5) << std::endl;
        }
        assert(upperCol<input.nCols());
        double c00=input[lowerRow][lowerCol];
        double c11=input[upperRow][upperCol];
        double c01=input[lowerRow][upperCol];
        double c10=input[upperRow][lowerCol];
        double a=(upperCol-readCol)*c00+(readCol-lowerCol)*c01;
        double b=(upperCol-readCol)*c10+(readCol-lowerCol)*c11;
        theValue=(upperRow-readRow)*a+(readRow-lowerRow)*b;
        break;
      }
      default:
        theValue=input[static_cast<int>(readRow)][static_cast<int>(readCol)];
        break;
      }
      assert(j>=0);
      assert(j<output.nRows());
      assert(i>=0);
      assert(i<output.nCols());
      output[j][i]=theValue;
    }
    progress=(1.0+j);
    progress/=output.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  gsl_rng_free(rangen);
}

template<class T> unsigned long int Filter2d::morphology(const Vector2d<T>& input, Vector2d<T>& output, const std::string& method, int dimX, int dimY, bool disc, double hThreshold)
{
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  double noDataValue=0;
  if(m_noDataValues.size())
    noDataValue=m_noDataValues[0];

  unsigned long int nchange=0;
  assert(dimX);
  assert(dimY);
  statfactory::StatFactory stat;
  Vector2d<T> inBuffer(dimY,input.nCols());
  output.clear();
  output.resize(input.nRows(),input.nCols());
  int indexI=0;
  int indexJ=0;
  //initialize last half of inBuffer
  for(int j=-(dimY-1)/2;j<=dimY/2;++j){
    for(int i=0;i<input.nCols();++i)
      inBuffer[indexJ][i]=input[abs(j)][i];
    ++indexJ;
  }
  for(int y=0;y<input.nRows();++y){
    if(y){//inBuffer already initialized for y=0
      //erase first line from inBuffer
      inBuffer.erase(inBuffer.begin());
      //read extra line and push back to inBuffer if not out of bounds
      if(y+dimY/2<input.nRows()){
        //allocate buffer
        inBuffer.push_back(inBuffer.back());
        for(int i=0;i<input.nCols();++i)
          inBuffer[inBuffer.size()-1][i]=input[y+dimY/2][i];
      }
      else{
        int over=y+dimY/2-input.nRows();
        int index=(inBuffer.size()-1)-over;
        assert(index>=0);
        assert(index<inBuffer.size());
        inBuffer.push_back(inBuffer[index]);
      }
    }
    for(int x=0;x<input.nCols();++x){
      output[y][x]=0;
      double currentValue=inBuffer[(dimY-1)/2][x];
      std::vector<double> statBuffer;
      bool currentMasked=false;
      for(int imask=0;imask<m_noDataValues.size();++imask){
        if(currentValue==m_noDataValues[imask]){
          currentMasked=true;
          break;
        }
      }
      output[y][x]=currentValue;//introduced due to hThreshold
      if(currentMasked){
        output[y][x]=currentValue;
      }
      else{
        for(int j=-(dimY-1)/2;j<=dimY/2;++j){
          for(int i=-(dimX-1)/2;i<=dimX/2;++i){
            double d2=i*i+j*j;//square distance
            if(disc&&(d2>(dimX/2)*(dimY/2)))
              continue;
            indexI=x+i;
            //check if out of bounds
            if(indexI<0)
              indexI=-indexI;
            else if(indexI>=input.nCols())
              indexI=input.nCols()-i;
            if(y+j<0)
              indexJ=-j;
            else if(y+j>=input.nRows())
                indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
              else
                indexJ=(dimY-1)/2+j;
            if(inBuffer[indexJ][indexI]==noDataValue)
              continue;
            bool masked=false;
            for(int imask=0;imask<m_noDataValues.size();++imask){
              if(inBuffer[indexJ][indexI]==m_noDataValues[imask]){
                masked=true;
                break;
              }
            }
            if(!masked){
              short binValue=0;
              for(int iclass=0;iclass<m_class.size();++iclass){
                if(inBuffer[indexJ][indexI]==m_class[iclass]){
                  binValue=1;
                  break;
                }
              }
              if(m_class.size())
                statBuffer.push_back(binValue);
              else
                statBuffer.push_back(inBuffer[indexJ][indexI]);
            }
          }
        }
        if(statBuffer.size()){
          switch(getFilterType(method)){
          case(filter2d::dilate):
            if(output[y][x]<stat.mymax(statBuffer)-hThreshold){
              output[y][x]=stat.mymax(statBuffer);
              ++nchange;
            }
            break;
          case(filter2d::erode):
            if(output[y][x]>stat.mymin(statBuffer)+hThreshold){
              output[y][x]=stat.mymin(statBuffer);
              ++nchange;
            }
            break;
          default:
            std::ostringstream ess;
            ess << "Error:  morphology method " << method << " not supported, choose " << filter2d::dilate << " (dilate) or " << filter2d::erode << " (erode)" << std::endl;
            throw(ess.str());
            break;
          }
          if(output[y][x]&&m_class.size())
            output[y][x]=m_class[0];
          // else{
          //   assert(m_noDataValues.size());
          //   output[x]=m_noDataValues[0];
          // }
        }
        else
          output[y][x]=noDataValue;
      }
    }
    progress=(1.0+y);
    progress/=output.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  return nchange;
}

 template<class T> unsigned long int Filter2d::dsm2dtm_nwse(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim)
{
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  Vector2d<T> tmpDSM(inputDSM);
  double noDataValue=0;
  if(m_noDataValues.size())
    noDataValue=m_noDataValues[0];

  unsigned long int nchange=0;
  int dimX=dim;
  int dimY=dim;
  assert(dimX);
  assert(dimY);
  statfactory::StatFactory stat;
  Vector2d<T> inBuffer(dimY,inputDSM.nCols());
  if(outputMask.size()!=inputDSM.nRows())
    outputMask.resize(inputDSM.nRows());
  int indexI=0;
  int indexJ=0;
  //initialize last half of inBuffer
  for(int j=-(dimY-1)/2;j<=dimY/2;++j){
    for(int i=0;i<inputDSM.nCols();++i)
      inBuffer[indexJ][i]=tmpDSM[abs(j)][i];
    ++indexJ;
  }
  for(int y=0;y<tmpDSM.nRows();++y){
    if(y){//inBuffer already initialized for y=0
      //erase first line from inBuffer
      inBuffer.erase(inBuffer.begin());
      //read extra line and push back to inBuffer if not out of bounds
      if(y+dimY/2<tmpDSM.nRows()){
        //allocate buffer
        inBuffer.push_back(inBuffer.back());
        for(int i=0;i<tmpDSM.nCols();++i)
          inBuffer[inBuffer.size()-1][i]=tmpDSM[y+dimY/2][i];
      }
      else{
        int over=y+dimY/2-tmpDSM.nRows();
        int index=(inBuffer.size()-1)-over;
        assert(index>=0);
        assert(index<inBuffer.size());
        inBuffer.push_back(inBuffer[index]);
      }
    }
    for(int x=0;x<tmpDSM.nCols();++x){
      double centerValue=inBuffer[(dimY-1)/2][x];
      short nmasked=0;
      std::vector<T> neighbors;
      for(int j=-(dimY-1)/2;j<=dimY/2;++j){
	for(int i=-(dimX-1)/2;i<=dimX/2;++i){
	  indexI=x+i;
	  //check if out of bounds
	  if(indexI<0)
	    indexI=-indexI;
	  else if(indexI>=tmpDSM.nCols())
	    indexI=tmpDSM.nCols()-i;
	  if(y+j<0)
	    indexJ=-j;
	  else if(y+j>=tmpDSM.nRows())
	    indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
	  else
	    indexJ=(dimY-1)/2+j;
	  double difference=(centerValue-inBuffer[indexJ][indexI]);
	  if(i||j)//skip centerValue
	    neighbors.push_back(inBuffer[indexJ][indexI]);
	  if(difference>hThreshold)
	    ++nmasked;
	}
      }
      if(nmasked<=nlimit){
	++nchange;
	//reset pixel in outputMask
	outputMask[y][x]=0;
      }
      else{
	//reset pixel height in tmpDSM
	sort(neighbors.begin(),neighbors.end());
	assert(neighbors.size()>1);
	inBuffer[(dimY-1)/2][x]=neighbors[1];
	/* inBuffer[(dimY-1)/2][x]=stat.mymin(neighbors); */
      }
    }
    progress=(1.0+y);
    progress/=outputMask.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  return nchange;
}

 template<class T> unsigned long int Filter2d::dsm2dtm_nesw(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim)
{
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  Vector2d<T> tmpDSM(inputDSM);
  double noDataValue=0;
  if(m_noDataValues.size())
    noDataValue=m_noDataValues[0];

  unsigned long int nchange=0;
  int dimX=dim;
  int dimY=dim;
  assert(dimX);
  assert(dimY);
  statfactory::StatFactory stat;
  Vector2d<T> inBuffer(dimY,inputDSM.nCols());
  if(outputMask.size()!=inputDSM.nRows())
    outputMask.resize(inputDSM.nRows());
  int indexI=0;
  int indexJ=0;
  //initialize last half of inBuffer
  for(int j=-(dimY-1)/2;j<=dimY/2;++j){
    for(int i=0;i<inputDSM.nCols();++i)
      inBuffer[indexJ][i]=tmpDSM[abs(j)][i];
    ++indexJ;
  }
  for(int y=0;y<tmpDSM.nRows();++y){
    if(y){//inBuffer already initialized for y=0
      //erase first line from inBuffer
      inBuffer.erase(inBuffer.begin());
      //read extra line and push back to inBuffer if not out of bounds
      if(y+dimY/2<tmpDSM.nRows()){
        //allocate buffer
        inBuffer.push_back(inBuffer.back());
        for(int i=0;i<tmpDSM.nCols();++i)
          inBuffer[inBuffer.size()-1][i]=tmpDSM[y+dimY/2][i];
      }
      else{
        int over=y+dimY/2-tmpDSM.nRows();
        int index=(inBuffer.size()-1)-over;
        assert(index>=0);
        assert(index<inBuffer.size());
        inBuffer.push_back(inBuffer[index]);
      }
    }
    for(int x=tmpDSM.nCols()-1;x>=0;--x){
      double centerValue=inBuffer[(dimY-1)/2][x];
      short nmasked=0;
      std::vector<T> neighbors;
      for(int j=-(dimY-1)/2;j<=dimY/2;++j){
	for(int i=-(dimX-1)/2;i<=dimX/2;++i){
	  indexI=x+i;
	  //check if out of bounds
	  if(indexI<0)
	    indexI=-indexI;
	  else if(indexI>=tmpDSM.nCols())
	    indexI=tmpDSM.nCols()-i;
	  if(y+j<0)
	    indexJ=-j;
	  else if(y+j>=tmpDSM.nRows())
	    indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
	  else
	    indexJ=(dimY-1)/2+j;
	  double difference=(centerValue-inBuffer[indexJ][indexI]);
	  if(i||j)//skip centerValue
	    neighbors.push_back(inBuffer[indexJ][indexI]);
	  if(difference>hThreshold)
	    ++nmasked;
	}
      }
      if(nmasked<=nlimit){
	++nchange;
	//reset pixel in outputMask
	outputMask[y][x]=0;
      }
      else{
	//reset pixel height in tmpDSM
	sort(neighbors.begin(),neighbors.end());
	assert(neighbors.size()>1);
	inBuffer[(dimY-1)/2][x]=neighbors[1];
	/* inBuffer[(dimY-1)/2][x]=stat.mymin(neighbors); */
      }
    }
    progress=(1.0+y);
    progress/=outputMask.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  return nchange;
}

 template<class T> unsigned long int Filter2d::dsm2dtm_senw(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim)
{
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  Vector2d<T> tmpDSM(inputDSM);
  double noDataValue=0;
  if(m_noDataValues.size())
    noDataValue=m_noDataValues[0];

  unsigned long int nchange=0;
  int dimX=dim;
  int dimY=dim;
  assert(dimX);
  assert(dimY);
  statfactory::StatFactory stat;
  Vector2d<T> inBuffer(dimY,inputDSM.nCols());
  if(outputMask.size()!=inputDSM.nRows())
    outputMask.resize(inputDSM.nRows());
  int indexI=0;
  int indexJ=inputDSM.nRows()-1;
  //initialize first half of inBuffer
  for(int j=inputDSM.nRows()-dimY/2;j<inputDSM.nRows();--j){
    for(int i=0;i<inputDSM.nCols();++i)
      inBuffer[indexJ][i]=tmpDSM[abs(j)][i];
    ++indexJ;
  }
  for(int y=tmpDSM.nRows()-1;y>=0;--y){
    if(y<tmpDSM.nRows()-1){//inBuffer already initialized for y=tmpDSM.nRows()-1
      //erase last line from inBuffer
      inBuffer.erase(inBuffer.end()-1);
      //read extra line and insert to inBuffer if not out of bounds
      if(y-dimY/2>0){
        //allocate buffer
        inBuffer.insert(inBuffer.begin(),inBuffer.back());
        for(int i=0;i<tmpDSM.nCols();++i)
          inBuffer[0][i]=tmpDSM[y-dimY/2][i];
      }
      else{
        inBuffer.insert(inBuffer.begin(),inBuffer[abs(y-dimY/2)]);
      }
    }
    for(int x=tmpDSM.nCols()-1;x>=0;--x){
      double centerValue=inBuffer[(dimY-1)/2][x];
      short nmasked=0;
      std::vector<T> neighbors;
      for(int j=-(dimY-1)/2;j<=dimY/2;++j){
	for(int i=-(dimX-1)/2;i<=dimX/2;++i){
	  indexI=x+i;
	  //check if out of bounds
	  if(indexI<0)
	    indexI=-indexI;
	  else if(indexI>=tmpDSM.nCols())
	    indexI=tmpDSM.nCols()-i;
	  if(y+j<0)
	    indexJ=-j;
	  else if(y+j>=tmpDSM.nRows())
	    indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
	  else
	    indexJ=(dimY-1)/2+j;
	  double difference=(centerValue-inBuffer[indexJ][indexI]);
	  if(i||j)//skip centerValue
	    neighbors.push_back(inBuffer[indexJ][indexI]);
	  if(difference>hThreshold)
	    ++nmasked;
	}
      }
      if(nmasked<=nlimit){
	++nchange;
	//reset pixel in outputMask
	outputMask[y][x]=0;
      }
      else{
	//reset pixel height in tmpDSM
	sort(neighbors.begin(),neighbors.end());
	assert(neighbors.size()>1);
	inBuffer[(dimY-1)/2][x]=neighbors[1];
	/* inBuffer[(dimY-1)/2][x]=stat.mymin(neighbors); */
      }
    }
    progress=(1.0+y);
    progress/=outputMask.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  return nchange;
}

 template<class T> unsigned long int Filter2d::dsm2dtm_swne(const Vector2d<T>& inputDSM, Vector2d<T>& outputMask, double hThreshold, int nlimit, int dim)
{
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  Vector2d<T> tmpDSM(inputDSM);
  double noDataValue=0;
  if(m_noDataValues.size())
    noDataValue=m_noDataValues[0];

  unsigned long int nchange=0;
  int dimX=dim;
  int dimY=dim;
  assert(dimX);
  assert(dimY);
  statfactory::StatFactory stat;
  Vector2d<T> inBuffer(dimY,inputDSM.nCols());
  if(outputMask.size()!=inputDSM.nRows())
    outputMask.resize(inputDSM.nRows());
  int indexI=0;
  int indexJ=0;
  //initialize first half of inBuffer
  for(int j=inputDSM.nRows()-dimY/2;j<inputDSM.nRows();--j){
    for(int i=0;i<inputDSM.nCols();++i)
      inBuffer[indexJ][i]=tmpDSM[abs(j)][i];
    ++indexJ;
  }
  for(int y=tmpDSM.nRows()-1;y>=0;--y){
    if(y<tmpDSM.nRows()-1){//inBuffer already initialized for y=0
      //erase last line from inBuffer
      inBuffer.erase(inBuffer.end()-1);
      //read extra line and insert to inBuffer if not out of bounds
      if(y-dimY/2>0){
        //allocate buffer
        inBuffer.insert(inBuffer.begin(),inBuffer.back());
        for(int i=0;i<tmpDSM.nCols();++i)
          inBuffer[0][i]=tmpDSM[y-dimY/2][i];
      }
      else{
        inBuffer.insert(inBuffer.begin(),inBuffer[abs(y-dimY/2)]);
      }
    }
    for(int x=0;x<tmpDSM.nCols();++x){
      double centerValue=inBuffer[(dimY-1)/2][x];
      short nmasked=0;
      std::vector<T> neighbors;
      for(int j=-(dimY-1)/2;j<=dimY/2;++j){
	for(int i=-(dimX-1)/2;i<=dimX/2;++i){
	  indexI=x+i;
	  //check if out of bounds
	  if(indexI<0)
	    indexI=-indexI;
	  else if(indexI>=tmpDSM.nCols())
	    indexI=tmpDSM.nCols()-i;
	  if(y+j<0)
	    indexJ=-j;
	  else if(y+j>=tmpDSM.nRows())
	    indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
	  else
	    indexJ=(dimY-1)/2+j;
	  double difference=(centerValue-inBuffer[indexJ][indexI]);
	  if(i||j)//skip centerValue
	    neighbors.push_back(inBuffer[indexJ][indexI]);
	  if(difference>hThreshold)
	    ++nmasked;
	}
      }
      if(nmasked<=nlimit){
	++nchange;
	//reset pixel in outputMask
	outputMask[y][x]=0;
      }
      else{
	//reset pixel height in tmpDSM
	sort(neighbors.begin(),neighbors.end());
	assert(neighbors.size()>1);
	inBuffer[(dimY-1)/2][x]=neighbors[1];
	/* inBuffer[(dimY-1)/2][x]=stat.mymin(neighbors); */
      }
    }
    progress=(1.0+y);
    progress/=outputMask.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  return nchange;
}

  template<class T> void Filter2d::shadowDsm(const Vector2d<T>& input, Vector2d<T>& output, double sza, double saa, double pixelSize, short shadowFlag)
{
  unsigned int ncols=input.nCols();
  output.clear();
  output.resize(input.nRows(),ncols);
  //do we need to initialize output?
  // for(int y=0;y<output.nRows();++y)
  //   for(int x=0;x<output.nCols();++x)
  //     output[y][x]=0;
  int indexI=0;
  int indexJ=0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int y=0;y<input.nRows();++y){
    for(int x=0;x<input.nCols();++x){
      double currentValue=input[y][x];
      int theDist=static_cast<int>(sqrt((currentValue*tan(DEG2RAD(sza))/pixelSize)*(currentValue*tan(DEG2RAD(sza))/pixelSize)));//in pixels
      double theDir=DEG2RAD(saa)+PI/2.0;
      if(theDir<0)
        theDir+=2*PI;
      for(int d=0;d<theDist;++d){//d in pixels
        indexI=x+d*cos(theDir);//in pixels
        indexJ=y+d*sin(theDir);//in pixels
        if(indexJ<0||indexJ>=input.size())
          continue;
        if(indexI<0||indexI>=input[indexJ].size())
          continue;
        if(input[indexJ][indexI]<currentValue-d*pixelSize/tan(DEG2RAD(sza))){//in m
          output[indexJ][indexI]=shadowFlag;
        }
      }
    }
    progress=(1.0+y);
    progress/=output.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

template<class T> void Filter2d::dwtForward(Vector2d<T>& theBuffer, const std::string& wavelet_type, int family){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  int nRow=theBuffer.size();
  assert(nRow);
  int nCol=theBuffer[0].size();
  assert(nCol);
  //make sure data size if power of 2
  while(theBuffer.size()&(theBuffer.size()-1))
    theBuffer.push_back(theBuffer.back());
  for(int irow=0;irow<theBuffer.size();++irow)
    while(theBuffer[irow].size()&(theBuffer[irow].size()-1))
      theBuffer[irow].push_back(theBuffer[irow].back());
  std::vector<double> vdata(theBuffer.size()*theBuffer[0].size());
  double* data=&(vdata[0]);
  for(int irow=0;irow<theBuffer.size();++irow){
    for(int icol=0;icol<theBuffer[0].size();++icol){
      int index=irow*theBuffer[0].size()+icol;
      data[index]=theBuffer[irow][icol];
    }
  }
  int nsize=theBuffer.size()*theBuffer[0].size();
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  assert(nsize);
  w=gsl_wavelet_alloc(filter::Filter::getWaveletType(wavelet_type),family);
  work=gsl_wavelet_workspace_alloc(nsize);
  gsl_wavelet2d_nstransform_forward (w, data, theBuffer.size(), theBuffer.size(),theBuffer[0].size(), work);
  theBuffer.erase(theBuffer.begin()+nRow,theBuffer.end());
  for(int irow=0;irow<theBuffer.size();++irow){
    theBuffer[irow].erase(theBuffer[irow].begin()+nCol,theBuffer[irow].end());
    for(int icol=0;icol<theBuffer[irow].size();++icol){
      int index=irow*theBuffer[irow].size()+icol;
      theBuffer[irow][icol]=data[index];
    }
    progress=(1.0+irow);
    progress/=theBuffer.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
}

template<class T> void Filter2d::dwtInverse(Vector2d<T>& theBuffer, const std::string& wavelet_type, int family){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  int nRow=theBuffer.size();
  assert(nRow);
  int nCol=theBuffer[0].size();
  assert(nCol);
  //make sure data size if power of 2
  while(theBuffer.size()&(theBuffer.size()-1))
    theBuffer.push_back(theBuffer.back());
  for(int irow=0;irow<theBuffer.size();++irow)
    while(theBuffer[irow].size()&(theBuffer[irow].size()-1))
      theBuffer[irow].push_back(theBuffer[irow].back());
  std::vector<double> vdata(theBuffer.size()*theBuffer[0].size());
  double* data=&(vdata[0]);
  //double data[theBuffer.size()*theBuffer[0].size()];
  for(int irow=0;irow<theBuffer.size();++irow){
    for(int icol=0;icol<theBuffer[0].size();++icol){
      int index=irow*theBuffer[0].size()+icol;
      data[index]=theBuffer[irow][icol];
    }
  }
  int nsize=theBuffer.size()*theBuffer[0].size();
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  assert(nsize);
  w=gsl_wavelet_alloc(filter::Filter::getWaveletType(wavelet_type),family);
  work=gsl_wavelet_workspace_alloc(nsize);
  gsl_wavelet2d_nstransform_inverse (w, data, theBuffer.size(), theBuffer.size(),theBuffer[0].size(), work);
  theBuffer.erase(theBuffer.begin()+nRow,theBuffer.end());
  for(int irow=0;irow<theBuffer.size();++irow){
    theBuffer[irow].erase(theBuffer[irow].begin()+nCol,theBuffer[irow].end());
    for(int icol=0;icol<theBuffer[irow].size();++icol){
      int index=irow*theBuffer[irow].size()+icol;
      theBuffer[irow][icol]=data[index];
    }
    progress=(1.0+irow);
    progress/=theBuffer.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);
}

template<class T> void Filter2d::dwtCut(Vector2d<T>& theBuffer, const std::string& wavelet_type, int family, double cut){
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);

  int nRow=theBuffer.size();
  assert(nRow);
  int nCol=theBuffer[0].size();
  assert(nCol);
  //make sure data size if power of 2
  while(theBuffer.size()&(theBuffer.size()-1))
    theBuffer.push_back(theBuffer.back());
  for(int irow=0;irow<theBuffer.size();++irow)
    while(theBuffer[irow].size()&(theBuffer[irow].size()-1))
      theBuffer[irow].push_back(theBuffer[irow].back());
  double* data=new double[theBuffer.size()*theBuffer[0].size()];
  double* abscoeff=new double[theBuffer.size()*theBuffer[0].size()];
  size_t* p=new size_t[theBuffer.size()*theBuffer[0].size()];
  for(int irow=0;irow<theBuffer.size();++irow){
    for(int icol=0;icol<theBuffer[0].size();++icol){
      int index=irow*theBuffer[0].size()+icol;
      assert(index<theBuffer.size()*theBuffer[0].size());
      data[index]=theBuffer[irow][icol];
    }
  }
  int nsize=theBuffer.size()*theBuffer[0].size();
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;
  assert(nsize);
  w=gsl_wavelet_alloc(filter::Filter::getWaveletType(wavelet_type),family);
  work=gsl_wavelet_workspace_alloc(nsize);
  gsl_wavelet2d_nstransform_forward (w, data, theBuffer.size(), theBuffer[0].size(),theBuffer[0].size(), work);
  for(int irow=0;irow<theBuffer.size();++irow){
    for(int icol=0;icol<theBuffer[0].size();++icol){
      int index=irow*theBuffer[0].size()+icol;
      abscoeff[index]=fabs(data[index]);
    }
  }
  int nc=(100-cut)/100.0*nsize;
  gsl_sort_index(p,abscoeff,1,nsize);
  for(int i=0;(i+nc)<nsize;i++)
   data[p[i]]=0;
  gsl_wavelet2d_nstransform_inverse (w, data, theBuffer.size(), theBuffer[0].size(),theBuffer[0].size(), work);
  for(int irow=0;irow<theBuffer.size();++irow){
    for(int icol=0;icol<theBuffer[irow].size();++icol){
      int index=irow*theBuffer[irow].size()+icol;
      theBuffer[irow][icol]=data[index];
    }
    progress=(1.0+irow);
    progress/=theBuffer.nRows();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  theBuffer.erase(theBuffer.begin()+nRow,theBuffer.end());
  for(int irow=0;irow<theBuffer.size();++irow)
    theBuffer[irow].erase(theBuffer[irow].begin()+nCol,theBuffer[irow].end());
  delete[] data;
  delete[] abscoeff;
  delete[] p;
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

}

}

#endif /* _MYFILTER_H_ */

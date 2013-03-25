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

#include <assert.h>
#include <math.h>
#include <limits>
#include <vector>
#include <string>
#include <map>
#include "base/Vector2d.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "algorithms/StatFactory.h"

using namespace std;
// using namespace cimg_library;
namespace filter2d
{
  enum FILTER_TYPE { median=0, var=1 , min=2, max=3, sum=4, mean=5, minmax=6, dilate=7, erode=8, close=9, open=10, homog=11, sobelx=12, sobely=13, sobelxy=14, sobelyx=-14, smooth=15, density=16, majority=17, mixed=18, smoothnodata=19, threshold=20, ismin=21, ismax=22, heterog=23, order=24, stdev=25, mrf=26};
  
class Filter2d
{
public:
  Filter2d(void);
  Filter2d(const Vector2d<double> &taps);
  virtual ~Filter2d(){};
  static FILTER_TYPE getFilterType(const std::string filterType){
    std::map<std::string, FILTER_TYPE> m_filterMap;
    initMap(m_filterMap);
    return m_filterMap[filterType];
  };

  void setTaps(const Vector2d<double> &taps);
  void setNoValue(double noValue=0){m_noValue=noValue;};
  void pushClass(short theClass=1){m_class.push_back(theClass);};
  void pushMask(short theMask=0){m_mask.push_back(theMask);};
  void pushThreshold(double theThreshold){m_threshold.push_back(theThreshold);};
  void setThresholds(const vector<double>& theThresholds){m_threshold=theThresholds;};
  void setClasses(const vector<short>& theClasses){m_class=theClasses;};
  void filter(const ImgReaderGdal& input, ImgWriterGdal& output, bool absolute=false, bool normalize=true, bool noData=false);
  void smooth(const ImgReaderGdal& input, ImgWriterGdal& output,int dim);
  void smooth(const ImgReaderGdal& input, ImgWriterGdal& output,int dimX, int dimY);
  void smoothNoData(const ImgReaderGdal& input, ImgWriterGdal& output,int dim);
  void smoothNoData(const ImgReaderGdal& input, ImgWriterGdal& output,int dimX, int dimY);
  template<class T1, class T2> void filter(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector);
  template<class T1, class T2> void smooth(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector,int dim);
  template<class T1, class T2> void smooth(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector,int dimX, int dimY);
  void majorVoting(const string& inputFilename, const string& outputFilename,int dim=0,const vector<int> &prior=vector<int>());
  /* void homogeneousSpatial(const string& inputFilename, const string& outputFilename, int dim, bool disc=false, int noValue=0); */
  void doit(const ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dim, short down=2, bool disc=false);
  void doit(const ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dimX, int dimY, short down=1, bool disc=false);
  void mrf(const ImgReaderGdal& input, ImgWriterGdal& output, int dimX, int dimY, double beta, bool eightConnectivity=true, short down=1);
  template<class T1, class T2> void doit(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector, const std::string& method, int dimX, int dimY, short down=1, bool disc=false);
  void median(const string& inputFilename, const string& outputFilename, int dim, bool disc=false);
  void var(const string& inputFilename, const string& outputFilename, int dim, bool disc=false);
  void morphology(const ImgReaderGdal& input, ImgWriterGdal& output, const std::string& method, int dimX, int dimY, bool disc=false, double angle=-190);
  template<class T> unsigned long int morphology(const Vector2d<T>& input, Vector2d<T>& output, const std::string& method, int dimX, int dimY, bool disc=false, double hThreshold=0);
  template<class T> void shadowDsm(const Vector2d<T>& input, Vector2d<T>& output, double sza, double saa, double pixelSize, short shadowFlag=1);
  void shadowDsm(const ImgReaderGdal& input, ImgWriterGdal& output, double sza, double saa, double pixelSize, short shadowFlag=1);
  void dwt_texture(const string& inputFilename, const string& outputFilename,int dim, int scale, int down=1, int iband=0, bool verbose=false);
  
private:
  static void initMap(std::map<std::string, FILTER_TYPE>& m_filterMap){
    //initialize selMap
    m_filterMap["mrf"]=filter2d::mrf;
    m_filterMap["stdev"]=filter2d::stdev;
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
    m_filterMap["majority"]=filter2d::majority;
    m_filterMap["mixed"]=filter2d::mixed;
    m_filterMap["smoothnodata"]=filter2d::smoothnodata;
    m_filterMap["threshold"]=filter2d::threshold;
    m_filterMap["ismin"]=filter2d::ismin;
    m_filterMap["ismax"]=filter2d::ismax;
    m_filterMap["heterog"]=filter2d::heterog;
    m_filterMap["order"]=filter2d::order;
    m_filterMap["median"]=filter2d::median;
  }

  Vector2d<double> m_taps;
  double m_noValue;
  vector<short> m_class;
  vector<short> m_mask;
  vector<double> m_threshold;
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
    vector<T2> outBuffer(inputVector[0].size());
    //initialize last half of inBuffer
    int indexI=0;
    int indexJ=0;
    for(int y=0;y<dimY;++y){
      if(y<dimY/2)
        continue;//skip first half
      inBuffer[y]=inputVector[indexJ++];
    }
    for(int y=0;y<inputVector.size();++y){
      if(y){//inBuffer already initialized for y=0
        //erase first line from inBuffer
        inBuffer.erase(inBuffer.begin());
        //read extra line and push back to inBuffer if not out of bounds
        if(y+dimY/2<inputVector.size())
          inBuffer.push_back(inputVector[y+dimY/2]);
      }
      for(int x=0;x<inputVector[0].size();++x){
        outBuffer[x]=0;
        for(int j=-dimY/2;j<(dimY+1)/2;++j){
          for(int i=-dimX/2;i<(dimX+1)/2;++i){
            indexI=x+i;
            indexJ=dimY/2+j;
            //check if out of bounds
            if(x<dimX/2)
              indexI=x+abs(i);
            else if(x>=inputVector[0].size()-dimX/2)
              indexI=x-abs(i);
            if(y<dimY/2)
              indexJ=dimY/2+abs(j);
            else if(y>=inputVector.size()-dimY/2)
              indexJ=dimY/2-abs(j);
            outBuffer[x]+=(m_taps[dimY/2+j][dimX/2+i]*inBuffer[indexJ][indexI]);
          }
        }
      }
      //copy outBuffer to outputVector
      outputVector[y]=outBuffer;
    }
  }

template<class T1, class T2> void Filter2d::doit(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector, const std::string& method, int dimX, int dimY, short down, bool disc)
{
  statfactory::StatFactory stat;
  outputVector.resize((inputVector.size()+down-1)/down);
  Vector2d<T1> inBuffer(dimY);
  vector<T2> outBuffer((inputVector[0].size()+down-1)/down);
  //initialize last half of inBuffer
  int indexI=0;
  int indexJ=0;
  for(int y=0;y<dimY;++y){
    if(y<dimY/2)
      continue;//skip first half
    inBuffer[y]=inputVector[indexJ++];
  }
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
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
      vector<double> windowBuffer;
      map<int,int> occurrence;
      for(int j=-dimY/2;j<(dimY+1)/2;++j){
	for(int i=-dimX/2;i<(dimX+1)/2;++i){
	  indexI=x+i;
	  indexJ=dimY/2+j;
	  //check if out of bounds
	  if(x<dimX/2)
	    indexI=x+abs(i);
	  else if(x>=inputVector[0].size()-dimX/2)
	    indexI=x-abs(i);
	  if(y<dimY/2)
	    indexJ=dimY/2+abs(j);
	  else if(y>=inputVector.size()-dimY/2)
	    indexJ=dimY/2-abs(j);
          bool masked=false;
          for(int imask=0;imask<m_mask.size();++imask){
            if(inBuffer[indexJ][indexI]==m_mask[imask]){
              masked=true;
              break;
            }
          }
          if(!masked){
            vector<short>::const_iterator vit=m_class.begin();
            //todo: test if this works (only add occurrence if within defined classes)!
            if(!m_class.size())
              ++occurrence[inBuffer[indexJ][indexI]];
            else{
              while(vit!=m_class.end()){
                if(inBuffer[indexJ][indexI]==*(vit++))
                  ++occurrence[inBuffer[indexJ][indexI]];
              }
            }
            windowBuffer.push_back(inBuffer[indexJ][indexI]);
          }
        }
      }
      switch(getFilterType(method)){
      case(filter2d::median):
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=stat.median(windowBuffer);
        break;
      case(filter2d::var):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=stat.var(windowBuffer);
        break;
      }
      case(filter2d::stdev):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=sqrt(stat.var(windowBuffer));
        break;
      }
      case(filter2d::mean):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=stat.mean(windowBuffer);
        break;
      }
      case(filter2d::min):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=stat.min(windowBuffer);
        break;
      }
      case(filter2d::ismin):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=(stat.min(windowBuffer)==windowBuffer[dimX*dimY/2])? 1:0;
        break;
      }
      case(filter2d::minmax):{
        double min=0;
        double max=0;
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else{
          stat.minmax(windowBuffer,windowBuffer.begin(),windowBuffer.end(),min,max);
          if(min!=max)
            outBuffer[x/down]=0;
          else
            outBuffer[x/down]=windowBuffer[dimX*dimY/2];//centre pixels
        }
        break;
      }
      case(filter2d::max):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=stat.max(windowBuffer);
        break;
      }
      case(filter2d::ismax):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=(stat.max(windowBuffer)==windowBuffer[dimX*dimY/2])? 1:0;
        break;
      }
      case(filter2d::order):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else{
          double lbound=0;
          double ubound=dimX*dimY;
          double theMin=stat.min(windowBuffer);
          double theMax=stat.max(windowBuffer);
          double scale=(ubound-lbound)/(theMax-theMin);
          outBuffer[x/down]=static_cast<short>(scale*(windowBuffer[dimX*dimY/2]-theMin)+lbound);
        }
        break;
      }
      case(filter2d::sum):{
        outBuffer[x/down]=stat.sum(windowBuffer);
        break;
      }
      case(filter2d::homog):
        if(occurrence.size()==1)//all values in window must be the same
          outBuffer[x/down]=inBuffer[dimY/2][x];
        else//favorize original value in case of ties
          outBuffer[x/down]=m_noValue;
        break;
      case(filter2d::heterog):{
        for(vector<double>::const_iterator wit=windowBuffer.begin();wit!=windowBuffer.end();++wit){
          if(wit==windowBuffer.begin()+windowBuffer.size()/2)
            continue;
          else if(*wit!=inBuffer[dimY/2][x])
            outBuffer[x/down]=1;
          else if(*wit==inBuffer[dimY/2][x]){//todo:wit mag niet central pixel zijn
            outBuffer[x/down]=m_noValue;
            break;
          }
        }
        break;
      }
      case(filter2d::density):{
        if(windowBuffer.size()){
          vector<short>::const_iterator vit=m_class.begin();
          while(vit!=m_class.end())
            outBuffer[x/down]+=100.0*occurrence[*(vit++)]/windowBuffer.size();
        }
        else
          outBuffer[x/down]=m_noValue;
        break;
      }
      case(filter2d::majority):{
        if(occurrence.size()){
          map<int,int>::const_iterator maxit=occurrence.begin();
          for(map<int,int>::const_iterator mit=occurrence.begin();mit!=occurrence.end();++mit){
            if(mit->second>maxit->second)
              maxit=mit;
          }
          if(occurrence[inBuffer[dimY/2][x]]<maxit->second)//
            outBuffer[x/down]=maxit->first;
          else//favorize original value in case of ties
            outBuffer[x/down]=inBuffer[dimY/2][x];
        }
        else
          outBuffer[x/down]=m_noValue;
        break;
      }
      case(filter2d::threshold):{
        assert(m_class.size()==m_threshold.size());
        if(windowBuffer.size()){
          outBuffer[x/down]=inBuffer[dimY/2][x];//initialize with original value (in case thresholds not met)
          for(int iclass=0;iclass<m_class.size();++iclass){
            if(100.0*(occurrence[m_class[iclass]])/windowBuffer.size()>m_threshold[iclass])
              outBuffer[x/down]=m_class[iclass];
          }
        }
        else
          outBuffer[x/down]=m_noValue;
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
    progress+=(outputVector.size());
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

template<class T> unsigned long int Filter2d::morphology(const Vector2d<T>& input, Vector2d<T>& output, const std::string& method, int dimX, int dimY, bool disc, double hThreshold)
{
  unsigned long int nchange=0;
  assert(dimX);
  assert(dimY);
  statfactory::StatFactory stat;
  Vector2d<T> inBuffer(dimY,input.nCols());
  output.clear();
  output.resize(input.nRows(),input.nCols());
  //initialize last half of inBuffer
  int indexI=0;
  int indexJ=0;
  for(int j=-dimY/2;j<(dimY+1)/2;++j){
    for(int i=0;i<input.nCols();++i)
      inBuffer[indexJ][i]=input[abs(j)][i];
    ++indexJ;
  }
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
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
    }
    for(int x=0;x<input.nCols();++x){
      output[y][x]=0;
      double currentValue=inBuffer[dimY/2][x];
      vector<double> statBuffer;
      bool currentMasked=false;
      for(int imask=0;imask<m_mask.size();++imask){
        if(currentValue==m_mask[imask]){
          currentMasked=true;
          break;
        }
      }
      output[y][x]=currentValue;//introduced due to hThrehold
      if(currentMasked){
        output[y][x]=currentValue;
      }
      else{
        for(int j=-dimY/2;j<(dimY+1)/2;++j){
          for(int i=-dimX/2;i<(dimX+1)/2;++i){
            if(disc&&(i*i+j*j>(dimX/2)*(dimY/2)))
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
              indexJ=dimY/2-j;//indexJ=inBuffer.size()-1-j;
            else
              indexJ=dimY/2+j;
            if(inBuffer[indexJ][indexI]==m_noValue)
              continue;
            bool masked=false;
            for(int imask=0;imask<m_mask.size();++imask){
              if(inBuffer[indexJ][indexI]==m_mask[imask]){
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
            if(output[y][x]<stat.max(statBuffer)-hThreshold){
              output[y][x]=stat.max(statBuffer);
              ++nchange;
            }
            break;
          case(filter2d::erode):
            if(output[y][x]>stat.min(statBuffer)+hThreshold){
              output[y][x]=stat.min(statBuffer);
              ++nchange;
            }
            break;
          default:
            ostringstream ess;
            ess << "Error:  morphology method " << method << " not supported, choose " << filter2d::dilate << " (dilate) or " << filter2d::erode << " (erode)" << endl;
            throw(ess.str());
            break;
          }
          if(output[y][x]&&m_class.size())
            output[y][x]=m_class[0];
          // else{
          //   assert(m_mask.size());
          //   output[x]=m_mask[0];
          // }
        }
        else
          output[y][x]=m_noValue;
      }
    }
    progress=(1.0+y);
    progress/=output.nRows();
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

}

#endif /* _MYFILTER_H_ */

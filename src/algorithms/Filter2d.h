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

using namespace std;
// using namespace cimg_library;
namespace Filter2d
{
  enum Type { MEDIAN=0, VAR=1 , MIN=2, MAX=3, SUM=4, MEAN=5, MINMAX=6, DILATE=7, ERODE=8, CLOSE=9, OPEN=10, HOMOG=11, SOBELX=12, SOBELY=13, SOBELXY=14, SOBELYX=-14, SMOOTH=15, DENSITY=16, MAJORITY=17, MIXED=18, SMOOTHNODATA=19, THRESHOLD=20, ISMIN=21, ISMAX=22, HETEROG=23, ORDER=24, STDEV=25};
  
class Filter2d
{
public:
  Filter2d(void);
  Filter2d(const Vector2d<double> &taps);
  virtual ~Filter2d(){};
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
  void doit(const ImgReaderGdal& input, ImgWriterGdal& output, int method, int dim, short down=2, bool disc=false);
  void doit(const ImgReaderGdal& input, ImgWriterGdal& output, int method, int dimX, int dimY, short down=1, bool disc=false);
  template<class T1, class T2> void doit(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector, int method, int dimX, int dimY, short down=1, bool disc=false);
  void median(const string& inputFilename, const string& outputFilename, int dim, bool disc=false);
  void var(const string& inputFilename, const string& outputFilename, int dim, bool disc=false);
  void morphology(const ImgReaderGdal& input, ImgWriterGdal& output, int method, int dimX, int dimY, bool disc=false, double angle=-190);
  template<class T> unsigned long int morphology(const Vector2d<T>& input, Vector2d<T>& output, int method, int dimX, int dimY, bool disc=false, double hThreshold=0);
  template<class T> void shadowDsm(const Vector2d<T>& input, Vector2d<T>& output, double sza, double saa, double pixelSize, short shadowFlag=1);
  void shadowDsm(const ImgReaderGdal& input, ImgWriterGdal& output, double sza, double saa, double pixelSize, short shadowFlag=1);
  void dwt_texture(const string& inputFilename, const string& outputFilename,int dim, int scale, int down=1, int iband=0, bool verbose=false);
  
private:
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

template<class T1, class T2> void Filter2d::doit(const Vector2d<T1>& inputVector, Vector2d<T2>& outputVector, int method, int dimX, int dimY, short down, bool disc)
{
  Histogram hist;
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
      switch(method){
      case(MEDIAN):
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=hist.median(windowBuffer);
        break;
      case(VAR):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=hist.var(windowBuffer);
        break;
      }
      case(STDEV):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=sqrt(hist.var(windowBuffer));
        break;
      }
      case(MEAN):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=hist.mean(windowBuffer);
        break;
      }
      case(MIN):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=hist.min(windowBuffer);
        break;
      }
      case(ISMIN):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=(hist.min(windowBuffer)==windowBuffer[dimX*dimY/2])? 1:0;
        break;
      }
      case(MINMAX):{
        double min=0;
        double max=0;
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else{
          hist.minmax(windowBuffer,windowBuffer.begin(),windowBuffer.end(),min,max);
          if(min!=max)
            outBuffer[x/down]=0;
          else
            outBuffer[x/down]=windowBuffer[dimX*dimY/2];//centre pixels
        }
        break;
      }
      case(MAX):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=hist.max(windowBuffer);
        break;
      }
      case(ISMAX):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else
          outBuffer[x/down]=(hist.max(windowBuffer)==windowBuffer[dimX*dimY/2])? 1:0;
        break;
      }
      case(ORDER):{
        if(windowBuffer.empty())
          outBuffer[x/down]=m_noValue;
        else{
          double lbound=0;
          double ubound=dimX*dimY;
          double theMin=hist.min(windowBuffer);
          double theMax=hist.max(windowBuffer);
          double scale=(ubound-lbound)/(theMax-theMin);
          outBuffer[x/down]=static_cast<short>(scale*(windowBuffer[dimX*dimY/2]-theMin)+lbound);
        }
        break;
      }
      case(SUM):{
        outBuffer[x/down]=hist.sum(windowBuffer);
        break;
      }
      case(HOMOG):
        if(occurrence.size()==1)//all values in window must be the same
          outBuffer[x/down]=inBuffer[dimY/2][x];
        else//favorize original value in case of ties
          outBuffer[x/down]=m_noValue;
        break;
      case(HETEROG):{
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
      case(DENSITY):{
        if(windowBuffer.size()){
          vector<short>::const_iterator vit=m_class.begin();
          while(vit!=m_class.end())
            outBuffer[x/down]+=100.0*occurrence[*(vit++)]/windowBuffer.size();
        }
        else
          outBuffer[x/down]=m_noValue;
        break;
      }
      case(MAJORITY):{
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
      case(THRESHOLD):{
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
      case(MIXED):{
        enum Type { BF=11, CF=12, MF=13, NF=20, W=30 };
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

template<class T> unsigned long int Filter2d::morphology(const Vector2d<T>& input, Vector2d<T>& output, int method, int dimX, int dimY, bool disc, double hThreshold)
{
  unsigned long int nchange=0;
  assert(dimX);
  assert(dimY);
  Histogram hist;
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
      vector<double> histBuffer;
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
                histBuffer.push_back(binValue);
              else
                histBuffer.push_back(inBuffer[indexJ][indexI]);
            }
          }
        }
        if(histBuffer.size()){
          switch(method){
          case(DILATE):
            if(output[y][x]<hist.max(histBuffer)-hThreshold){
              output[y][x]=hist.max(histBuffer);
              ++nchange;
            }
            break;
          case(ERODE):
            if(output[y][x]>hist.min(histBuffer)+hThreshold){
              output[y][x]=hist.min(histBuffer);
              ++nchange;
            }
            break;
          default:
            ostringstream ess;
            ess << "Error:  morphology method " << method << " not supported, choose " << DILATE << " (dilate) or " << ERODE << " (erode)" << endl;
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

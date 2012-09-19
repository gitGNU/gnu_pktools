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

#include <assert.h>
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
  enum Type { MEDIAN=0, VAR=1 , MIN=2, MAX=3, SUM=4, MEAN=5, MINMAX=6, DILATE=7, ERODE=8, CLOSE=9, OPEN=10, HOMOG=11, SOBELX=12, SOBELY=13, SOBELXY=14, SOBELYX=-14, SMOOTH=15, DENSITY=16, MAJORITY=17, MIXED=18, SMOOTHNODATA=19, THRESHOLD=20, ISMIN=21, ISMAX=22, HETEROG=23, ORDER=24};
  
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
  void median(const string& inputFilename, const string& outputFilename, int dim, bool disc=false);
  void var(const string& inputFilename, const string& outputFilename, int dim, bool disc=false);
  void morphology(const ImgReaderGdal& input, ImgWriterGdal& output, int method, int dimX, int dimY, bool disc=false, double angle=-190);
  template<class T> unsigned long int morphology(const Vector2d<T>& input, Vector2d<T>& output, int method, int dimX, int dimY, bool disc=false, double hThreshold=0);
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

}

#endif /* _MYFILTER_H_ */

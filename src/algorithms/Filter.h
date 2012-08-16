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
#include "Histogram.h"
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

private:
  vector<double> m_taps;
  vector<short> m_class;
  vector<short> m_mask;
};

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
  Histogram hist;
  vector<T> histBuffer;
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
      histBuffer.push_back(binValue);
    else
      histBuffer.push_back(input[i]);
    for(int t=1;t<=dim/2;++t){
      binValue=0;
      for(int iclass=0;iclass<m_class.size();++iclass){
        if(input[i+t]==m_class[iclass]){
          binValue=m_class[0];
          break;
        }
      }
      if(m_class.size()){
        histBuffer.push_back(binValue);
        histBuffer.push_back(binValue);
      }
      else{
        histBuffer.push_back(input[i+t]);
        histBuffer.push_back(input[i+t]);
      }
    }
    assert(histBuffer.size()==dim);
    if((i-offset)%down){
      histBuffer.clear();
      continue;
    }
    switch(method){
    case(DILATE):
      output[(i-offset+down-1)/down]=hist.max(histBuffer);
      break;
    case(ERODE):
      output[(i-offset+down-1)/down]=hist.min(histBuffer);
      break;
    default:
      string errorString="method not supported";
      throw(errorString);
      break;
    }
    if(verbose){
      cout << "buffer: ";
      for(int ibuf=0;ibuf<histBuffer.size();++ibuf)
        cout << histBuffer[ibuf] << " ";
      cout << "->" << output[(i-offset+down-1)/down] << endl;
    }
  }
  //main
  histBuffer.clear();
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
        histBuffer.push_back(binValue);
      else
        histBuffer.push_back(input[i-dim/2+t]);
    }
    assert(histBuffer.size()==dim);
    if((i-offset)%down){
      histBuffer.clear();
      continue;
    }
    switch(method){
    case(DILATE):
      output[(i-offset+down-1)/down]=hist.max(histBuffer);
      break;
    case(ERODE):
      output[(i-offset+down-1)/down]=hist.min(histBuffer);
      break;
    default:
      string errorString="method not supported";
      throw(errorString);
      break;
    }
    if(verbose){
      cout << "buffer: ";
      for(int ibuf=0;ibuf<histBuffer.size();++ibuf)
        cout << histBuffer[ibuf] << " ";
      cout << "->" << output[(i-offset+down-1)/down] << endl;
    }
    histBuffer.clear();
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
        histBuffer.push_back(binValue);
      else
        histBuffer.push_back(input[i]);
      for(int t=1;t<=dim/2;++t){
        binValue=0;
        for(int iclass=0;iclass<m_class.size();++iclass){
          if(input[i-t]==m_class[iclass]){
            binValue=m_class[0];
            break;
          }
        }
        if(m_class.size()){
          histBuffer.push_back(binValue);
          histBuffer.push_back(binValue);
        }
        else{
          histBuffer.push_back(input[i-t]);
          histBuffer.push_back(input[i-t]);
        }
      }
    if((i-offset)%down){
      histBuffer.clear();
      continue;
    }
    switch(method){
    case(DILATE):
      output[(i-offset+down-1)/down]=hist.max(histBuffer);
      break;
    case(ERODE):
      output[(i-offset+down-1)/down]=hist.min(histBuffer);
      break;
    default:
      string errorString="method not supported";
      throw(errorString);
      break;
    }
    if(verbose){
      cout << "buffer: ";
      for(int ibuf=0;ibuf<histBuffer.size();++ibuf)
        cout << histBuffer[ibuf] << " ";
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
using namespace filter;

#endif /* _MYFILTER_H_ */

/**********************************************************************
ImgReaderGdal.h: class to read raster files using GDAL API library
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
#ifndef _IMGREADERGDAL_H_
#define _IMGREADERGDAL_H_

#include "ImgRasterGdal.h"
#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "base/Vector2d.h"

//--------------------------------------------------------------------------
class ImgReaderGdal : public virtual ImgRasterGdal
{
public:
  ImgReaderGdal(void);
  //memory (in MB) indicates maximum memory read in memory for this image
  ImgReaderGdal(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly, unsigned long int memory=0){open(filename, readMode, memory);};
  ImgReaderGdal(void* dataPointer, unsigned int ncol, unsigned int nrow, unsigned short nband, const GDALDataType& dataType){open(dataPointer,ncol,nrow,nband,dataType);};
  void open(void* dataPointer, unsigned int ncol, unsigned int nrow, unsigned short nband, const GDALDataType& dataType);
  void open(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly, unsigned long int memory=0);
  void setMemory(unsigned long int memory=0){initMem(memory);};
  ~ImgReaderGdal(void);
  //  ImgReaderGdal(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly){open(filename, readMode);};
  void setScale(double theScale, int band=0){
    if(m_scale.size()!=nrOfBand()){//initialize
      m_scale.resize(nrOfBand());
      for(int iband=0;iband<nrOfBand();++iband)
       m_scale[iband]=1.0;
    }
    m_scale[band]=theScale;
  }
  void setOffset(double theOffset, int band=0){
    if(m_offset.size()!=nrOfBand()){
      m_offset.resize(nrOfBand());
      for(int iband=0;iband<nrOfBand();++iband)
       m_offset[iband]=0.0;
    }
      m_offset[band]=theOffset;
  }

  void close(void);
  template<typename T> void readData(T& value,  int col, int row, int band=0);
  template<typename T> void readData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  template<typename T> void readData(std::vector<T>& buffer, int minCol, int maxCol, double row, int band=0, RESAMPLE resample=NEAR);
  template<typename T> void readDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  template<typename T> void readDataBlock(std::vector<T>& buffer , int minCol, int maxCol, int minRow, int maxRow, int band=0);
  template<typename T> void readData(std::vector<T>& buffer, int row, int band=0);
  template<typename T> void readData(std::vector<T>& buffer, double row, int band=0, RESAMPLE resample=NEAR);
  void getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue);
  void getMinMax(double& minValue, double& maxValue, int band=0);
  double getMin(int& col, int& row, int band=0);
  double getHistogram(std::vector<double>& histvector, double& min, double& max,unsigned int& nbin, int theBand=0, bool kde=false);
  double getMax(int& col, int& row, int band=0);
  void getRefPix(double& refX, double &refY, int band=0);
  void getRange(std::vector<short>& range, int Band=0);
  unsigned long int getNvalid(int band);

protected:
  void setCodec(const GDALAccess& readMode=GA_ReadOnly);
  std::vector<double> m_scale;
  std::vector<double> m_offset;
  unsigned int m_blockSize;
  std::vector<void *> m_data;
  std::vector<unsigned int> m_begin;//first line that has been read
  std::vector<unsigned int> m_end;//beyond last line read

private:
  void initMem(unsigned long int memory);
  bool m_deletePointer;
  bool readNewBlock(int row, int band);
};

//     adfGeoTransform[0] /* top left x */
//     adfGeoTransform[1] /* w-e pixel resolution */
//     adfGeoTransform[2] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[3] /* top left y */
//     adfGeoTransform[4] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[5] /* n-s pixel resolution */

//todo: make sure m_deletePointer set to false if memory = 0

template<typename T> void ImgReaderGdal::readData(T& value, int col, int row, int band)
{
  assert(band<nrOfBand()+1);
  assert(col<nrOfCol());
  assert(col>=0);
  assert(row<nrOfRow());
  assert(row>=0);
  if(m_data.size()){
    //only support random access reading if entire image is in memory for performance reasons
    if(m_blockSize!=nrOfRow()){
      std::ostringstream s;
      s << "Error: increase memory to support random access reading (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
      throw(s.str());
    }
    if(row<m_begin[band]||row>=m_end[band])
      readNewBlock(row,band);
    int index=(row-m_begin[band])*nrOfCol()+col;
    double theScale=1;
    double theOffset=0;
    if(m_scale.size()>band||m_offset.size()>band){
      if(m_scale.size()>band)
	theScale=m_scale[band];
      if(m_offset.size()>band)
	theOffset=m_offset[band];
    }
    double dvalue=0;
    switch(getDataType()){
    case(GDT_Byte):
      dvalue=theScale*(static_cast<unsigned char*>(m_data[band])[index])+theOffset;
      break;
    case(GDT_Int16):
      dvalue=theScale*(static_cast<short*>(m_data[band])[index])+theOffset;
      break;
    case(GDT_UInt16):
      dvalue=theScale*(static_cast<unsigned short*>(m_data[band])[index])+theOffset;
      break;
    case(GDT_Int32):
      dvalue=theScale*(static_cast<int*>(m_data[band])[index])+theOffset;
      break;
    case(GDT_UInt32):
      dvalue=theScale*(static_cast<unsigned int*>(m_data[band])[index])+theOffset;
      break;
    case(GDT_Float32):
      dvalue=theScale*(static_cast<float*>(m_data[band])[index])+theOffset;
      break;
    case(GDT_Float64):
      dvalue=theScale*(static_cast<double*>(m_data[band])[index])+theOffset;
      break;
    default:
      std::string errorString="Error: data type not supported";
      throw(errorString);
      break;
    }
    value=static_cast<T>(dvalue);
  }
  else{
    //fetch raster band
    GDALRasterBand  *poBand;
    poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
    poBand->RasterIO(GF_Read,col,row,1,1,&value,1,1,getGDALDataType<T>(),0,0);
    if(m_scale.size()>band)
      value=static_cast<T>(value)*m_scale[band];
    if(m_offset.size()>band)
      value=static_cast<T>(value)+m_offset[band];
  }
}

template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band)
{
  assert(band<nrOfBand()+1);
  assert(minCol<nrOfCol());
  assert(minCol>=0);
  assert(maxCol<nrOfCol());
  assert(minCol<=maxCol);
  assert(row<nrOfRow());
  assert(row>=0);
  if(m_data.size()){
    if(row<m_begin[band]||row>=m_end[band]){
      readNewBlock(row,band);
    }
    if(buffer.size()!=maxCol-minCol+1)
      buffer.resize(maxCol-minCol+1);
    int index=(row-m_begin[band])*nrOfCol();
    int minindex=(index+minCol);
    int maxindex=(index+maxCol);
    double theScale=1;
    double theOffset=0;
    if(m_scale.size()>band||m_offset.size()>band){
      if(m_scale.size()>band)
	theScale=m_scale[band];
      if(m_offset.size()>band)
	theOffset=m_offset[band];
    }
    if(getGDALDataType<T>()==getDataType()){//no conversion needed
      buffer.assign(static_cast<T*>(m_data[band])+minindex,static_cast<T*>(m_data[band])+maxindex);
      typename std::vector<T>::iterator bufit=buffer.begin();
      while(bufit!=buffer.end()){
	double dvalue=theScale*(*bufit)+theOffset;
	*(bufit++)=static_cast<T>(dvalue);
      }
    }
    else{
      typename std::vector<T>::iterator bufit=buffer.begin();
      for(index=minindex;index<=maxindex;++index,++bufit){
        double dvalue=0;
        switch(getDataType()){
        case(GDT_Byte):
          dvalue=theScale*(static_cast<unsigned char*>(m_data[band])[index])+theOffset;
          break;
        case(GDT_Int16):
          dvalue=theScale*(static_cast<short*>(m_data[band])[index])+theOffset;
          break;
        case(GDT_UInt16):
          dvalue=theScale*(static_cast<unsigned short*>(m_data[band])[index])+theOffset;
          break;
        case(GDT_Int32):
          dvalue=theScale*(static_cast<int*>(m_data[band])[index])+theOffset;
          break;
        case(GDT_UInt32):
          dvalue=theScale*(static_cast<unsigned int*>(m_data[band])[index])+theOffset;
          break;
        case(GDT_Float32):
          dvalue=theScale*(static_cast<float*>(m_data[band])[index])+theOffset;
          break;
        case(GDT_Float64):
          dvalue=theScale*(static_cast<double*>(m_data[band])[index])+theOffset;
          break;
        default:
          std::string errorString="Error: data type not supported";
          throw(errorString);
	  break;
        }
        // double dvalue=theScale*(*(static_cast<double*>(m_data[band])+index))+theOffset;
        *(bufit)=static_cast<T>(dvalue);
      }
    }
  }
  else{
    //fetch raster band
    GDALRasterBand  *poBand;
    poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
    if(buffer.size()!=maxCol-minCol+1)
      buffer.resize(maxCol-minCol+1);
    poBand->RasterIO(GF_Read,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,getGDALDataType<T>(),0,0);
    if(m_scale.size()>band||m_offset.size()>band){
      double theScale=1;
      double theOffset=0;
      if(m_scale.size()>band)
	theScale=m_scale[band];
      if(m_offset.size()>band)
	theOffset=m_offset[band];
      for(int index=0;index<buffer.size();++index)
	buffer[index]=theScale*static_cast<double>(buffer[index])+theOffset;
    }
  }
}

//deprecated: there is a new interpolation argument from GDAL 2.0 (see http://www.gdal.org/structGDALRasterIOExtraArg.html)
template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, int minCol, int maxCol, double row, int band, RESAMPLE resample)
{
  std::vector<T> readBuffer_upper;
  std::vector<T> readBuffer_lower;
  if(buffer.size()!=maxCol-minCol+1)
    buffer.resize(maxCol-minCol+1);
  double upperRow=row-0.5;
  upperRow=static_cast<int>(upperRow);
  double lowerRow=row+0.5;
  lowerRow=static_cast<int>(lowerRow);
  switch(resample){
  case(BILINEAR):
    if(lowerRow>=nrOfRow())
      lowerRow=nrOfRow()-1;
    if(upperRow<0)
      upperRow=0;
    readData(readBuffer_upper,minCol,maxCol,static_cast<int>(upperRow),band);
    readData(readBuffer_lower,minCol,maxCol,static_cast<int>(lowerRow),band);
    //do interpolation in y
    for(int icol=0;icol<maxCol-minCol+1;++icol){
      buffer[icol]=(lowerRow-row+0.5)*readBuffer_upper[icol]+(1-lowerRow+row-0.5)*readBuffer_lower[icol];
    }
    break;
  default:
    readData(buffer,minCol,maxCol,static_cast<int>(row),band);
    break;
  }
}

template<typename T> void ImgReaderGdal::readDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  buffer2d.resize(maxRow-minRow+1);
  typename std::vector<T> buffer;
  readDataBlock(buffer,minCol,maxCol,minRow,maxRow,band);
  typename std::vector<T>::const_iterator startit=buffer.begin();
  typename std::vector<T>::const_iterator endit=startit;
  for(int irow=minRow;irow<=maxRow;++irow){
    buffer2d[irow-minRow].resize(maxCol-minCol+1);
    endit+=maxCol-minCol+1;
    buffer2d[irow-minRow].assign(startit,endit);
    startit+=maxCol-minCol+1;
  }
}
  
template<typename T> void ImgReaderGdal::readDataBlock(std::vector<T>& buffer, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band)
    theScale=m_scale[band];
  if(m_offset.size()>band)
    theOffset=m_offset[band];
  if(minCol>=nrOfCol() ||
     (minCol<0) ||
     (maxCol>=nrOfCol()) ||
     (minCol>maxCol) ||
     (minRow>=nrOfRow()) ||
     (minRow<0) ||
     (maxRow>=nrOfRow()) ||
     (minRow>maxRow)){
    std::string errorString="block not within image boundaries";
    throw(errorString);
  }
  if(buffer.size()!=(maxRow-minRow+1)*(maxCol-minCol+1))
    buffer.resize((maxRow-minRow+1)*(maxCol-minCol+1));
  if(m_data.size()){
    typename std::vector<T>::iterator bufit=buffer.begin();
    for(int irow=minRow;irow<=maxRow;++irow){
      if(irow<m_begin[band]||irow>=m_end[band])
	readNewBlock(irow,band);
      int index=(irow-m_begin[band])*nrOfCol();
      int minindex=(index+minCol);//*(GDALGetDataTypeSize(getDataType())>>3);
      int maxindex=(index+maxCol);//*(GDALGetDataTypeSize(getDataType())>>3);

      if(getGDALDataType<T>()==getDataType()){//no conversion needed
	//assign will replace current contents and modify its size accordingly
	buffer.assign(static_cast<T*>(m_data[band])+minindex,static_cast<T*>(m_data[band])+maxindex);
      }
      else{
	for(index=minindex;index<=maxindex;++index,++bufit){
	  double dvalue=0;
	  switch(getDataType()){
	  case(GDT_Byte):
	    dvalue=theScale*(static_cast<unsigned char*>(m_data[band])[index])+theOffset;
	    break;
	  case(GDT_Int16):
	    dvalue=theScale*(static_cast<short*>(m_data[band])[index])+theOffset;
	    break;
	  case(GDT_UInt16):
	    dvalue=theScale*(static_cast<unsigned short*>(m_data[band])[index])+theOffset;
	    break;
	  case(GDT_Int32):
	    dvalue=theScale*(static_cast<int*>(m_data[band])[index])+theOffset;
	    break;
	  case(GDT_UInt32):
	    dvalue=theScale*(static_cast<unsigned int*>(m_data[band])[index])+theOffset;
	    break;
	  case(GDT_Float32):
	    dvalue=theScale*(static_cast<float*>(m_data[band])[index])+theOffset;
	    break;
	  case(GDT_Float64):
	    dvalue=theScale*(static_cast<double*>(m_data[band])[index])+theOffset;
	    break;
	  default:
	    std::string errorString="Error: data type not supported";
	    throw(errorString);
	    break;
	  }
	  *(bufit)=static_cast<T>(dvalue);
	}//for index
      }//else
      if(getGDALDataType<T>()==getDataType()){
	if(m_scale.size()>band||m_offset.size()>band){
	  for(bufit=buffer.begin();bufit!=buffer.end();++bufit){
	    double dvalue=theScale*(*bufit)+theOffset;
	    *(bufit)=static_cast<T>(dvalue);
	  }
	}
      }
    }
  }
  else{
    //fetch raster band
    GDALRasterBand  *poBand;
    assert(band<nrOfBand()+1);
    poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
    poBand->RasterIO(GF_Read,minCol,minRow,maxCol-minCol+1,maxRow-minRow+1,&(buffer[0]),(maxCol-minCol+1),(maxRow-minRow+1),getGDALDataType<T>(),0,0);
    if(m_scale.size()>band||m_offset.size()>band){
      for(int index=0;index<buffer.size();++index)
	buffer[index]=theScale*buffer[index]+theOffset;
    }
  }
}

template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, int row, int band)
{
  readData(buffer,0,nrOfCol()-1,row,band);
}

template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, double row, int band, RESAMPLE resample)
{
  readData(buffer,0,nrOfCol()-1,row,band,resample);
}

#endif // _IMGREADERGDAL_H_

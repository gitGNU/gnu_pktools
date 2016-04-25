/**********************************************************************
ImgReaderMem.h: class to read raster files using GDAL API library
Copyright (C) 2008-2016 Pieter Kempeneers

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
#ifndef _IMGREADERMEM_H_
#define _IMGREADERMEM_H_

#include "ImgReaderGdal.h"
#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "base/Vector2d.h"

//--------------------------------------------------------------------------
class ImgReaderMem : public ImgReaderGdal
{
public:
  ImgReaderMem(void) : m_deletePointer(true) {};
  //memory (in MB) indicates maximum memory read in memory for this image
  ImgReaderMem(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly, unsigned long int memory=0) : ImgReaderGdal(filename, readMode), m_deletePointer(true) {initMem(memory);};
  ImgReaderMem(void* dataPointer, unsigned int ncol, unsigned int nrow, unsigned short nband, const GDALDataType& dataType){open(dataPointer,ncol,nrow,nband,dataType);};
  //memory (in MB) indicates maximum memory read in memory for this image
  void open(void* dataPointer, unsigned int ncol, unsigned int nrow, unsigned short nband, const GDALDataType& dataType);
  void open(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly, unsigned long int memory=0){ImgReaderGdal::open(filename, readMode);initMem(memory);};
  void setMemory(unsigned long int memory=0){initMem(memory);};
  ~ImgReaderMem(void){if(m_deletePointer);for(int iband=0;iband<m_nband;++iband) free(m_data[iband]);};
  template<typename T1> void readData(T1& value, int col, int row, int band=0);
  template<typename T1> void readData(std::vector<T1>& buffer, int minCol, int maxCol, int row, int band=0);
  template<typename T1> void readData(std::vector<T1>& buffer, int minCol, int maxCol, double row, int band=0, RESAMPLE resample=NEAR);
  template<typename T1> void readDataBlock(Vector2d<T1>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  template<typename T1> void readDataBlock(std::vector<T1>& buffer, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  template<typename T1> void readData(std::vector<T1>& buffer, int row, int band=0);
  template<typename T1> void readData(std::vector<T1>& buffer, double row, int band=0, RESAMPLE resample=NEAR);

protected:
  unsigned int m_blockSize;
  std::vector<void *> m_data;
  std::vector<unsigned int> m_begin;//first line that has been read
  std::vector<unsigned int> m_end;//beyond last line read
private:
  void initMem(unsigned long int memory);
  bool m_deletePointer;
  bool readNewBlock(int row, int band);
};

//not tested yet!!!
//open image in memory (passing pointer to allocated memory). This will allow in place image processing in memory (streaming)
//what about projection and geotransform data (must be set from outside!)
void ImgReaderMem::open(void* dataPointer, unsigned int ncol, unsigned int nrow, unsigned short nband, const GDALDataType& dataType)
{
  m_deletePointer=false;//we are not the owner
  m_nband=nband;
  m_ncol=ncol;
  m_nrow=nrow;
  m_data.resize(nband);
  m_begin.resize(nband);
  m_end.resize(nband);
  for(int iband=0;iband<nband;++iband){
    m_data[iband]=dataPointer+iband*m_ncol*m_nrow*(GDALGetDataTypeSize(getDataType())>>3);
    m_begin[iband]=0;
    m_blockSize=nrow;//memory contains entire image
    m_end[iband]=nrow;//and has been read already
  }
}

void ImgReaderMem::initMem(unsigned long int memory)
{
  m_blockSize=static_cast<unsigned int>(memory*1000000/nrOfBand()/nrOfCol());
  if(m_blockSize<1)
    m_blockSize=1;
  if(m_blockSize>nrOfRow())
    m_blockSize=nrOfRow();
  m_data.resize(nrOfBand());
  m_begin.resize(nrOfBand());
  m_end.resize(nrOfBand());
  for(int iband=0;iband<m_nband;++iband){
    m_data[iband]=(void *) CPLMalloc((GDALGetDataTypeSize(getDataType())>>3)*nrOfCol()*m_blockSize);
    m_begin[iband]=0;
    m_end[iband]=0;
  }
  // readNewBlock(0);
}

bool ImgReaderMem::readNewBlock(int row, int band)
{
  if(m_end[band]<m_blockSize)//first time
    m_end[band]=m_blockSize;
  while(row>=m_end[band]&&m_begin[band]<nrOfRow()){
    m_begin[band]+=m_blockSize;
    m_end[band]=m_begin[band]+m_blockSize;
  }
  if(m_end[band]>nrOfRow())
    m_end[band]=nrOfRow();
  for(int iband=0;iband<m_nband;++iband){
    //fetch raster band
    GDALRasterBand  *poBand;
    assert(iband<nrOfBand()+1);
    poBand = m_gds->GetRasterBand(iband+1);//GDAL uses 1 based index
    poBand->RasterIO(GF_Read,0,m_begin[iband],nrOfCol(),m_end[iband]-m_begin[iband],m_data[iband],nrOfCol(),m_end[iband]-m_begin[iband],getDataType(),0,0);
  }
  return true;//new block was read
}

template<typename T1> void ImgReaderMem::readData(T1& value, int col, int row, int band)
{
  //only support random access reading if entire image is in memory for performance reasons
  if(m_blockSize!=nrOfRow()){
    std::ostringstream s;
    s << "Error: increase memory to support random access reading (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
    throw(s.str());
  }
  if(row<m_begin[band]||row>=m_end[band])
    readNewBlock(row,band);
  assert(band<nrOfBand()+1);
  assert(col<nrOfCol());
  assert(col>=0);
  assert(row<nrOfRow());
  assert(row>=0);
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
  value=static_cast<T1>(dvalue);
}

template<typename T1> void ImgReaderMem::readData(std::vector<T1>& buffer, int minCol, int maxCol, int row, int band)
{
  if(row<m_begin[band]||row>=m_end[band]){
    readNewBlock(row,band);
  }
  assert(band<nrOfBand()+1);
  assert(minCol<nrOfCol());
  assert(minCol>=0);
  assert(maxCol<nrOfCol());
  assert(minCol<=maxCol);
  assert(row<nrOfRow());
  assert(row>=0);
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
  if(getGDALDataType<T1>()==getDataType()){//no conversion needed
    buffer.assign(static_cast<T1*>(m_data[band])+minindex,static_cast<T1*>(m_data[band])+maxindex);
    typename std::vector<T1>::iterator bufit=buffer.begin();
    while(bufit!=buffer.end()){
      double dvalue=theScale*(*bufit)+theOffset;
      *(bufit++)=static_cast<T1>(dvalue);
    }
  }
  else{
    typename std::vector<T1>::iterator bufit=buffer.begin();
    for(index=minindex;index<maxindex;++index,++bufit){
      if(m_scale.size()>band||m_offset.size()>band){
        double theScale=1;
        double theOffset=0;
        if(m_scale.size()>band)
          theScale=m_scale[band];
        if(m_offset.size()>band)
          theOffset=m_offset[band];
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
        *(bufit)=static_cast<T1>(dvalue);
      }
    }
  }
}

//deprecated: there is a new interpolation argument from GDAL 2.0 (see http://www.gdal.org/structGDALRasterIOExtraArg.html)
template<typename T1> void ImgReaderMem::readData(std::vector<T1>& buffer, int minCol, int maxCol, double row, int band, RESAMPLE resample)
{
  std::vector<T1> readBuffer_upper;
  std::vector<T1> readBuffer_lower;
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

template<typename T1> void ImgReaderMem::readDataBlock(Vector2d<T1>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  buffer2d.resize(maxRow-minRow+1);
  typename std::vector<T1> buffer;
  readDataBlock(buffer,minCol,maxCol,minRow,maxRow,band);
  typename std::vector<T1>::const_iterator startit=buffer.begin();
  typename std::vector<T1>::const_iterator endit=startit;
  for(int irow=minRow;irow<=maxRow;++irow){
    buffer2d[irow-minRow].resize(maxCol-minCol+1);
    endit+=maxCol-minCol+1;
    buffer2d[irow-minRow].assign(startit,endit);
    startit+=maxCol-minCol+1;
  }
}

template<typename T1> void ImgReaderMem::readDataBlock(std::vector<T1>& buffer, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band)
    theScale=m_scale[band];
  if(m_offset.size()>band)
    theOffset=m_offset[band];
  assert(band<nrOfBand()+1);
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
  typename std::vector<T1>::iterator bufit=buffer.begin();
  for(int irow=minRow;irow<=maxRow;++irow){
    if(irow<m_begin[band]||irow>=m_end[band])
      readNewBlock(irow,band);
    int index=(irow-m_begin[band])*nrOfCol();
    int minindex=(index+minCol);//*(GDALGetDataTypeSize(getDataType())>>3);
    int maxindex=(index+maxCol);//*(GDALGetDataTypeSize(getDataType())>>3);

    if(getGDALDataType<T1>()==getDataType()){//no conversion needed
      //assign will replace current contents and modify its size accordingly
      buffer.assign(static_cast<T1*>(m_data[band])+minindex,static_cast<T1*>(m_data[band])+maxindex);
    }
    else{
      for(index=minindex;index<maxindex;++index,++bufit){
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
        *(bufit)=static_cast<T1>(dvalue);
      }//for index
    }//else
    if(getGDALDataType<T1>()==getDataType()){
      if(m_scale.size()>band||m_offset.size()>band){
	for(bufit=buffer.begin();bufit!=buffer.end();++bufit){
	  double dvalue=theScale*(*bufit)+theOffset;
	  *(bufit)=static_cast<T1>(dvalue);
	}
      }
    }
  }
}
  
template<typename T1> void ImgReaderMem::readData(std::vector<T1>& buffer, int row, int band)
{
  readData(buffer,0,nrOfCol()-1,row,band);
}

template<typename T1> void ImgReaderMem::readData(std::vector<T1>& buffer, double row, int band, RESAMPLE resample)
{
  readData(buffer,0,nrOfCol()-1,row,band,resample);
}


#endif // _IMGREADERMEM_H_

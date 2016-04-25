/**********************************************************************
ImgWriterMem.h: class to write raster files using GDAL API library
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
#ifndef _IMGWRITERMEM_H_
#define _IMGWRITERMEM_H_

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "ImgRasterGdal.h"
#include "ImgReaderGdal.h"
#include "ImgWriterGdal.h"
#include "ImgReaderOgr.h"

//--------------------------------------------------------------------------
class ImgWriterMem : public virtual ImgWriterGdal
{
public:
  ImgWriterMem(void) : m_deletePointer(true) {};
  //memory (in MB) indicates maximum memory read in memory for this image
  ImgWriterMem(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>(), unsigned long int memory=0) : ImgWriterGdal(filename, imgSrc, options), m_deletePointer(true) {initMem(memory);};
  ImgWriterMem(void* dataPointer, const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>(), unsigned long int memory=0){open(dataPointer,filename, imgSrc, options, memory);};
  //memory (in MB) indicates maximum memory read in memory for this image
  ImgWriterMem(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>()){open(dataPointer, filename, ncol, nrow, nband, dataType, imageType, options);};
  ImgWriterMem(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>(), unsigned long int memory=0) : m_deletePointer(true) {ImgWriterGdal::open(filename, ncol, nrow, nband, dataType, imageType, options);initMem(memory);};
  void open(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>(), unsigned long int memory=0){ImgWriterGdal::open(filename, ncol, nrow, nband, dataType, imageType, options);initMem(memory);};
  void open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>(), unsigned long int memory=0){ImgWriterGdal::open(filename, imgSrc, options),initMem(memory);};
  void open(void* dataPointer, const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>(), unsigned long int memory=0){open(dataPointer,filename, imgSrc, options, memory);};
  ~ImgWriterMem(void){
    if(m_deletePointer){
      for(int iband=0;iband<m_nband;++iband) 
        free(m_data[iband]);
    }
  };
  void close(void){
    for(int iband=0;iband<nrOfBand();++iband) 
      writeNewBlock(nrOfRow(),iband);
    ImgRasterGdal::close();
  };
  template<typename T> bool writeData(T& value, int col, int row, int band=0);
  template<typename T> bool writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  template<typename T> bool writeData(std::vector<T>& buffer, int row, int band=0);
  bool writeData(void* pdata, const GDALDataType& dataType, int band=0);
  template<typename T> bool writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);

protected:
  unsigned int m_blockSize;
  std::vector<void *> m_data;
  std::vector<unsigned int> m_begin;//[band] first line in block that will be written next
  std::vector<unsigned int> m_end;//[band] beyond last line in block that will be written next
private:
  void initMem(unsigned long int memory);
  bool m_deletePointer;
  bool writeNewBlock(int row, int band);
};

//not tested yet!!!
//open image in memory (passing pointer to allocated memory). This will allow in place image processing in memory (streaming)
void ImgWriterMem::open(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options){
  ImgWriterGdal::open(filename, ncol, nrow, nband, dataType, imageType, options);
  m_deletePointer=false;//we are not the owner
  m_data.resize(nband);
  m_begin.resize(nband);
  m_end.resize(nband);
  for(int iband=0;iband<nband;++iband){
    m_data[iband]=dataPointer+iband*ncol*nrow*(GDALGetDataTypeSize(getDataType())>>3);
    m_begin[iband]=0;
    m_end[iband]=nrow;
  }
  m_blockSize=nrow;//memory contains entire image and has been read already
}


void ImgWriterMem::initMem(unsigned long int memory)
{
  //we will only write after processing a block
  m_blockSize=static_cast<unsigned int>(memory*1000000/nrOfBand()/nrOfCol());
  if(m_blockSize<1)
    m_blockSize=1;
  if(m_blockSize>nrOfRow())
    m_blockSize=nrOfRow();
  //allocate block of memory
  m_data.resize(m_nband);
  m_begin.resize(m_nband);
  m_end.resize(m_nband);
  for(int iband=0;iband<m_nband;++iband){
    m_data[iband]=(void *) CPLMalloc((GDALGetDataTypeSize(getDataType())>>3)*nrOfCol()*m_blockSize);
    m_begin[iband]=0;
    m_end[iband]=m_blockSize;
  }
}

bool ImgWriterMem::writeNewBlock(int row, int band)
{
  //assert(row==m_end)
  if(m_end[band]>nrOfRow())
    m_end[band]=nrOfRow();
  //fetch raster band
  GDALRasterBand  *poBand;
  assert(band<nrOfBand()+1);
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  poBand->RasterIO(GF_Write,0,m_begin[band],nrOfCol(),m_end[band]-m_begin[band],m_data[band],nrOfCol(),m_end[band]-m_begin[band],getDataType(),0,0);
  m_begin[band]+=m_blockSize;//m_begin points to first line in block that will be written next
  m_end[band]=m_begin[band]+m_blockSize;//m_end points to last line in block that will be written next
  return true;//new block was read
}

template<typename T> bool ImgWriterMem::writeData(T& value, int col, int row, int band)
{
  //only support random access writing if entire image is in memory
  //if m_blockSize!=nrOfRow() we risk not to write all data if not written in sequence (e.g. line per line)
  if(m_blockSize!=nrOfRow()){
    std::ostringstream s;
    s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
    throw(s.str());
  }
  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "Error: band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  if(col>=nrOfCol()){
    std::ostringstream s;
    s << "Error: col (" << col << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(col<0){
    std::ostringstream s;
    s << "Error: col (" << col << ") is negative";
    throw(s.str());
  }
  if(row>=nrOfRow()){
    std::ostringstream s;
    s << "Error: row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  if(row<0){
    std::ostringstream s;
    s << "Error: row (" << row << ") is negative";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  int index=(row-m_begin[band])*nrOfCol()+col;
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band||m_offset.size()>band){
    if(m_scale.size()>band)
      theScale=m_scale[band];
    if(m_offset.size()>band)
      theOffset=m_offset[band];
  }
  double dvalue=theScale*value+theOffset;
  switch(getDataType()){
  case(GDT_Byte):
    *(static_cast<unsigned char*>((m_data[band])+index))=static_cast<unsigned char>(dvalue);
    break;
  case(GDT_Int16):
    *(static_cast<short*>((m_data[band])+index))=static_cast<short>(dvalue);
    break;
  case(GDT_UInt16):
    *(static_cast<unsigned short*>((m_data[band])+index))=static_cast<unsigned short>(dvalue);
    break;
  case(GDT_Int32):
    *(static_cast<int*>((m_data[band])+index))=static_cast<int>(dvalue);
    break;
  case(GDT_UInt32):
    *(static_cast<unsigned int*>((m_data[band])+index))=static_cast<unsigned int>(dvalue);
    break;
  case(GDT_Float32):
    *(static_cast<float*>((m_data[band])+index))=static_cast<float>(dvalue);
    break;
  case(GDT_Float64):
    *(static_cast<double*>((m_data[band])+index))=static_cast<double>(dvalue);
    break;
  default:
    std::string errorString="Error: data type not supported";
    throw(errorString);
    break;
  }
  return true;
}

//todo: make buffer const?
template<typename T> bool ImgWriterMem::writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band)
{
  if(buffer.size()!=maxCol-minCol+1){
    std::string errorstring="invalid buffer size";
    throw(errorstring);
  }
  if(row>=nrOfRow()){
    std::ostringstream s;
    s << "row (" << row << ") exceeds nrOfRow (" << nrOfRow() << ")";
    throw(s.str());
  }
  if(row<0){
    std::ostringstream s;
    s << "row (" << row << ") is negative";
    throw(s.str());
  }
  if(row<m_begin[band]){
    std::ostringstream s;
    s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
    throw(s.str());
  }
  if(row>=m_end[band]){
    if(row>=m_end[band]+m_blockSize){
      std::ostringstream s;
      s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
      throw(s.str());
    }
    else
      writeNewBlock(row,band);
  }
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  if(minCol>=nrOfCol()){
    std::ostringstream s;
    s << "minCol (" << minCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(minCol<0){
    std::ostringstream s;
    s << "mincol (" << minCol << ") is negative";
    throw(s.str());
  }
  if(minCol>0){
    std::ostringstream s;
    s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
    throw(s.str());
  }
  if(maxCol>=nrOfCol()){
    std::ostringstream s;
    s << "maxCol (" << maxCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(maxCol<minCol){
    std::ostringstream s;
    s << "maxCol (" << maxCol << ") is less than minCol (" << minCol << ")";
    throw(s.str());
  }
  //fetch raster band
  GDALRasterBand  *poBand;
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  int index=(row-m_begin[band])*nrOfCol();
  int minindex=(index+minCol);
  int maxindex=(index+maxCol);
  typename std::vector<T>::iterator bufit=buffer.begin();
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band)
    theScale=m_scale[band];
  if(m_offset.size()>band)
    theOffset=m_offset[band];
  for(index=minindex;index<maxindex;++index,++bufit){
    double dvalue=theScale*(*(bufit))+theOffset;
    switch(getDataType()){
    case(GDT_Byte):
      //test
      if(row<10&&index-minindex<10)
	;
	/* std::cout << buffer[0] << " "; */
	/* std::cout << dvalue << " "; */
      if(row<10&&index==maxindex-1)
	;/* std::cout << std::endl; */
      *(static_cast<unsigned char*>((m_data[band])+index))=static_cast<unsigned char>(dvalue);
      break;
    case(GDT_Int16):
      *(static_cast<short*>((m_data[band])+index))=static_cast<short>(dvalue);
      break;
    case(GDT_UInt16):
      *(static_cast<unsigned short*>((m_data[band])+index))=static_cast<unsigned short>(dvalue);
      break;
    case(GDT_Int32):
      *(static_cast<int*>((m_data[band])+index))=static_cast<int>(dvalue);
      break;
    case(GDT_UInt32):
      *(static_cast<unsigned int*>((m_data[band])+index))=static_cast<unsigned int>(dvalue);
      break;
    case(GDT_Float32):
      *(static_cast<float*>((m_data[band])+index))=static_cast<float>(dvalue);
      break;
    case(GDT_Float64):
      *(static_cast<double*>((m_data[band])+index))=static_cast<double>(dvalue);
      break;
    default:
      std::string errorString="Error: data type not supported";
      throw(errorString);
      break;
    }
  }
}

//todo: make buffer const?
template<typename T> bool ImgWriterMem::writeData(std::vector<T>& buffer, int row, int band)
{
  return writeData(buffer,0,nrOfCol()-1,row,band);
}

//todo: make buffer2d const?
template<typename T> bool ImgWriterMem::writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band)
    theScale=m_scale[band];
  if(m_offset.size()>band)
    theOffset=m_offset[band];
  if(buffer2d.size()!=maxRow-minRow+1){
    std::string errorstring="invalid buffer size";
    throw(errorstring);
  }
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  if(minCol>=nrOfCol()){
    std::ostringstream s;
    s << "minCol (" << minCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(minCol<0){
    std::ostringstream s;
    s << "mincol (" << minCol << ") is negative";
    throw(s.str());
  }
  if(minCol>0){
    std::ostringstream s;
    s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
    throw(s.str());
  }
  if(maxCol>=nrOfCol()){
    std::ostringstream s;
    s << "maxCol (" << maxCol << ") exceeds nrOfCol (" << nrOfCol() << ")";
    throw(s.str());
  }
  if(maxCol<minCol){
    std::ostringstream s;
    s << "maxCol (" << maxCol << ") is less than minCol (" << minCol << ")";
    throw(s.str());
  }

  //fetch raster band
  GDALRasterBand  *poBand;
  if(band>=nrOfBand()+1){
    std::ostringstream s;
    s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
    throw(s.str());
  }
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  for(int irow=minRow;irow<=maxRow;++irow){
    if(irow>=nrOfRow()){
      std::ostringstream s;
      s << "row (" << irow << ") exceeds nrOfRow (" << nrOfRow() << ")";
      throw(s.str());
    }
    if(irow<0){
      std::ostringstream s;
      s << "row (" << irow << ") is negative";
      throw(s.str());
    }
    if(irow<m_begin[band]){
      std::ostringstream s;
      s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
      throw(s.str());
    }
    if(irow>=m_end[band]){
      if(irow>=m_end[band]+m_blockSize){
        std::ostringstream s;
        s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
        throw(s.str());
      }
      else
        writeNewBlock(irow,band);
    }
    int index=(irow-m_begin[band])*nrOfCol();
    int minindex=index+minCol;
    int maxindex=index+maxCol;
    typename std::vector<T>::iterator bufit=buffer2d[irow-minRow].begin();
    for(index=minindex;index<maxindex;++index,++bufit){
      double dvalue=theScale*(*(bufit))+theOffset;
      switch(getDataType()){
      case(GDT_Byte):
	*(static_cast<unsigned char*>((m_data[band])+index))=static_cast<unsigned char>(dvalue);
	break;
      case(GDT_Int16):
	*(static_cast<short*>((m_data[band])+index))=static_cast<short>(dvalue);
	break;
      case(GDT_UInt16):
	*(static_cast<unsigned short*>((m_data[band])+index))=static_cast<unsigned short>(dvalue);
	break;
      case(GDT_Int32):
	*(static_cast<int*>((m_data[band])+index))=static_cast<int>(dvalue);
	break;
      case(GDT_UInt32):
	*(static_cast<unsigned int*>((m_data[band])+index))=static_cast<unsigned int>(dvalue);
	break;
      case(GDT_Float32):
	*(static_cast<float*>((m_data[band])+index))=static_cast<float>(dvalue);
	break;
      case(GDT_Float64):
	*(static_cast<double*>((m_data[band])+index))=static_cast<double>(dvalue);
	break;
      default:
	std::string errorString="Error: data type not supported";
	throw(errorString);
	break;
      }
    }
  }
  return true;
}

#endif // _IMGWRITERMEM_H_

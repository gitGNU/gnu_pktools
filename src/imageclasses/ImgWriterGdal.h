/**********************************************************************
ImgWriterGdal.h: class to write raster files using GDAL API library
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
#ifndef _IMGWRITERGDAL_H_
#define _IMGWRITERGDAL_H_

#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "ImgRasterGdal.h"
#include "ImgReaderGdal.h"
#include "ImgReaderOgr.h"

//--------------------------------------------------------------------------
class ImgWriterGdal : public virtual ImgRasterGdal
{
public:
  ImgWriterGdal(void);
  ~ImgWriterGdal(void);
  ImgWriterGdal(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, imgSrc, options);};
  ImgWriterGdal(const std::string& filename, const ImgReaderGdal& imgSrc, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, imgSrc, memory, options);};
  ImgWriterGdal(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, ncol, nrow, nband, dataType, imageType, options);};
  ImgWriterGdal(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, ncol, nrow, nband, dataType, imageType, options);};
  ImgWriterGdal(void* dataPointer, const std::string& filename, const ImgReaderGdal& imgSrc, unsigned long int memory=0, const std::vector<std::string>& options=std::vector<std::string>()){open(dataPointer,filename, imgSrc, options);};
  ImgWriterGdal(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>()){open(dataPointer, filename, ncol, nrow, nband, dataType, imageType, options);};
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
  // void open(const std::string& filename);//not needed?
  void open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>());
  void open(const std::string& filename, const ImgReaderGdal& imgSrc, unsigned int memory, const std::vector<std::string>& options=std::vector<std::string>());
  void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory, const std::vector<std::string>& options=std::vector<std::string>());
  void open(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  void open(void* dataPointer, const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>());
  
  void close(void);//definition in ImgWritergdal.cc

  void copyGeoTransform(const ImgReaderGdal& imgSrc);
  void setProjection(const std::string& projection);
  std::string setProjectionProj4(const std::string& projection);

  void setImageDescription(const std::string& imageDescription){m_gds->SetMetadataItem( "TIFFTAG_IMAGEDESCRIPTION",imageDescription.c_str());};
  void setGeoTransform(double* gt);

  template<typename T> bool writeData(T& value, int col, int row, int band=0);
  template<typename T> bool writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  template<typename T> bool writeData(std::vector<T>& buffer, int row, int band=0);
  bool writeData(void* pdata, const GDALDataType& dataType, int band=0);
  template<typename T> bool writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  // std::string getInterleave(){return m_interleave;};
  // std::string getCompression(){return m_compression;};
  void setColorTable(const std::string& filename, int band=0);
  void setColorTable(GDALColorTable* colorTable, int band=0);
  void setMetadata(char** metadata);
  void rasterizeOgr(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues=std::vector<double>(), const std::vector<std::string>& layernames=std::vector<std::string>());

protected:
  virtual void setCodec(const GDALDataType& dataType, const std::string& imageType);
  virtual void setCodec(const ImgReaderGdal& ImgSrc);
  std::vector<double> m_scale;
  std::vector<double> m_offset;

  std::vector<std::string> m_options;
  unsigned int m_blockSize;
  std::vector<void *> m_data;
  std::vector<unsigned int> m_begin;//[band] first line in block that will be written next
  std::vector<unsigned int> m_end;//[band] beyond last line in block that will be written next
private:
  void initMem(unsigned long int memory);
  bool m_deletePointer;
  bool writeNewBlock(int row, int band);
};

template<typename T> bool ImgWriterGdal::writeData(T& value, int col, int row, int band)
{
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
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band||m_offset.size()>band){
    if(m_scale.size()>band)
      theScale=m_scale[band];
    if(m_offset.size()>band)
      theOffset=m_offset[band];
  }
  if(m_data.size()){
    //only support random access writing if entire image is in memory
    //if m_blockSize!=nrOfRow() we risk not to write all data if not written in sequence (e.g. line per line)
    if(m_blockSize!=nrOfRow()){
      std::ostringstream s;
      s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
      throw(s.str());
    }
    int index=(row-m_begin[band])*nrOfCol()+col;
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
  }
  else{
    //fetch raster band
    GDALRasterBand  *poBand;
    T dvalue=theScale*value+theOffset;
    poBand->RasterIO(GF_Write,col,row,1,1,&dvalue,1,1,getGDALDataType<T>(),0,0);
  }
  return true;
}

//todo: make buffer const?
template<typename T> bool ImgWriterGdal::writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band)
{
  if(buffer.size()!=maxCol-minCol+1){
    std::string errorstring="invalid size of buffer";
    throw(errorstring);
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
  if(m_data.size()){
    if(minCol>0){
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
    for(index=minindex;index<=maxindex;++index,++bufit){
      double dvalue=theScale*(*(bufit))+theOffset;
      switch(getDataType()){
      case(GDT_Byte):
	*(static_cast<unsigned char*>((m_data[band])+index))=static_cast<unsigned char>(dvalue);
	//test
	/* if(band==0&&row==511){ */
	/*   std::cout << dvalue << " "; */
	/* if(maxindex-minindex+1==512) */
	/*   std::cout << std::endl; */
	/* } */
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
  else{
    //todo: scaling and offset!
    //fetch raster band
    GDALRasterBand  *poBand;
    if(band>=nrOfBand()+1){
      std::ostringstream s;
      s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
      throw(s.str());
    }
    poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
    poBand->RasterIO(GF_Write,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,getGDALDataType<T>(),0,0);
  }
  return true;
}

//todo: make buffer const?
template<typename T> bool ImgWriterGdal::writeData(std::vector<T>& buffer, int row, int band)
{
  return writeData(buffer,0,nrOfCol()-1,row,band);
}

//todo: make buffer2d const?
template<typename T> bool ImgWriterGdal::writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band)
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
  if(m_data.size()){
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
      for(index=minindex;index<=maxindex;++index,++bufit){
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
  }
  else{
    //todo: apply scaling and offset!
    typename std::vector<T> buffer((maxRow-minRow+1)*(maxCol-minCol+1));
    //fetch raster band
    GDALRasterBand  *poBand;
    // typename std::vector<T>::iterator startit=buffer.begin();
    for(int irow=minRow;irow<=maxRow;++irow){
      buffer.insert(buffer.begin()+(maxCol-minCol+1)*(irow-minRow),buffer2d[irow-minRow].begin(),buffer2d[irow-minRow].end());
    }
    poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
    poBand->RasterIO(GF_Write,minCol,minRow,maxCol-minCol+1,maxRow-minRow+1,&(buffer[0]),(maxCol-minCol+1),(maxRow-minRow+1),getGDALDataType<T>(),0,0);
  }
  return true;
}

#endif // _IMGWRITERGDAL_H_


//     adfGeoTransform[0] /* top left x */
//     adfGeoTransform[1] /* w-e pixel resolution */
//     adfGeoTransform[2] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[3] /* top left y */
//     adfGeoTransform[4] /* rotation, 0 if image is "north up" */
//     adfGeoTransform[5] /* n-s pixel resolution */

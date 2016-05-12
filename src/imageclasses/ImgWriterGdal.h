/**********************************************************************
ImgWriterGdal.h: class to write raster files using GDAL API library
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

/** Class to write a raster dataset in a format supported by GDAL. Currently only those formats where the drivers support the Create method can be written. Data are cached in memory for a number of rows (if memory>0) before written to file.

   This class inherits from ImgRasterGdal, a general raster class to store e.g., filename, number of columns, rows and bands of the dataset. 

   If memory is set (in MB) to 0 (default), the raster is written line by line directly from file. A scale and offset can be set when writing the raster data values. The scaling and offset are applied on a per band basis. 

   For random access writing (not in sequential order line by line), set memory to 0 or a value sufficiently large to write the entire image to memory.
 **/

class ImgWriterGdal : public virtual ImgRasterGdal
{
public:
  ///default constructor. Image needs to be opened later with one of the open methods.
  ImgWriterGdal(void);
  ///constructor opening an image for writing, copying image attributes from a source image. Image is directly writen to file. Use the constructor with memory>0 to support caching
  ImgWriterGdal(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, imgSrc, options);};
  ///constructor opening an image for writing, copying image attributes from a source image. Caching is supported when memory>0
  ImgWriterGdal(const std::string& filename, const ImgReaderGdal& imgSrc, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, imgSrc, memory, options);};
  ///constructor opening an image for writing in memory, copying image attributes from a source image.
  ImgWriterGdal(const ImgReaderGdal& imgSrc){open(imgSrc);};
  ///constructor opening an image for writing, defining all image attributes. Image is directly written to file. Use the constructor with memory>0 to support caching
  ImgWriterGdal(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, ncol, nrow, nband, dataType, imageType, options);};
  ///constructor opening an image for writing, defining all image attributes. Caching is supported when memory>0
  ImgWriterGdal(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>()){open(filename, ncol, nrow, nband, dataType, imageType, options);};
  ///constructor opening an image for writing in memory, defining all image attributes
  ImgWriterGdal(int ncol, int nrow, int nband, const GDALDataType& dataType){open(ncol, nrow, nband, dataType);};
  ///constructor opening an image for writing using an external data pointer (not tested yet)
  ImgWriterGdal(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType){open(dataPointer, filename, ncol, nrow, nband, dataType);};
  ///constructor opening an image for writing in memory using an external data pointer (not tested yet)
  ImgWriterGdal(void* dataPointer, int ncol, int nrow, int nband, const GDALDataType& dataType){open(dataPointer, ncol, nrow, nband, dataType);};
  ///destructor
  ~ImgWriterGdal(void);

  ///Open an image for writing, copying image attributes from a source image. Image is directly written to file. Use the constructor with memory>0 to support caching
  void open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing, copying image attributes from a source image. Caching is supported when memory>0
  void open(const std::string& filename, const ImgReaderGdal& imgSrc, unsigned int memory, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing, defining all image attributes. Image is directly written to file. Use the constructor with memory>0 to support caching
  void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing, defining all image attributes. Caching is supported when memory>0
  void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory, const std::vector<std::string>& options=std::vector<std::string>());
  void open(int ncol, int nrow, int nband, const GDALDataType& dataType);
  ///Open an image for writing, defining all image attributes. Caching is supported when memory>0
  ///Open an image for writing in memory, copying image attributes from a source image.
  void open(const ImgReaderGdal& imgSrc);
  ///Open an image for writing using an external data pointer (not tested yet)
  void open(void* dataPointer, const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType);
  ///Open an image for writing in memory using an external data pointer (not tested yet)
  void open(void* dataPointer, int ncol, int nrow, int nband, const GDALDataType& dataType);
  ///Close the raster dataset
  void close(void);//definition in ImgWritergdal.cc
  ///Set the image description (only for GeoTiff format: TIFFTAG_IMAGEDESCRIPTION)
  void setImageDescription(const std::string& imageDescription){m_gds->SetMetadataItem( "TIFFTAG_IMAGEDESCRIPTION",imageDescription.c_str());};
  ///Copy geotransform information from another georeferenced image
  void copyGeoTransform(const ImgReaderGdal& imgSrc);

  ///Write a single pixel cell value at a specific column and row for a specific band (all indices start counting from 0)
  template<typename T> bool writeData(T& value, int col, int row, int band=0);
  ///Write pixel cell values for a range of columns for a specific row and band (all indices start counting from 0)
  template<typename T> bool writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  ///Write pixel cell values for an entire row for a specific band (all indices start counting from 0)
  template<typename T> bool writeData(std::vector<T>& buffer, int row, int band=0);
  // deprecated? Write an entire image from memory to file
  // bool writeData(void* pdata, const GDALDataType& dataType, int band=0);
  ///Write pixel cell values for a range of columns and rows for a specific band (all indices start counting from 0). The buffer is a two dimensional vector (stl vector of stl vector) representing [row][col].
  template<typename T> bool writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  ///Write pixel cell values for the entire image from memory to file
  void setFile(const std::string& filename, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  ///Set the color table using an (ASCII) file with 5 columns (value R G B alpha)
  void setColorTable(const std::string& filename, int band=0);
  ///Set the color table using the GDAL class GDALColorTable
  void setColorTable(GDALColorTable* colorTable, int band=0);
  ///Set specific metadata (driver specific)
  void setMetadata(char** metadata);
  ///Rasterize an OGR vector dataset using the gdal algorithm "GDALRasterizeLayers"
  void rasterizeOgr(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues, const std::vector<std::string>& controlOptions=std::vector<std::string>(), const std::vector<std::string>& layernames=std::vector<std::string>());
  ///Rasterize an OGR vector dataset in memory using the gdal algorithm "GDALRasterizeLayersBuf"
  void rasterizeBuf(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues, const std::vector<std::string>& controlOptions=std::vector<std::string>(), const std::vector<std::string>& layernames=std::vector<std::string>());

protected:
  ///Register GDAL driver, setting the datatype, imagetype and some metadata
  virtual void setCodec(const std::string& imageType);
  ///Register GDAL driver, setting the datatype, imagetype and some metadata
  virtual void setCodec(const ImgReaderGdal& ImgSrc);
  std::vector<std::string> m_options;

private:
  ///Write new block from cache (defined by m_begin and m_end)
  bool writeNewBlock(int row, int band);
  ///Initialize the memory for read/write image in cache
  void initMem(unsigned long int memory);
};

/**
 * @param[in] value The cell value to write
 * @param[in] col The column number to write (counting starts from 0)
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
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
    if(m_blockSize!=nrOfRow()){
      std::ostringstream s;
      s << "Error: increase memory to support random access writing (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
      throw(s.str());
    }
    int index=(row-m_begin[band])*nrOfCol()+col;
    double dvalue=theScale*value+theOffset;
    switch(getDataType()){
    case(GDT_Byte):
      static_cast<unsigned char*>(m_data[band])[index]=static_cast<unsigned char>(dvalue);
      break;
    case(GDT_Int16):
      static_cast<short*>(m_data[band])[index]=static_cast<short>(dvalue);
      break;
    case(GDT_UInt16):
      static_cast<unsigned short*>(m_data[band])[index]=static_cast<unsigned short>(dvalue);
      break;
    case(GDT_Int32):
      static_cast<int*>(m_data[band])[index]=static_cast<int>(dvalue);
      break;
    case(GDT_UInt32):
      static_cast<unsigned int*>(m_data[band])[index]=static_cast<unsigned int>(dvalue);
      break;
    case(GDT_Float32):
      static_cast<float*>(m_data[band])[index]=static_cast<float>(dvalue);
      break;
    case(GDT_Float64):
      static_cast<double*>(m_data[band])[index]=static_cast<double>(dvalue);
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

/**
 * @param[in] buffer The vector with all cell values to write
 * @param[in] minCol First column from where to start writing (counting starts from 0)
 * @param[in] maxCol Last column that must be written (counting starts from 0)
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
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
      else if(m_filename.size())
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
	static_cast<unsigned char*>(m_data[band])[index]=static_cast<unsigned char>(dvalue);
	break;
      case(GDT_Int16):
	static_cast<short*>(m_data[band])[index]=static_cast<short>(dvalue);
	break;
      case(GDT_UInt16):
	static_cast<unsigned short*>(m_data[band])[index]=static_cast<unsigned short>(dvalue);
	break;
      case(GDT_Int32):
	static_cast<int*>(m_data[band])[index]=static_cast<int>(dvalue);
	break;
      case(GDT_UInt32):
	static_cast<unsigned int*>(m_data[band])[index]=static_cast<unsigned int>(dvalue);
	break;
      case(GDT_Float32):
	static_cast<float*>(m_data[band])[index]=static_cast<float>(dvalue);
	break;
      case(GDT_Float64):
	static_cast<double*>(m_data[band])[index]=static_cast<double>(dvalue);
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

/**
 * @param[in] buffer The vector with all cell values to write
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
template<typename T> bool ImgWriterGdal::writeData(std::vector<T>& buffer, int row, int band)
{
  return writeData(buffer,0,nrOfCol()-1,row,band);
}

/**
 * @param[in] buffer2d Two dimensional vector of type Vector2d (stl vector of stl vector) representing [row][col]. This vector contains all cell values that must be written
 * @param[in] minCol First column from where to start writing (counting starts from 0)
 * @param[in] maxCol Last column that must be written (counting starts from 0)
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
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
	else if(m_filename.size())
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
	  static_cast<unsigned char*>(m_data[band])[index]=static_cast<unsigned char>(dvalue);
	  break;
	case(GDT_Int16):
	  static_cast<short*>(m_data[band])[index]=static_cast<short>(dvalue);
	  break;
	case(GDT_UInt16):
	  static_cast<unsigned short*>(m_data[band])[index]=static_cast<unsigned short>(dvalue);
	  break;
	case(GDT_Int32):
	  static_cast<int*>(m_data[band])[index]=static_cast<int>(dvalue);
	  break;
	case(GDT_UInt32):
	  static_cast<unsigned int*>(m_data[band])[index]=static_cast<unsigned int>(dvalue);
	  break;
	case(GDT_Float32):
	  static_cast<float*>(m_data[band])[index]=static_cast<float>(dvalue);
	  break;
	case(GDT_Float64):
	  static_cast<double*>(m_data[band])[index]=static_cast<double>(dvalue);
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

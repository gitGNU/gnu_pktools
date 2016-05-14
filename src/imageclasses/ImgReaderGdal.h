/**********************************************************************
ImgReaderGdal.h: class to read raster files using GDAL API library
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
#ifndef _IMGREADERGDAL_H_
#define _IMGREADERGDAL_H_

#include "ImgRasterGdal.h"
#include <assert.h>
#include <fstream>
#include <string>
#include <sstream>
#include "gdal_priv.h"
#include "base/Vector2d.h"

/**
   Class to read a raster dataset in a format supported by GDAL. Data are cached in memory for a number of rows (if memory>0) before read from file.

   This class inherits from ImgRasterGdal, a general raster class to store e.g., filename, number of columns, rows and bands of the dataset. 
**/
class ImgReaderGdal : public virtual ImgRasterGdal
{
public:
  ///default constructor. Image needs to be opened later with one of the open methods.
  ImgReaderGdal(void);
  ///constructor opening an image. Set memory (in MB) to cache a number of rows in memory
  ImgReaderGdal(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly){open(filename, readMode);};
  ///destructor
  ~ImgReaderGdal(void);
  ///Open an image.
  void open(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly);
  ///Set the memory (in MB) to cache a number of rows in memory

  ///Close the image.
  void close(void);
  ///Read a single pixel cell value at a specific column and row for a specific band (all indices start counting from 0)
  template<typename T> void readData(T& value, int col, int row, int band=0);
  ///Read pixel cell values for a range of columns for a specific row and band (all indices start counting from 0)
  template<typename T> void readData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  ///Read pixel cell values for a range of columns for a specific row and band (all indices start counting from 0). The row counter can be floating, in which case a resampling is applied at the row level. You still must apply the resampling at column level. This function will be deprecated, as the GDAL API now supports rasterIO resampling (see http://www.gdal.org/structGDALRasterIOExtraArg.html)
  template<typename T> void readData(std::vector<T>& buffer, int minCol, int maxCol, double row, int band=0, RESAMPLE resample=NEAR);
  ///Read pixel cell values for a range of columns and rows for a specific band (all indices start counting from 0). The buffer is a two dimensional vector (stl vector of stl vector) representing [row][col].
  template<typename T> void readDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  ///Read pixel cell values for a range of columns and rows for a specific band (all indices start counting from 0). The buffer is a one dimensional stl vector representing all pixel values read starting from upper left to lower right.
  template<typename T> void readDataBlock(std::vector<T>& buffer , int minCol, int maxCol, int minRow, int maxRow, int band=0);
  ///Read pixel cell values for an entire row for a specific band (all indices start counting from 0)
  template<typename T> void readData(std::vector<T>& buffer, int row, int band=0);
  ///Read pixel cell values for an entire row for a specific band (all indices start counting from 0). The row counter can be floating, in which case a resampling is applied at the row level. You still must apply the resampling at column level. This function will be deprecated, as the GDAL API now supports rasterIO resampling (see http://www.gdal.org/structGDALRasterIOExtraArg.html)
  template<typename T> void readData(std::vector<T>& buffer, double row, int band=0, RESAMPLE resample=NEAR);
  ///Get the minimum and maximum cell values for a specific band in a region of interest defined by startCol, endCol, startRow and endRow (all indices start counting from 0).
  void getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue);
  ///Get the minimum and maximum cell values for a specific band (all indices start counting from 0).
  void getMinMax(double& minValue, double& maxValue, int band=0);
  ///Get the minimum cell values for a specific band and report the column and row in which the minimum value was found (all indices start counting from 0).
  double getMin(int& col, int& row, int band=0);
  ///Get the maximum cell values for a specific band and report the column and row in which the maximum value was found (all indices start counting from 0).
  double getMax(int& col, int& row, int band=0);
  ///Calculate the image histogram for a specific band using a defined number of bins and constrained   by a minimum and maximum value. A kernel density function can also be applied (default is false).
  double getHistogram(std::vector<double>& histvector, double& min, double& max,unsigned int& nbin, int theBand=0, bool kde=false);
  ///Calculate the reference pixel as the centre of gravity pixel (weighted average of all values not taking into account no data values) for a specific band (start counting from 0).
  void getRefPix(double& refX, double &refY, int band=0);
  ///Calculate the range of cell values in the image for a specific band (start counting from 0).
  void getRange(std::vector<short>& range, int Band=0);
  ///Calculate the number of valid pixels (with a value not defined as no data).
  unsigned long int getNvalid(int band);

protected:
  ///Set GDAL dataset number of columns, rows, bands and geotransform.
  void setCodec(const GDALAccess& readMode=GA_ReadOnly);
  ///Block size to cache pixel cell values in memory (calculated from user provided memory size in MB)
private:
};

/**
 * @param[out] value The cell value that was read
 * @param[in] col The column number to read (counting starts from 0)
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> void ImgReaderGdal::readData(T& value, int col, int row, int band)
{
  assert(band<nrOfBand()+1);
  assert(col<nrOfCol());
  assert(col>=0);
  assert(row<nrOfRow());
  assert(row>=0);
  double dvalue=0;
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band||m_offset.size()>band){
    if(m_scale.size()>band)
      theScale=m_scale[band];
    if(m_offset.size()>band)
      theOffset=m_offset[band];
  }
  //fetch raster band
  GDALRasterBand  *poBand;
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  poBand->RasterIO(GF_Read,col,row,1,1,&value,1,1,getGDALDataType<T>(),0,0);
  dvalue=theScale*value+theOffset;
  value=static_cast<T>(dvalue);
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band)
{
  assert(band<nrOfBand()+1);
  assert(minCol<nrOfCol());
  assert(minCol>=0);
  assert(maxCol<nrOfCol());
  assert(minCol<=maxCol);
  assert(row<nrOfRow());
  assert(row>=0);
  double theScale=1;
  double theOffset=0;
  if(m_scale.size()>band||m_offset.size()>band){
    if(m_scale.size()>band)
      theScale=m_scale[band];
    if(m_offset.size()>band)
      theOffset=m_offset[band];
  }
  //fetch raster band
  GDALRasterBand  *poBand;
  poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
  if(buffer.size()!=maxCol-minCol+1)
    buffer.resize(maxCol-minCol+1);
  poBand->RasterIO(GF_Read,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,getGDALDataType<T>(),0,0);
  if(m_scale.size()>band||m_offset.size()>band){
    for(int index=0;index<buffer.size();++index)
      buffer[index]=theScale*static_cast<double>(buffer[index])+theOffset;
  }
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 * @param[in] resample The resampling method (currently only BILINEAR and NEAR are supported)
 **/
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

/**
 * @param[out] buffer2d Two dimensional vector of type Vector2d (stl vector of stl vector) representing [row][col]. This vector contains all cell values that were read
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] minRow First row from where to start reading (counting starts from 0)
 * @param[in] maxRow Last row that must be read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
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
  
/**
 * @param[out] buffer One dimensional vector representing all pixel values read starting from upper left to lower right.
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] minRow First row from where to start reading (counting starts from 0)
 * @param[in] maxRow Last row that must be read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
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

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, int row, int band)
{
  readData(buffer,0,nrOfCol()-1,row,band);
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 * @param[in] resample The resampling method (currently only BILINEAR and NEAR are supported).
 **/
template<typename T> void ImgReaderGdal::readData(std::vector<T>& buffer, double row, int band, RESAMPLE resample)
{
  readData(buffer,0,nrOfCol()-1,row,band,resample);
}

#endif // _IMGREADERGDAL_H_

/**********************************************************************
ImgRasterGdal.h: class to read raster files using GDAL API library
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
#ifndef _IMGRASTERGDAL_H_
#define _IMGRASTERGDAL_H_

#include <typeinfo>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <assert.h>
#include "gdal_priv.h"

enum RESAMPLE { NEAR = 0, BILINEAR = 1, BICUBIC = 2 };

/**
 * @param C++ data type to be converted to GDAL data type
 * @return the GDAL data type that corresponds to the given C++ data type
 **/
template<typename T1> GDALDataType getGDALDataType(){
  if (typeid(T1) == typeid(char))
    return GDT_Byte;
  else if (typeid(T1) == typeid(unsigned char))
    return GDT_Byte;
  else if (typeid(T1) == typeid(unsigned short))
    return GDT_UInt16;
  else if (typeid(T1) == typeid(short))
    return GDT_Int16;
  else if (typeid(T1) == typeid(int))
    return GDT_Int32;
  else if (typeid(T1) == typeid(unsigned int))
    return GDT_UInt32;
  else if (typeid(T1) == typeid(long))
    return GDT_Int32;
  else if (typeid(T1) == typeid(unsigned long))
    return GDT_UInt32;
  else if (typeid(T1) == typeid(float))
    return GDT_Float32;
  else if (typeid(T1) == typeid(double))
    return GDT_Float64;
  else
    return GDT_Byte;
};

/**
   Base class for raster dataset (read and write) in a format supported by GDAL. This general raster class is used to store e.g., filename, number of columns, rows and bands of the dataset. 
**/
class ImgRasterGdal
{
public:
  ///default constructor
  ImgRasterGdal(void);
  ///destructor
  virtual ~ImgRasterGdal(void);
  ///Set scale for a specific band when writing the raster data values. The scaling and offset are applied on a per band basis. You need to set the scale for each band. If the image data are cached (class was created with memory>0), the scaling is applied on the cached memory.
  void setScale(double theScale, int band=0){
    if(m_scale.size()!=nrOfBand()){//initialize
      m_scale.resize(nrOfBand());
      for(int iband=0;iband<nrOfBand();++iband)
       m_scale[iband]=1.0;
    }
    m_scale[band]=theScale;
  };
  ///Set offset for a specific band when writing the raster data values. The scaling and offset are applied on a per band basis. You need to set the offset for each band. If the image data are cached (class was created with memory>0), the offset is applied on the cached memory.
  void setOffset(double theOffset, int band=0){
    if(m_offset.size()!=nrOfBand()){
      m_offset.resize(nrOfBand());
      for(int iband=0;iband<nrOfBand();++iband)
       m_offset[iband]=0.0;
    }
      m_offset[band]=theOffset;
  };

  ///Close the image.
  virtual void close(void);
  ///Get the filename of this dataset
  std::string getFileName() const {return m_filename;};
  ///Get the number of columns of this dataset
  int nrOfCol(void) const { return m_ncol;};
  ///Get the number of rows of this dataset
  int nrOfRow(void) const { return m_nrow;};
  ///Get the number of bands of this dataset
  int nrOfBand(void) const { return m_nband;};
  ///Is this dataset georeferenced (pixel size in y must be negative) ?
  bool isGeoRef() const {double gt[6];getGeoTransform(gt);if(gt[5]<0) return true;else return false;};
  ///Get the projection string (deprecated, use getProjectionRef instead)
  std::string getProjection(void) const;
  ///Get the projection reference
  std::string getProjectionRef(void) const;
  ///Get the geotransform data for this dataset as a string
  std::string getGeoTransform() const;
  ///Get the geotransform data for this dataset
  void getGeoTransform(double* gt) const;
  ///Set the geotransform data for this dataset
  CPLErr setGeoTransform(double* gt);
  ///Copy geotransform information from another georeferenced image
  void copyGeoTransform(const ImgRasterGdal& imgSrc);
  ///Set the projection for this dataset in well known text (wkt) format
  CPLErr setProjection(const std::string& projection);
  ///Set the projection for this dataset from user input (supports epsg:<number> format)
  CPLErr setProjectionProj4(const std::string& projection);
  ///Get the bounding box of this dataset in georeferenced coordinates
  bool getBoundingBox (double& ulx, double& uly, double& lrx, double& lry) const;
  ///Get the center position of this dataset in georeferenced coordinates
  bool getCenterPos(double& x, double& y) const;
  ///Get the upper left corner x (georeferenced) coordinate of this dataset
  double getUlx() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(ulx);};
  ///Get the upper left corner y (georeferenced) coordinate of this dataset
  double getUly() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(uly);};
  ///Get the lower right corner x (georeferenced) coordinate of this dataset
  double getLrx() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(lrx);};
  ///Get the lower right corner y (georeferenced) coordinate of this dataset
  double getLry() const {double ulx, uly, lrx,lry;getBoundingBox(ulx,uly,lrx,lry);return(lry);};
  ///Get the no data values of this dataset as a standard template library (stl) vector
  int getNoDataValues(std::vector<double>& noDataValues) const;
  ///Check if value is nodata in this dataset
  bool isNoData(double value) const{if(m_noDataValues.empty()) return false;else return find(m_noDataValues.begin(),m_noDataValues.end(),value)!=m_noDataValues.end();};
  ///Push a no data value for this dataset
  int pushNoDataValue(double noDataValue);
  ///Set the no data values of this dataset using a standard template library (stl) vector as input
  int setNoData(const std::vector<double> nodata){m_noDataValues=nodata; return(m_noDataValues.size());};
  ///Set the GDAL (internal) no data value for this data set. Only a single no data value per band is supported.
  CPLErr GDALSetNoDataValue(double noDataValue, int band=0) {return getRasterBand(band)->SetNoDataValue(noDataValue);};
  ///Check if a geolocation is covered by this dataset. Only the bounding box is checked, irrespective of no data values.
  bool covers(double x, double y) const;
  ///Check if a region of interest is (partially) covered by this dataset. Only the bounding box is checked, irrespective of no data values.
  bool covers(double ulx, double  uly, double lrx, double lry) const;
  ///Convert georeferenced coordinates (x and y) to image coordinates (column and row)
  bool geo2image(double x, double y, double& i, double& j) const;
  ///Convert image coordinates (column and row) to georeferenced coordinates (x and y)
  bool image2geo(double i, double j, double& x, double& y) const;
  ///Get the pixel cell spacing in x
  double getDeltaX(void) const {double gt[6];getGeoTransform(gt);return gt[1];};
  ///Get the pixel cell spacing in y
  double getDeltaY(void) const {double gt[6];getGeoTransform(gt);return -gt[5];};
  ///Get the GDAL datatype for this dataset
  GDALDataType getDataType(int band=0) const;
  ///Get the GDAL rasterband for this dataset
  GDALRasterBand* getRasterBand(int band=0) const;
  ///Get the GDAL color table for this dataset as an instance of the GDALColorTable class
  GDALColorTable* getColorTable(int band=0) const;
  ///Get the GDAL driver description of this dataset
  std::string getDriverDescription() const;
  ///Get the image type (implemented as the driver description)
  std::string getImageType() const{return getDriverDescription();};
  ///Get the band coding (interleave)
  std::string getInterleave() const;
  ///Get the compression from the metadata of this dataset
  std::string getCompression() const;
  //Get a pointer to the GDAL dataset
  GDALDataset* getDataset(){return m_gds;};
  ///Get the metadata of this dataset
  char** getMetadata();
  ///Get the metadata of this dataset (const version)
  char** getMetadata() const;
  ///Get the metadata of this dataset in the form of a list of strings (const version)
  void getMetadata(std::list<std::string>& metadata) const;
  ///Get the image description from the driver of this dataset
  std::string getDescription() const;
  ///Get metadata item of this dataset
  std::string getMetadataItem() const;
  ///Get the image description from the metadata of this dataset
  std::string getImageDescription() const;
  unsigned int getBlockSize() const{return m_blockSize;};
  int getBlockSizeX(int band=0)
  {
    int blockSizeX, blockSizeY;
    getRasterBand(band)->GetBlockSize( &blockSizeX, &blockSizeY );
    return blockSizeX;
  }
  int getBlockSizeY(int band=0)
  {
    int blockSizeX, blockSizeY;
    getRasterBand(band)->GetBlockSize( &blockSizeX, &blockSizeY );
    return blockSizeY;
  }
  int nrOfBlockX(int band=0)
  {
    int nXBlockSize, nYBlockSize;
    getRasterBand(band)->GetBlockSize( &nXBlockSize, &nYBlockSize );
    int nXBlocks = (nrOfCol() + nXBlockSize - 1) / nXBlockSize;
    return nXBlocks;
  }
  int nrOfBlockY(int band=0)
  {
    int nXBlockSize, nYBlockSize;
    getRasterBand(band)->GetBlockSize( &nXBlockSize, &nYBlockSize );
    int nYBlocks = (nrOfRow() + nYBlockSize - 1) / nYBlockSize;
    return nYBlocks;
  }

  friend class ImgReaderGdal;
  friend class ImgWriterGdal;

protected:
  ///filename of this dataset
  std::string m_filename;
  ///instance of the GDAL dataset of this dataset
  GDALDataset *m_gds;
  ///number of columns in this dataset
  int m_ncol;
  ///number of rows in this dataset
  int m_nrow;
  ///number of bands in this dataset
  int m_nband;
  ///GDAL data type for this dataset
  GDALDataType m_dataType;
  ///geotransform information of this dataset
  double m_gt[6];
  //projection string in wkt format
  std::string m_projection;
  ///no data values for this dataset
  std::vector<double> m_noDataValues;
  ///Vector containing the scale factor to be applied (one scale value for each band)
  std::vector<double> m_scale;
  ///Vector containing the offset factor to be applied (one offset value for each band)
  std::vector<double> m_offset;

  ///Block size to cache pixel cell values in memory (calculated from user provided memory size in MB)
  unsigned int m_blockSize;
  ///The cached pixel cell values for a certain block: a vector of void pointers (one void pointer for each band)
  std::vector<void *> m_data;
  ///first line that has been read in cache for a specific band
  std::vector<unsigned int> m_begin;
  ///beyond last line read in cache for a specific band
  std::vector<unsigned int> m_end;


private:
};

#endif // _IMGRASTERGDAL_H_

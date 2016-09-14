/**********************************************************************
ImgRaster.h: class to read/write raster files using GDAL API library
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
#ifndef _IMGRASTER_H_
#define _IMGRASTER_H_

#include <typeinfo>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <list>
#include <algorithm>
#include <utility>
#include <memory>
#include <assert.h>
#include "gdal_priv.h"
#include "base/Vector2d.h"
#include "ImgReaderOgr.h"
#include "apps/AppFactory.h"

namespace app{
class AppFactory;
}

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
class ImgRaster : public std::enable_shared_from_this<ImgRaster>
{
public:
  ///default constructor
  ImgRaster();
  ///reset all member variables
  void reset(void);
  ///constructor opening an image in memory using an external data pointer (not tested yet)
  ImgRaster(void* dataPointer, int ncol, int nrow, const GDALDataType& dataType);
  //from Reader
  ImgRaster(const std::string& filename, unsigned int memory=0) : ImgRaster() {
    open(filename, memory);
  };
  // ImgRaster(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly, unsigned int memory=0) : m_writeMode(false) {open(filename, readMode, memory);};
  //from Writer
  ///constructor opening an image for writing, copying image attributes from a source image. Caching is supported when memory>0
  ImgRaster(const std::string& filename, const ImgRaster& imgSrc, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>()) : ImgRaster() {open(filename, imgSrc, memory, options);
};
  ///copy constructor opening an image for writing in memory, copying image attributes from a source image.
  ImgRaster(ImgRaster& imgSrc, bool copyData=true) : ImgRaster() {
    open(imgSrc,copyData);
  };
  ///copy constructor opening an image for writing in memory, copying image attributes from a source image.
  ImgRaster(std::shared_ptr<ImgRaster> imgSrc, bool copyData=true) : ImgRaster() {
    open(imgSrc,copyData);
  };
  ///constructor opening an image for writing, defining all image attributes. Caching is supported when memory>0
  ImgRaster(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>()) : ImgRaster() {open(filename, ncol, nrow, nband, dataType, imageType, memory, options);
  };
  ///constructor opening an image for writing in memory, defining all image attributes
  ImgRaster(int ncol, int nrow, int nband, const GDALDataType& dataType) : ImgRaster() {open(ncol, nrow, nband, dataType);
  };

  ///destructor
  virtual ~ImgRaster(void){freeMem();};

  ///Initialize the memory for read/write image in cache
  void initMem(unsigned int memory);
  ///assignment operator
  ImgRaster& operator=(ImgRaster& imgSrc);
  ///get write mode
  bool writeMode(){return m_writeMode;};
  ///check if data pointer has been initialized
  bool isInit(){return(m_data.size()>0);};
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

  ///Open image from allocated memory instead of from file. This will allow in place image processing in memory (streaming). Notice that an extra call must be made to set the geotranform and projection. This function has not been tested yet!
  //void open(void* dataPointer, int ncol, int nrow, int nband, const GDALDataType& dataType);
  ///Close the image.
  void close(void);
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
  void copyGeoTransform(const ImgRaster& imgSrc);
  ///Copy geotransform information from another georeferenced image pointer
  void copyGeoTransform(const std::shared_ptr<ImgRaster>& imgSrc);
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
  ///Get the scale as a standard template library (stl) vector
  void getScale(std::vector<double>& scale) const {scale=m_scale;};
  ///Get the offset as a standard template library (stl) vector
  void getOffset(std::vector<double>& offset) const {offset=m_offset;};
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
  ///Get the datapointer
  void* getDataPointer(int band=0){return(m_data[band]);};
  ///free memory os data pointer
  void freeMem();
  ///Copy data 
  void copyData(void* data, int band=0);
  ///Copy data 
  // void copyData(ImgRaster& imgRaster, int band=0);
  //todo: introduce smart pointer instead of void*
  // std::unique_ptr<void> getDataPointer(int band=0){return(m_data[band]);};
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
  char** getMetadata() const;
  // Get the metadata of this dataset in the form of a list of strings (const version)
  void getMetadata(std::list<std::string>& metadata) const;
  ///Get the image description from the driver of this dataset
  std::string getDescription() const;
  ///Get metadata item of this dataset
  std::string getMetadataItem() const;
  ///Get the image description from the metadata of this dataset
  std::string getImageDescription() const;
  int getBlockSize() const{return m_blockSize;};
  int getBlockSizeX(int band=0)
  {
    int blockSizeX=0;
    int blockSizeY=0;
    if(getRasterBand(band))
      getRasterBand(band)->GetBlockSize( &blockSizeX, &blockSizeY );
    return blockSizeX;
  }
  int getBlockSizeY(int band=0)
  {
    int blockSizeX=0;
    int blockSizeY=0;
    if(getRasterBand(band))
      getRasterBand(band)->GetBlockSize( &blockSizeX, &blockSizeY );
    return blockSizeY;
  }
  int nrOfBlockX(int band=0)
  {
    int blockSizeX=0;
    int blockSizeY=0;
    int nXBlocks=0;
    if(getRasterBand(band)){
      getRasterBand(band)->GetBlockSize( &blockSizeX, &blockSizeY );
      nXBlocks = (nrOfCol() + blockSizeX - 1) / blockSizeX;
    }
    return nXBlocks;
  }
  int nrOfBlockY(int band=0)
  {
    int blockSizeX=0;
    int blockSizeY=0;
    int nYBlocks=0;
    if(getRasterBand(band)){
      getRasterBand(band)->GetBlockSize( &blockSizeX, &blockSizeY );
      nYBlocks = (nrOfRow() + blockSizeY - 1) / blockSizeY;
    }
    return nYBlocks;
  }
  ///Clone as new shared pointer to ImgRaster object
  /** 
   * 
   * @return shared pointer to new ImgRaster object
   */  
  virtual std::shared_ptr<ImgRaster> clone() {
    return(cloneImpl());
  };
  ///Create new shared pointer to ImgRaster object
  /** 
   * 
   * @return shared pointer to new ImgRaster object
   */  
  static std::shared_ptr<ImgRaster> createImg() {
    return(std::make_shared<ImgRaster>());
  };
  //From Reader
  ///Open an image. Set memory (in MB) to cache a number of rows in memory
  CPLErr open(const std::string& filename, unsigned int memory=0);
  // void open(const std::string& filename, const GDALAccess& readMode=GA_ReadOnly, unsigned int memory=0);
  ///Read all pixels from image in memory for specific band
  CPLErr readData(int band);
  ///Read all pixels from image in memory
  CPLErr readData();
  ///Read a single pixel cell value at a specific column and row for a specific band (all indices start counting from 0)
  template<typename T> CPLErr readData(T& value, int col, int row, int band=0);
  ///Return a single pixel cell value at a specific column and row for a specific band (all indices start counting from 0)
  double readData(int col, int row, int band=0){
    double value;
    readData(value, col, row, band);
    return(value);
  };
  ///Read pixel cell values for a range of columns for a specific row and band (all indices start counting from 0)
  template<typename T> CPLErr readData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  ///Read pixel cell values for a range of columns for a specific row and band (all indices start counting from 0). The row counter can be floating, in which case a resampling is applied at the row level. You still must apply the resampling at column level. This function will be deprecated, as the GDAL API now supports rasterIO resampling (see http://www.gdal.org/structGDALRasterIOExtraArg.html)
  template<typename T> CPLErr readData(std::vector<T>& buffer, int minCol, int maxCol, double row, int band, RESAMPLE resample);
  ///Read pixel cell values for a range of columns and rows for a specific band (all indices start counting from 0). The buffer is a two dimensional vector (stl vector of stl vector) representing [row][col].
  template<typename T> CPLErr readDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  ///Read pixel cell values for a range of columns and rows for a specific band (all indices start counting from 0). The buffer is a one dimensional stl vector representing all pixel values read starting from upper left to lower right.
  template<typename T> CPLErr readDataBlock(std::vector<T>& buffer , int minCol, int maxCol, int minRow, int maxRow, int band=0);
  ///Read pixel cell values for an entire row for a specific band (all indices start counting from 0)
  template<typename T> CPLErr readData(std::vector<T>& buffer, int row, int band=0);
  ///Read pixel cell values for an entire row for a specific band (all indices start counting from 0). The row counter can be floating, in which case a resampling is applied at the row level. You still must apply the resampling at column level. This function will be deprecated, as the GDAL API now supports rasterIO resampling (see http://www.gdal.org/structGDALRasterIOExtraArg.html)
  template<typename T> CPLErr readData(std::vector<T>& buffer, double row, int band, RESAMPLE resample);
  ///Get the minimum and maximum cell values for a specific band in a region of interest defined by startCol, endCol, startRow and endRow (all indices start counting from 0).
  void getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue);
  ///Get the minimum and maximum cell values for a specific band (all indices start counting from 0).
  void getMinMax(double& minValue, double& maxValue, int band=0);
  ///Get the minimum cell values for a specific band and report the column and row in which the minimum value was found (all indices start counting from 0).
  double getMin(int& col, int& row, int band=0);
  ///Get the maximum cell values for a specific band and report the column and row in which the maximum value was found (all indices start counting from 0).
  double getMax(int& col, int& row, int band=0);
  ///Calculate the image histogram for a specific band using a defined number of bins and constrained   by a minimum and maximum value. A kernel density function can also be applied (default is false).
  double getHistogram(std::vector<double>& histvector, double& min, double& max,int& nbin, int theBand=0, bool kde=false);
  ///Calculate the reference pixel as the centre of gravity pixel (weighted average of all values not taking into account no data values) for a specific band (start counting from 0).
  void getRefPix(double& refX, double &refY, int band=0);
  ///Calculate the range of cell values in the image for a specific band (start counting from 0).
  void getRange(std::vector<short>& range, int band=0);
  ///Calculate the number of valid pixels (with a value not defined as no data).
  unsigned long int getNvalid(int band);
  ///Calculate the number of invalid pixels (with a value defined as no data).
  unsigned long int getNinvalid(int band);

  //From Writer
  ///Open an image for writing, copying image attributes from a source image. Image is directly written to file. Use the constructor with memory>0 to support caching
  CPLErr open(const std::string& filename, const ImgRaster& imgSrc, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing, copying image attributes from a source image. Caching is supported when memory>0
  CPLErr open(const std::string& filename, const ImgRaster& imgSrc, unsigned int memory, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing, defining all image attributes. Image is directly written to file. Use the constructor with memory>0 to support caching
  // void open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing, defining all image attributes. Caching is supported when memory>0
  CPLErr open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>());
  ///Open an image for writing in memory, defining image attributes.
  CPLErr open(int ncol, int nrow, int nband, const GDALDataType& dataType);
  ///Open an image for writing in memory, copying image attributes from a source image.
  CPLErr open(ImgRaster& imgSrc,  bool copyData=true);
  ///Open an image for writing in memory, copying image attributes from a source image.
  CPLErr open(std::shared_ptr<ImgRaster> imgSrc,  bool copyData=true);
  ///Open an image for writing using an external data pointer (not tested yet)
  /* void open(void* dataPointer, int ncol, int nrow, const GDALDataType& dataType); */
  CPLErr open(void* dataPointer, int ncol, int nrow, const GDALDataType& dataType);
  ///Set the image description (only for GeoTiff format: TIFFTAG_IMAGEDESCRIPTION)
  void setImageDescription(const std::string& imageDescription){m_gds->SetMetadataItem( "TIFFTAG_IMAGEDESCRIPTION",imageDescription.c_str());};

  ///Write a single pixel cell value at a specific column and row for a specific band (all indices start counting from 0)
  template<typename T> CPLErr writeData(const T& value, int col, int row, int band=0);
  ///Write pixel cell values for a range of columns for a specific row and band (all indices start counting from 0)
  template<typename T> CPLErr writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band=0);
  ///Write pixel cell values for an entire row for a specific band (all indices start counting from 0)
  template<typename T> CPLErr writeData(std::vector<T>& buffer, int row, int band=0);
  // deprecated? Write an entire image from memory to file
  // CPLErr writeData(void* pdata, const GDALDataType& dataType, int band=0);
  ///Write pixel cell values for a range of columns and rows for a specific band (all indices start counting from 0). The buffer is a two dimensional vector (stl vector of stl vector) representing [row][col].
  template<typename T> CPLErr writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band=0);
  ///Prepare image writer to write to file
  CPLErr setFile(const std::string& filename, const std::string& imageType, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>());
  ///Prepare image writer to write to file
  // void setFile(const std::string& filename, const ImgRaster& imgSrc, unsigned int memory=0, const std::vector<std::string>& options=std::vector<std::string>());
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
  ///Apply thresholds: set to no data if not within thresholds t1 and t2
  CPLErr setThreshold(double t1, double t2, double value);

  //lib functions
  ///extract pixel values from raster image from a vector sample
  CPLErr extractOgr(const app::AppFactory& app);
  ///extract pixel values from raster image from a raster sample
  CPLErr extractImg(const app::AppFactory& app);
  ///calculate statistics profile based on multiband raster dataset
  CPLErr statProfile(std::shared_ptr<ImgRaster> imgWriter, const app::AppFactory& app);
  ///calculate statistics profile based on multiband raster dataset only for in memory
  std::shared_ptr<ImgRaster> statProfile(const app::AppFactory& app);
  ///filter raster dataset
  CPLErr filter(std::shared_ptr<ImgRaster> imgWriter, const app::AppFactory& app);
  ///filter raster dataset only for in memory
  std::shared_ptr<ImgRaster> filter(const app::AppFactory& app);
  ///check the difference between two images (validate in case of classification image)
  CPLErr diff(std::shared_ptr<ImgRaster> imgReference, const app::AppFactory& app);
  ///svm raster dataset
  CPLErr svm(std::shared_ptr<ImgRaster> imgWriter, const app::AppFactory& app);
  ///svm raster dataset only for in memory
  std::shared_ptr<ImgRaster> svm(const app::AppFactory& app);
  ///create shared pointer to ImgRaster with random values
  static CPLErr createImg(std::shared_ptr<ImgRaster> imgWriter, const app::AppFactory& app);
  ///create shared pointer to ImgRaster with random values only for in memory
  static std::shared_ptr<ImgRaster> createImg(const app::AppFactory& app);
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
  ///image type for this dataset
  std::string m_imageType;
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
  int m_blockSize;
  ///The cached pixel cell values for a certain block: a vector of void pointers (one void pointer for each band)
  std::vector<void*> m_data;
  //todo: use smart pointer
  //std::vector<std::unique_ptr<void[]> > m_data;
  ///first line that has been read in cache for a specific band
  std::vector<int> m_begin;
  ///beyond last line read in cache for a specific band
  std::vector<int> m_end;

  //From Reader
  ///register driver for GDAl
  void registerDriver();
  ///Create options
  std::vector<std::string> m_options;
  ///We are writing a physical file
  bool m_writeMode;

  ///Read new block in cache (defined by m_begin and m_end)
  CPLErr readNewBlock(int row, int band);
  ///Write new block from cache (defined by m_begin and m_end)
  CPLErr writeNewBlock(int row, int band);

 private:
  virtual std::shared_ptr<ImgRaster> cloneImpl() {
    /* return(std::shared_ptr<ImgRaster>(new ImgRaster(*this,false) )); */
    return(std::make_shared<ImgRaster>(*this,false));
  };
};

/**
 * @param[out] value The cell value that was read
 * @param[in] col The column number to read (counting starts from 0)
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> CPLErr ImgRaster::readData(T& value, int col, int row, int band)
{
  CPLErr returnValue=CE_None;
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
  if(m_data.size()){
    //only support random access reading if entire image is in memory for performance reasons
    if(m_blockSize!=nrOfRow()){
      std::ostringstream s;
      s << "Error: increase memory to support random access reading (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
      throw(s.str());
    }
    if(row<m_begin[band]||row>=m_end[band]){
      if(m_filename.size())
        readNewBlock(row,band);
    }
    int index=(row-m_begin[band])*nrOfCol()+col;
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
    returnValue=poBand->RasterIO(GF_Read,col,row,1,1,&value,1,1,getGDALDataType<T>(),0,0);
    dvalue=theScale*value+theOffset;
    value=static_cast<T>(dvalue);
  }
  return(returnValue);
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> CPLErr ImgRaster::readData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band)
{
  CPLErr returnValue=CE_None;
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
  if(m_data.size()){
    if(row<m_begin[band]||row>=m_end[band]){
      if(m_filename.size())
        returnValue=readNewBlock(row,band);
    }
    if(buffer.size()!=maxCol-minCol+1)
      buffer.resize(maxCol-minCol+1);
    int index=(row-m_begin[band])*nrOfCol();
    int minindex=(index+minCol);
    int maxindex=(index+maxCol);
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
    returnValue=poBand->RasterIO(GF_Read,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,getGDALDataType<T>(),0,0);
    if(m_scale.size()>band||m_offset.size()>band){
      for(int index=0;index<buffer.size();++index)
	buffer[index]=theScale*static_cast<double>(buffer[index])+theOffset;
    }
  }
  return(returnValue);
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 * @param[in] resample The resampling method (currently only BILINEAR and NEAR are supported)
 **/
template<typename T> CPLErr ImgRaster::readData(std::vector<T>& buffer, int minCol, int maxCol, double row, int band, RESAMPLE resample)
{
  CPLErr returnValue=CE_None;
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
    returnValue=readData(readBuffer_upper,minCol,maxCol,static_cast<int>(upperRow),band);
    returnValue=readData(readBuffer_lower,minCol,maxCol,static_cast<int>(lowerRow),band);
    //do interpolation in y
    for(int icol=0;icol<maxCol-minCol+1;++icol){
      buffer[icol]=(lowerRow-row+0.5)*readBuffer_upper[icol]+(1-lowerRow+row-0.5)*readBuffer_lower[icol];
    }
    break;
  default:
    returnValue=readData(buffer,minCol,maxCol,static_cast<int>(row),band);
    break;
  }
  return(returnValue);
}

/**
 * @param[out] buffer2d Two dimensional vector of type Vector2d (stl vector of stl vector) representing [row][col]. This vector contains all cell values that were read
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] minRow First row from where to start reading (counting starts from 0)
 * @param[in] maxRow Last row that must be read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> CPLErr ImgRaster::readDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  CPLErr returnValue=CE_None;
  buffer2d.resize(maxRow-minRow+1);
  typename std::vector<T> buffer;
  returnValue=readDataBlock(buffer,minCol,maxCol,minRow,maxRow,band);
  typename std::vector<T>::const_iterator startit=buffer.begin();
  typename std::vector<T>::const_iterator endit=startit;
  for(int irow=minRow;irow<=maxRow;++irow){
    buffer2d[irow-minRow].resize(maxCol-minCol+1);
    endit+=maxCol-minCol+1;
    buffer2d[irow-minRow].assign(startit,endit);
    startit+=maxCol-minCol+1;
  }
  return(returnValue);
}
  
/**
 * @param[out] buffer One dimensional vector representing all pixel values read starting from upper left to lower right.
 * @param[in] minCol First column from where to start reading (counting starts from 0)
 * @param[in] maxCol Last column that must be read (counting starts from 0)
 * @param[in] minRow First row from where to start reading (counting starts from 0)
 * @param[in] maxRow Last row that must be read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> CPLErr ImgRaster::readDataBlock(std::vector<T>& buffer, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  CPLErr returnValue=CE_None;
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
      if(irow<m_begin[band]||irow>=m_end[band]){
        if(m_filename.size())
          returnValue=readNewBlock(irow,band);
      }
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
    returnValue=poBand->RasterIO(GF_Read,minCol,minRow,maxCol-minCol+1,maxRow-minRow+1,&(buffer[0]),(maxCol-minCol+1),(maxRow-minRow+1),getGDALDataType<T>(),0,0);
    if(m_scale.size()>band||m_offset.size()>band){
      for(int index=0;index<buffer.size();++index)
	buffer[index]=theScale*buffer[index]+theOffset;
    }
  }
  return(returnValue);
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 **/
template<typename T> CPLErr ImgRaster::readData(std::vector<T>& buffer, int row, int band)
{
  return(readData(buffer,0,nrOfCol()-1,row,band));
}

/**
 * @param[out] buffer The vector with all cell values that were read
 * @param[in] row The row number to read (counting starts from 0)
 * @param[in] band The band number to read (counting starts from 0)
 * @param[in] resample The resampling method (currently only BILINEAR and NEAR are supported).
 **/
template<typename T> CPLErr ImgRaster::readData(std::vector<T>& buffer, double row, int band, RESAMPLE resample)
{
  return(readData(buffer,0,nrOfCol()-1,row,band,resample));
}

//From Writer
/**
 * @param[in] value The cell value to write
 * @param[in] col The column number to write (counting starts from 0)
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
template<typename T> CPLErr ImgRaster::writeData(const T& value, int col, int row, int band)
{
  CPLErr returnValue=CE_None;
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
    returnValue=poBand->RasterIO(GF_Write,col,row,1,1,&dvalue,1,1,getGDALDataType<T>(),0,0);
  }
  return(returnValue);
}

/**
 * @param[in] buffer The vector with all cell values to write
 * @param[in] minCol First column from where to start writing (counting starts from 0)
 * @param[in] maxCol Last column that must be written (counting starts from 0)
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
template<typename T> CPLErr ImgRaster::writeData(std::vector<T>& buffer, int minCol, int maxCol, int row, int band)
{
  CPLErr returnValue=CE_None;
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
	returnValue=writeNewBlock(row,band);
    }
    int index=(row-m_begin[band])*nrOfCol();
    int minindex=(index+minCol);
    int maxindex=(index+maxCol);
    typename std::vector<T>::const_iterator bufit=buffer.begin();
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
    returnValue=poBand->RasterIO(GF_Write,minCol,row,buffer.size(),1,&(buffer[0]),buffer.size(),1,getGDALDataType<T>(),0,0);
  }
  return(returnValue);
}

/**
 * @param[in] buffer The vector with all cell values to write
 * @param[in] row The row number to write (counting starts from 0)
 * @param[in] band The band number to write (counting starts from 0)
 * @return true if write successful
 **/
template<typename T> CPLErr ImgRaster::writeData(std::vector<T>& buffer, int row, int band)
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
template<typename T> CPLErr ImgRaster::writeDataBlock(Vector2d<T>& buffer2d, int minCol, int maxCol, int minRow, int maxRow, int band)
{
  CPLErr returnValue=CE_None;
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
          returnValue=writeNewBlock(irow,band);
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
    returnValue=poBand->RasterIO(GF_Write,minCol,minRow,maxCol-minCol+1,maxRow-minRow+1,&(buffer[0]),(maxCol-minCol+1),(maxRow-minRow+1),getGDALDataType<T>(),0,0);
  }
  return(returnValue);
}

#endif // _IMGRASTER_H_

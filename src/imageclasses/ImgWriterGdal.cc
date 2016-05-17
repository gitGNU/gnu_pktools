/**********************************************************************
ImgWriterGdal.cc: class to write raster files using GDAL API library
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
#include <iostream>
#include <iomanip>
#include <time.h>
#include <algorithm>
extern "C" {
#include "gdal_alg.h"
}
#include "ImgWriterGdal.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

ImgWriterGdal::ImgWriterGdal(void){};

ImgWriterGdal::~ImgWriterGdal(void){};

/**
 * @param row Write a new block for caching this row (if needed)
 * @param band Band that must be written in cache
 * @return true if write was successful
 **/
bool ImgWriterGdal::writeNewBlock(int row, int band)
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
  return true;//new block was written
}

// /**
//  * @param filename Open a raster dataset with this filename
//  * @param imgSrc Use this source image as a template to copy image attributes
//  * @param options Creation options
//  **/
// void ImgWriterGdal::open(const std::string& filename, const ImgReaderGdal& imgSrc, const std::vector<std::string>& options)
// {
//   m_ncol=imgSrc.nrOfCol();
//   m_nrow=imgSrc.nrOfRow();
//   m_nband=imgSrc.nrOfBand();
//   m_dataType=imgSrc.getDataType();
//   setFile(filename,imgSrc,options);
//   // m_filename=filename;
//   // m_options=options;
//   // setCodec(imgSrc);
// }

/**
 * @param filename Open a raster dataset with this filename
 * @param imgSrc Use this source image as a template to copy image attributes
 * @param memory Available memory to cache image raster data (in MB)
 * @param options Creation options
 **/
void ImgWriterGdal::open(const std::string& filename, const ImgReaderGdal& imgSrc, unsigned int memory, const std::vector<std::string>& options)
{
  m_ncol=imgSrc.nrOfCol();
  m_nrow=imgSrc.nrOfRow();
  m_nband=imgSrc.nrOfBand();
  m_dataType=imgSrc.getDataType();
  setFile(filename,imgSrc,memory,options);
}

/**
 * @param imgSrc Use this source image as a template to copy image attributes
 * @param options Creation options
 **/
void ImgWriterGdal::open(const ImgReaderGdal& imgSrc)
{
  m_ncol=imgSrc.nrOfCol();
  m_nrow=imgSrc.nrOfRow();
  m_nband=imgSrc.nrOfBand();
  m_dataType=imgSrc.getDataType();
  initMem(0);
  for(int iband=0;iband<m_nband;++iband){
    m_begin[iband]=0;
    m_end[iband]=m_begin[iband]+m_blockSize;
  }
}

// /**
//  * @param filename Open a raster dataset with this filename
//  * @param ncol Number of columns in image
//  * @param nrow Number of rows in image
//  * @param nband Number of bands in image
//  * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
//  * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
//  * @param options Creation options
//  **/
// void ImgWriterGdal::open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, const std::vector<std::string>& options)
// {
//   m_ncol = ncol;
//   m_nrow = nrow;
//   m_nband = nband;
//   m_dataType = dataType;
//   setFile(filename,imageType,options);
//   // m_filename = filename;
//   // m_options=options;
//   // setCodec(imageType);
// }

/**
 * @param filename Open a raster dataset with this filename
 * @param ncol Number of columns in image
 * @param nrow Number of rows in image
 * @param nband Number of bands in image
 * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 * @param memory Available memory to cache image raster data (in MB)
 * @param options Creation options
 **/
void ImgWriterGdal::open(const std::string& filename, int ncol, int nrow, int nband, const GDALDataType& dataType, const std::string& imageType, unsigned int memory, const std::vector<std::string>& options)
{
  m_ncol = ncol;
  m_nrow = nrow;
  m_nband = nband;
  m_dataType = dataType;
  setFile(filename,imageType,memory,options);
}

/**
 * @param ncol Number of columns in image
 * @param nrow Number of rows in image
 * @param nband Number of bands in image
 * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
 **/
void ImgWriterGdal::open(int ncol, int nrow, int nband, const GDALDataType& dataType)
{
  m_ncol = ncol;
  m_nrow = nrow;
  m_nband = nband;
  m_dataType = dataType;
  initMem(0);
  for(int iband=0;iband<m_nband;++iband){
    m_begin[iband]=0;
    m_end[iband]=m_begin[iband]+m_blockSize;
  }
}

void ImgWriterGdal::close(void)
{
  if(m_data.size()&&m_filename.size()){
    for(int iband=0;iband<nrOfBand();++iband) 
      writeNewBlock(nrOfRow(),iband);
  }
  ImgRasterGdal::close();
  char **papszOptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  if(papszOptions)
    CSLDestroy(papszOptions);
}


/**
 * @param imgSrc Use this source image as a template to copy image attributes
 **/
void ImgWriterGdal::setCodec(const ImgReaderGdal& imgSrc){
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imgSrc.getDriverDescription().c_str());
  if( poDriver == NULL ){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
  char **papszMetadata;
  papszMetadata = poDriver->GetMetadata();
  //todo: try and catch if CREATE is not supported (as in PNG)
  assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
  char **papszOptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());

  m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,imgSrc.getDataType(),papszOptions);
  double gt[6];
  imgSrc.getGeoTransform(gt);
  if(setGeoTransform(gt)!=CE_None)
    std::cerr << "Warning: could not write geotransform information in " << m_filename << std::endl;
  setProjection(imgSrc.getProjection());
  if(setProjection(imgSrc.getProjection())!=CE_None)
    std::cerr << "Warning: could not write projection information in " << m_filename << std::endl;

  if(m_noDataValues.size()){
    for(int iband=0;iband<nrOfBand();++iband)
      GDALSetNoDataValue(m_noDataValues[0],iband);
  }

  m_gds->SetMetadata(imgSrc.getMetadata() ); 
  m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
  std::string versionString="pktools ";
  versionString+=VERSION;
  versionString+=" by Pieter Kempeneers";
  m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", versionString.c_str());
  time_t rawtime;
  time ( &rawtime );

  time_t tim=time(NULL);
  tm *now=localtime(&tim);
  std::ostringstream datestream;
  //date std::string must be 20 characters long...
  datestream << now->tm_year+1900;
  if(now->tm_mon+1<10)
    datestream << ":0" << now->tm_mon+1;
  else
    datestream << ":" << now->tm_mon+1;
  if(now->tm_mday<10)
    datestream << ":0" << now->tm_mday;
  else
    datestream << ":" << now->tm_mday;
  if(now->tm_hour<10)
    datestream << " 0" << now->tm_hour;
  else
    datestream << " " << now->tm_hour;
  if(now->tm_min<10)
    datestream << ":0" << now->tm_min;
  else
    datestream << ":" << now->tm_min;
  if(now->tm_sec<10)
    datestream << ":0" << now->tm_sec;
  else
    datestream << ":" << now->tm_sec;
  m_gds->SetMetadataItem( "TIFFTAG_DATETIME", datestream.str().c_str());
  if(imgSrc.getColorTable()!=NULL)
    setColorTable(imgSrc.getColorTable());
}

/**
 * @param dataType The data type of the image (one of the GDAL supported datatypes: GDT_Byte, GDT_[U]Int[16|32], GDT_Float[32|64])
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 **/
void ImgWriterGdal::setCodec(const std::string& imageType)
{
  GDALAllRegister();
  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(imageType.c_str());
  if( poDriver == NULL ){
    std::ostringstream s;
    s << "FileOpenError (" << imageType << ")";
    throw(s.str());
  }
  char **papszMetadata;
  papszMetadata = poDriver->GetMetadata();
  //todo: try and catch if CREATE is not supported (as in PNG)
  assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
  char **papszOptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=m_options.begin();optionIt!=m_options.end();++optionIt)
    papszOptions=CSLAddString(papszOptions,optionIt->c_str());
  m_gds=poDriver->Create(m_filename.c_str(),m_ncol,m_nrow,m_nband,m_dataType,papszOptions);
  double gt[6];
  getGeoTransform(gt);
  if(setGeoTransform(gt)!=CE_None)
    std::cerr << "Warning: could not write geotransform information in " << m_filename << std::endl;
  if(setProjection(m_projection)!=CE_None)
    std::cerr << "Warning: could not write projection information in " << m_filename << std::endl;
  m_gds->SetMetadataItem( "TIFFTAG_DOCUMENTNAME", m_filename.c_str());
  std::string versionString="pktools ";
  versionString+=VERSION;
  versionString+=" by Pieter Kempeneers";
  m_gds->SetMetadataItem( "TIFFTAG_SOFTWARE", versionString.c_str());
  time_t rawtime;
  time ( &rawtime );

  time_t tim=time(NULL);
  tm *now=localtime(&tim);
  std::ostringstream datestream;
  //date std::string must be 20 characters long...
  datestream << now->tm_year+1900;
  if(now->tm_mon+1<10)
    datestream << ":0" << now->tm_mon+1;
  else
    datestream << ":" << now->tm_mon+1;
  if(now->tm_mday<10)
    datestream << ":0" << now->tm_mday;
  else
    datestream << ":" << now->tm_mday;
  if(now->tm_hour<10)
    datestream << " 0" << now->tm_hour;
  else
    datestream << " " << now->tm_hour;
  if(now->tm_min<10)
    datestream << ":0" << now->tm_min;
  else
    datestream << ":" << now->tm_min;
  if(now->tm_sec<10)
    datestream << ":0" << now->tm_sec;
  else
    datestream << ":" << now->tm_sec;
  m_gds->SetMetadataItem( "TIFFTAG_DATETIME", datestream.str().c_str());
}

/**
 * @param filename Open a raster dataset with this filename
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 **/
void ImgWriterGdal::setFile(const std::string& filename, const std::string& imageType, unsigned long int memory, const std::vector<std::string>& options)
{
  m_filename=filename;
  m_options=options;
  setCodec(imageType);
  if(m_data.empty())
    initMem(memory);
  for(int iband=0;iband<m_nband;++iband){
    m_begin[iband]=0;
    m_end[iband]=m_begin[iband]+m_blockSize;
  }
}

/**
 * @param filename Open a raster dataset with this filename
 * @param imageType Image type. Currently only those formats where the drivers support the Create method can be written
 **/
void ImgWriterGdal::setFile(const std::string& filename, const ImgReaderGdal& imgSrc, unsigned long int memory, const std::vector<std::string>& options)
{
  m_filename=filename;
  m_options=options;
  setCodec(imgSrc);
  if(m_data.empty())
    initMem(memory);
  for(int iband=0;iband<m_nband;++iband){
    m_begin[iband]=0;
    m_end[iband]=m_begin[iband]+m_blockSize;
  }
}

/**
 * @param metadata Set this metadata when writing the image (if supported byt the driver)
 **/
void ImgWriterGdal::setMetadata(char** metadata)
{
  if(m_gds)
    m_gds->SetMetadata(metadata); 
}

//default projection: ETSR-LAEA
// std::string ImgWriterGdal::setProjection(void)
// {
//   std::string theProjection;
//   OGRSpatialReference oSRS;
//   char *pszSRS_WKT = NULL;
//   //// oSRS.importFromEPSG(3035);
//   oSRS.SetGeogCS("ETRS89","European_Terrestrial_Reference_System_1989","GRS 1980",6378137,298.2572221010042,"Greenwich",0,"degree",0.0174532925199433);
//   // cout << setprecision(16) << "major axis: " << oSRS.GetSemiMajor(NULL) << endl;//notice that major axis can be set to a different value than the default to the well known standard corresponding to the name (European_Terrestrial_Reference_System_1989), but that new value, while recognized by GetSemiMajor, will not be written in the geotiff tag!
//   oSRS.SetProjCS( "ETRS89 / ETRS-LAEA" );
//   oSRS.SetLAEA(52,10,4321000,3210000);
//   oSRS.exportToWkt( &pszSRS_WKT );
//   theProjection=pszSRS_WKT;
//   CPLFree( pszSRS_WKT );
//   assert(m_gds);
//   m_gds->SetProjection(theProjection.c_str());
//   return(theProjection);
// }


/**
 * @param filename ASCII file containing 5 columns: index R G B ALFA (0:transparent, 255:solid)
 * @param band band number to set color table (starting counting from 0)
 **/
void ImgWriterGdal::setColorTable(const std::string& filename, int band)
{
  //todo: fool proof table in file (no checking currently done...)
  std::ifstream ftable(filename.c_str(),std::ios::in);
  std::string line;
//   poCT=new GDALColorTable();
  GDALColorTable colorTable;
  short nline=0;
  while(getline(ftable,line)){
    ++nline;
    std::istringstream ist(line);
    GDALColorEntry sEntry;
    short id;
    ist >> id >> sEntry.c1 >> sEntry.c2 >> sEntry.c3 >> sEntry.c4;
    colorTable.SetColorEntry(id,&sEntry);
  }
  if(m_gds)
    (m_gds->GetRasterBand(band+1))->SetColorTable(&colorTable);
}

/**
 * @param colorTable Instance of the GDAL class GDALColorTable
 * @param band band number to set color table (starting counting from 0)
 **/
void ImgWriterGdal::setColorTable(GDALColorTable* colorTable, int band)
{
  if(m_gds)
    (m_gds->GetRasterBand(band+1))->SetColorTable(colorTable);
}

// //write an entire image from memory to file
// bool ImgWriterGdal::writeData(void* pdata, const GDALDataType& dataType, int band){
//   //fetch raster band
//   GDALRasterBand  *poBand;
//   if(band>=nrOfBand()+1){
//     std::ostringstream s;
//     s << "band (" << band << ") exceeds nrOfBand (" << nrOfBand() << ")";
//     throw(s.str());
//   }
//   poBand = m_gds->GetRasterBand(band+1);//GDAL uses 1 based index
//   poBand->RasterIO(GF_Write,0,0,nrOfCol(),nrOfRow(),pdata,nrOfCol(),nrOfRow(),dataType,0,0);
//   return true;
// }  

/**
 * @param ogrReader Vector dataset as an instance of the ImgReaderOgr that must be rasterized
 * @param burnValues Values to burn into raster cells (one value for each band)
 * @param controlOptions special options controlling rasterization (ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG)
 * "ATTRIBUTE":
 * Identifies an attribute field on the features to be used for a burn in value. The value will be burned into all output bands. If specified, padfLayerBurnValues will not be used and can be a NULL pointer. 
 * "CHUNKYSIZE":
 * The height in lines of the chunk to operate on. The larger the chunk size the less times we need to make a pass through all the shapes. If it is not set or set to zero the default chunk size will be used. Default size will be estimated based on the GDAL cache buffer size using formula: cache_size_bytes/scanline_size_bytes, so the chunk will not exceed the cache.
 * "ALL_TOUCHED":
 * May be set to TRUE to set all pixels touched by the line or polygons, not just those whose center is within the polygon or that are selected by brezenhams line algorithm. Defaults to FALSE. 
 "BURN_VALUE_
 * May be set to "Z" to use the Z values of the geometries. The value from padfLayerBurnValues or the attribute field value is added to this before burning. In default case dfBurnValue is burned as it is. This is implemented properly only for points and lines for now. Polygons will be burned using the Z value from the first point. The M value may be supported in the future. 
 * "MERGE_ALG":
 * May be REPLACE (the default) or ADD. REPLACE results in overwriting of value, while ADD adds the new value to the existing raster, suitable for heatmaps for instance. 
 * @param layernames Names of the vector dataset layers to process. Leave empty to process all layers
 **/
void ImgWriterGdal::rasterizeOgr(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues, const std::vector<std::string>& controlOptions, const std::vector<std::string>& layernames ){
  std::vector<int> bands;
  if(burnValues.empty()&&controlOptions.empty()){
    std::string errorString="Error: either burn values or control options must be provided";
    throw(errorString);
  }
  for(int iband=0;iband<nrOfBand();++iband)
    bands.push_back(iband+1);
  std::vector<OGRLayerH> layers;
  int nlayer=0;

  std::vector<double> burnBands;//burn values for all bands in a single layer
  std::vector<double> burnLayers;//burn values for all bands and all layers
  if(burnValues.size()){
    burnBands=burnValues;
    while(burnBands.size()<nrOfBand())
      burnBands.push_back(burnValues[0]);
  }
  for(int ilayer=0;ilayer<ogrReader.getLayerCount();++ilayer){
    std::string currentLayername=ogrReader.getLayer(ilayer)->GetName();
    if(layernames.size())
      if(find(layernames.begin(),layernames.end(),currentLayername)==layernames.end())
	continue;
    std::cout << "processing layer " << currentLayername << std::endl;
    layers.push_back((OGRLayerH)ogrReader.getLayer(ilayer));
    ++nlayer;
    if(burnValues.size()){
      for(int iband=0;iband<nrOfBand();++iband)
        burnLayers.insert(burnLayers.end(),burnBands.begin(),burnBands.end());
    }
  }
  void* pTransformArg;
  GDALProgressFunc pfnProgress=NULL;
  void* pProgressArg=NULL;

  char **coptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=controlOptions.begin();optionIt!=controlOptions.end();++optionIt)
    coptions=CSLAddString(coptions,optionIt->c_str());

  if(controlOptions.size()){
    if(GDALRasterizeLayers( (GDALDatasetH)m_gds,nrOfBand(),&(bands[0]),layers.size(),&(layers[0]),NULL,pTransformArg,NULL,coptions,pfnProgress,pProgressArg)!=CE_None){
      std::string errorString(CPLGetLastErrorMsg());
      throw(errorString);
    }
  }
  else if(burnValues.size()){
    if(GDALRasterizeLayers( (GDALDatasetH)m_gds,nrOfBand(),&(bands[0]),layers.size(),&(layers[0]),NULL,pTransformArg,&(burnLayers[0]),NULL,pfnProgress,pProgressArg)!=CE_None){
      std::string errorString(CPLGetLastErrorMsg());
      throw(errorString);
    }
  }
  else{
    std::string errorString="Error: either attribute fieldname or burn values must be set to rasterize vector dataset";
    throw(errorString);
  }
}

/**
 * @param ogrReader Vector dataset as an instance of the ImgReaderOgr that must be rasterized
 * @param burnValues Values to burn into raster cells (one value for each band)
 * @param controlOptions special options controlling rasterization (ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG)
 * "ATTRIBUTE":
 * Identifies an attribute field on the features to be used for a burn in value. The value will be burned into all output bands. If specified, padfLayerBurnValues will not be used and can be a NULL pointer. 
 * "ALL_TOUCHED":
 * May be set to TRUE to set all pixels touched by the line or polygons, not just those whose center is within the polygon or that are selected by brezenhams line algorithm. Defaults to FALSE. 
 "BURN_VALUE_FROM":
 * May be set to "Z" to use the Z values of the geometries. The value from padfLayerBurnValues or the attribute field value is added to this before burning. In default case dfBurnValue is burned as it is. This is implemented properly only for points and lines for now. Polygons will be burned using the Z value from the first point. The M value may be supported in the future. 
 * "MERGE_ALG":
 * May be REPLACE (the default) or ADD. REPLACE results in overwriting of value, while ADD adds the new value to the existing raster, suitable for heatmaps for instance. 
 * @param layernames Names of the vector dataset layers to process. Leave empty to process all layers
 **/
void ImgWriterGdal::rasterizeBuf(ImgReaderOgr& ogrReader, const std::vector<double>& burnValues, const std::vector<std::string>& controlOptions, const std::vector<std::string>& layernames ){
  if(m_blockSize<nrOfRow()){
    std::ostringstream s;
    s << "Error: increase memory to perform rasterize in entirely in buffer (now at " << 100.0*m_blockSize/nrOfRow() << "%)";
    throw(s.str());
  }
  if(burnValues.empty()&&controlOptions.empty()){
    std::string errorString="Error: either burn values or control options must be provided";
    throw(errorString);
  }

  std::vector<OGRLayerH> layers;
  int nlayer=0;

  std::vector<double> burnBands;//burn values for all bands in a single layer
   if(burnValues.size()){
    burnBands=burnValues;
    while(burnBands.size()<nrOfBand())
      burnBands.push_back(burnValues[0]);
  }
  for(int ilayer=0;ilayer<ogrReader.getLayerCount();++ilayer){
    std::string currentLayername=ogrReader.getLayer(ilayer)->GetName();
    if(layernames.size())
      if(find(layernames.begin(),layernames.end(),currentLayername)==layernames.end())
	continue;
    std::cout << "processing layer " << currentLayername << std::endl;
    layers.push_back((OGRLayerH)ogrReader.getLayer(ilayer));
    ++nlayer;
  }
  void* pTransformArg;
  GDALProgressFunc pfnProgress=NULL;
  void* pProgressArg=NULL;

  char **coptions=NULL;
  for(std::vector<std::string>::const_iterator optionIt=controlOptions.begin();optionIt!=controlOptions.end();++optionIt)
    coptions=CSLAddString(coptions,optionIt->c_str());

  for(int iband=0;iband<nrOfBand();++iband){
    double gt[6];
    getGeoTransform(gt);
    if(controlOptions.size()){
      if(GDALRasterizeLayersBuf(m_data[iband],nrOfCol(),nrOfRow(),getDataType(),0,0,layers.size(),&(layers[0]), getProjectionRef().c_str(),gt,NULL, pTransformArg, NULL,coptions,pfnProgress,pProgressArg)!=CE_None){
        std::string errorString(CPLGetLastErrorMsg());
        throw(errorString);
      }
    }
    else if(burnValues.size()){
      if(GDALRasterizeLayersBuf(m_data[iband],nrOfCol(),nrOfRow(),getDataType(),0,0,layers.size(),&(layers[0]), getProjectionRef().c_str(),gt,NULL, pTransformArg, burnBands[iband],NULL,pfnProgress,pProgressArg)!=CE_None){
        std::string errorString(CPLGetLastErrorMsg());
        throw(errorString);
      }
    }
    else{
      std::string errorString="Error: either attribute fieldname or burn values must be set to rasterize vector dataset";
        throw(errorString);
    }
  }
}

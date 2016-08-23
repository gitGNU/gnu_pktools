/**********************************************************************
ImgCollection.cc: class to read raster files using GDAL API library
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
#include <vector>
#include <string>
#include <iostream>
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "algorithms/Egcs.h"
#include "apps/AppFactory.h"
#include "ImgCollection.h"

using namespace std;
using namespace app;

//   ///time period covering the image collection (check http://www.boost.org/doc/libs/1_55_0/doc/html/date_time/examples.html#date_time.examples.time_periods for how to use boost period)
// boost::posix_time::time_period ImgCollection::getTimePeriod(){
//   if(m_time.size()){
//     std::vector<boost::posix_time::time_period>::iterator tit=m_time.begin();
//     boost::posix_time::time_period timePeriod=*(tit++);
//     while(tit!=m_time.end()){
//       timePeriod.span(*(tit++));
//     }
//     return(timePeriod);
//   }
// }    

/**
 * @param ulx upper left coordinate in x
 * @param uly upper left coordinate in y
 * @param lrx lower left coordinate in x
 * @param lry lower left coordinate in y
 **/
void ImgCollection::getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const{
  std::vector<std::shared_ptr<ImgRaster> >::const_iterator it=begin();
  if(it!=end())
    (*(it++))->getBoundingBox(ulx,uly,lrx,lry);
  while(it!=end()){
    double imgulx,imguly,imglrx,imglry;
    (*(it++))->getBoundingBox(imgulx,imguly,imglrx,imglry);
    ulx=(ulx>imgulx)? imgulx : ulx;
    uly=(uly<imguly)? imguly : uly;
    lrx=(lrx<imglrx)? imglrx : lrx;
    lry=(lry>imglry)? imglry : lry;
  }
}    

/**
 * @param x,y georeferenced coordinates in x and y
 * @return true if image covers the georeferenced location
 **/
bool ImgCollection::covers(double x, double  y) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((x > theULX)&&
         (x < theLRX)&&
         (y < theULY)&&
         (y >theLRY));
}

/**
 * @param ulx upper left coordinate in x
 * @param uly upper left coordinate in y
 * @param lrx lower left coordinate in x
 * @param lry lower left coordinate in y
 * @return true if image (partially) covers the bounding box
 **/
bool ImgCollection::covers(double ulx, double  uly, double lrx, double lry) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((ulx < theLRX)&&(lrx > theULX)&&(lry < theULY)&&(uly > theLRY));
}

/**
 * @param noDataValues standard template library (stl) vector containing no data values
 * @return number of no data values in this dataset
 **/
unsigned int ImgCollection::getNoDataValues(std::vector<double>& noDataValues) const
{
  if(m_noDataValues.size()){
    noDataValues=m_noDataValues;
    return m_noDataValues.size();
  }
  else
    return 0;
}

/**
 * @param noDataValue no data value to be pushed for this dataset
 * @return number of no data values in this dataset
 **/
unsigned int ImgCollection::pushNoDataValue(double noDataValue)
{
  if(find(m_noDataValues.begin(),m_noDataValues.end(),noDataValue)==m_noDataValues.end())
    m_noDataValues.push_back(noDataValue);
  return(m_noDataValues.size());
}

/**
 * @param imgRaster output raster crop dataset
 * @return output image
 **/
shared_ptr<ImgRaster> ImgCollection::crop(const AppFactory& app){
  shared_ptr<ImgRaster> imgWriter(new ImgRaster);
  crop(imgWriter, app);
  return(imgWriter);
}

/**
 * @param imgRaster output raster crop dataset
 * @return 0 if successful
 **/
int ImgCollection::crop(shared_ptr<ImgRaster> imgWriter, const AppFactory& app){
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  //todo: support layer names
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file");
  Optionpk<bool> cut_opt("cut", "crop_to_cutline", "Crop the extent of the target dataset to the extent of the cutline.",false);
  Optionpk<string> eoption_opt("eo","eo", "special extent options controlling rasterization: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG, e.g., -eo ATTRIBUTE=fieldname");
  Optionpk<string> mask_opt("m", "mask", "Use the the specified file as a validity mask (0 is nodata).");
  Optionpk<double> msknodata_opt("msknodata", "msknodata", "Mask value not to consider for crop.", 0);
  Optionpk<unsigned int> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed)", 0);
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box", 0.0);
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box", 0.0);
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box", 0.0);
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box", 0.0);
  Optionpk<double>  dx_opt("dx", "dx", "Output resolution in x (in meter) (empty: keep original resolution)");
  Optionpk<double>  dy_opt("dy", "dy", "Output resolution in y (in meter) (empty: keep original resolution)");
  Optionpk<double> cx_opt("x", "x", "x-coordinate of image center to crop (in meter)");
  Optionpk<double> cy_opt("y", "y", "y-coordinate of image center to crop (in meter)");
  Optionpk<double> nx_opt("nx", "nx", "image size in x to crop (in meter)");
  Optionpk<double> ny_opt("ny", "ny", "image size in y to crop (in meter)");
  Optionpk<unsigned int> ns_opt("ns", "ns", "number of samples  to crop (in pixels)");
  Optionpk<unsigned int> nl_opt("nl", "nl", "number of lines to crop (in pixels)");
  Optionpk<unsigned int>  band_opt("b", "band", "band index to crop (leave empty to retain all bands)");
  Optionpk<unsigned int> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned int> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<double> autoscale_opt("as", "autoscale", "scale output to min and max, e.g., --autoscale 0 --autoscale 255");
  Optionpk<double> scale_opt("scale", "scale", "output=scale*input+offset");
  Optionpk<double> offset_opt("offset", "offset", "output=scale*input+offset");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<double>  nodata_opt("nodata", "nodata", "Nodata value to put in image if out of bounds.");
  Optionpk<string>  resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation).", "near");
  Optionpk<string>  description_opt("d", "description", "Set image description");
  Optionpk<bool>  align_opt("align", "align", "Align output bounding box to input image",false);
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);
  Optionpk<short>  verbose_opt("v", "verbose", "verbose", 0,2);

  extent_opt.setHide(1);
  cut_opt.setHide(1);
  eoption_opt.setHide(1);
  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  mask_opt.setHide(1);
  msknodata_opt.setHide(1);
  mskband_opt.setHide(1);
  option_opt.setHide(1);
  cx_opt.setHide(1);
  cy_opt.setHide(1);
  nx_opt.setHide(1);
  ny_opt.setHide(1);
  ns_opt.setHide(1);
  nl_opt.setHide(1);
  scale_opt.setHide(1);
  offset_opt.setHide(1);
  nodata_opt.setHide(1);
  description_opt.setHide(1);
  memory_opt.setHide(1);
  
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=projection_opt.retrieveOption(app.getArgc(),app.getArgv());
    ulx_opt.retrieveOption(app.getArgc(),app.getArgv());
    uly_opt.retrieveOption(app.getArgc(),app.getArgv());
    lrx_opt.retrieveOption(app.getArgc(),app.getArgv());
    lry_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    bstart_opt.retrieveOption(app.getArgc(),app.getArgv());
    bend_opt.retrieveOption(app.getArgc(),app.getArgv());
    autoscale_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    oformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    colorTable_opt.retrieveOption(app.getArgc(),app.getArgv());
    dx_opt.retrieveOption(app.getArgc(),app.getArgv());
    dy_opt.retrieveOption(app.getArgc(),app.getArgv());
    resample_opt.retrieveOption(app.getArgc(),app.getArgv());
    extent_opt.retrieveOption(app.getArgc(),app.getArgv());
    cut_opt.retrieveOption(app.getArgc(),app.getArgv());
    eoption_opt.retrieveOption(app.getArgc(),app.getArgv());
    mask_opt.retrieveOption(app.getArgc(),app.getArgv());
    msknodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    mskband_opt.retrieveOption(app.getArgc(),app.getArgv());
    option_opt.retrieveOption(app.getArgc(),app.getArgv());
    cx_opt.retrieveOption(app.getArgc(),app.getArgv());
    cy_opt.retrieveOption(app.getArgc(),app.getArgv());
    nx_opt.retrieveOption(app.getArgc(),app.getArgv());
    ny_opt.retrieveOption(app.getArgc(),app.getArgv());
    ns_opt.retrieveOption(app.getArgc(),app.getArgv());
    nl_opt.retrieveOption(app.getArgc(),app.getArgv());
    scale_opt.retrieveOption(app.getArgc(),app.getArgv());
    offset_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    description_opt.retrieveOption(app.getArgc(),app.getArgv());
    align_opt.retrieveOption(app.getArgc(),app.getArgv());
    memory_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    std::ostringstream helpStream;
    helpStream << "Usage: pkcrop -i input -o output" << std::endl;
    helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    throw(helpStream.str());//help was invoked, stop processing
  }
  // if(input_opt.empty()){
  if(empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(0);
  }
  // if(output_opt.empty()){
  //   std::cerr << "No output file provided (use option -o). Use --help for help information" << std::endl;
  //   exit(0);
  // }

  double nodataValue=nodata_opt.size()? nodata_opt[0] : 0;
  RESAMPLE theResample;
  if(resample_opt[0]=="near"){
    theResample=NEAR;
    if(verbose_opt[0])
      cout << "resampling: nearest neighbor" << endl;
  }
  else if(resample_opt[0]=="bilinear"){
    theResample=BILINEAR;
    if(verbose_opt[0])
      cout << "resampling: bilinear interpolation" << endl;
  }
  else{
    std::cout << "Error: resampling method " << resample_opt[0] << " not supported" << std::endl;
    exit(1);
  }

  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  // ImgReaderGdal imgReader;
  // ImgWriterGdal imgWriter;
  //open input images to extract number of bands and spatial resolution
  int ncropband=0;//total number of bands to write
  double dx=0;
  double dy=0;
  if(dx_opt.size())
    dx=dx_opt[0];
  if(dy_opt.size())
    dy=dy_opt[0];

  //convert start and end band options to vector of band indexes
  try{
    if(bstart_opt.size()){
      if(bend_opt.size()!=bstart_opt.size()){
	string errorstring="Error: options for start and end band indexes must be provided as pairs, missing end band";
	throw(errorstring);
      }
      band_opt.clear();
      for(int ipair=0;ipair<bstart_opt.size();++ipair){
	if(bend_opt[ipair]<=bstart_opt[ipair]){
	  string errorstring="Error: index for end band must be smaller then start band";
	  throw(errorstring);
	}
	for(unsigned int iband=bstart_opt[ipair];iband<=bend_opt[ipair];++iband)
	  band_opt.push_back(iband);
      }
    }
  }
  catch(string error){
    cerr << error << std::endl;
    exit(1);
  }

  bool isGeoRef=false;
  string projectionString;
  // for(int iimg=0;iimg<input_opt.size();++iimg){

  std::vector<std::shared_ptr<ImgRaster> >::const_iterator imit=begin();

  for(imit=begin();imit!=end();++imit){
  // while((imgReader=getNextImage())){
    if(verbose_opt[0])
      cout << "m_index: " << m_index << endl;
  // for(int iimg=0;iimg<imgReader.size();++iimg){
    // try{
    // imgReader.open(input_opt[iimg],GA_ReadOnly,memory_opt[0]);
    // }
    // catch(string error){
    //   cerr << "Error: could not open file " << input_opt[iimg] << ": " << error << std::endl;
    //   exit(1);
    // }
    if(!isGeoRef)
      isGeoRef=(*imit)->isGeoRef();
    if((*imit)->isGeoRef()&&projection_opt.empty())
      projectionString=(*imit)->getProjection();
    if(dx_opt.empty()){
      if(m_index<=1||(*imit)->getDeltaX()<dx)
        dx=(*imit)->getDeltaX();
    }
    
    if(dy_opt.empty()){
      if(m_index<=1||(*imit)->getDeltaY()<dy)
        dy=(*imit)->getDeltaY();
    }
    if(band_opt.size())
      ncropband+=band_opt.size();
    else
      ncropband+=(*imit)->nrOfBand();
    // (*imit)->close();
  }

  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0])
    cout << "possible output data types: ";
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(verbose_opt[0]){
    cout << endl;
    if(theType==GDT_Unknown)
      cout << "Unknown output pixel type: " << otype_opt[0] << endl;
    else
      cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
  }
  //bounding box of cropped image
  double cropulx=ulx_opt[0];
  double cropuly=uly_opt[0];
  double croplrx=lrx_opt[0];
  double croplry=lry_opt[0];
  //get bounding box from extentReader if defined
  ImgReaderOgr extentReader;

  if(extent_opt.size()){
    double e_ulx;
    double e_uly;
    double e_lrx;
    double e_lry;
    for(int iextent=0;iextent<extent_opt.size();++iextent){
      try{
        extentReader.open(extent_opt[iextent]);
        if(!(extentReader.getExtent(e_ulx,e_uly,e_lrx,e_lry))){
          ostringstream os;
          os << "Error: could not get extent from " << extent_opt[0] << endl;
          throw(os.str());
        }
      }
      catch(string error){
        cerr << error << std::endl;
        exit(1);
      }
      if(!iextent){
	ulx_opt[0]=e_ulx;
	uly_opt[0]=e_uly;
	lrx_opt[0]=e_lrx;
	lry_opt[0]=e_lry;
      }
      else{
	if(e_ulx<ulx_opt[0])
	  ulx_opt[0]=e_ulx;
	if(e_uly>uly_opt[0])
	  uly_opt[0]=e_uly;
	if(e_lrx>lrx_opt[0])
	  lrx_opt[0]=e_lrx;
	if(e_lry<lry_opt[0])
	  lry_opt[0]=e_lry;
      }
      extentReader.close();
    }
    if(croplrx>cropulx&&cropulx>ulx_opt[0])
      ulx_opt[0]=cropulx;
    if(croplrx>cropulx&&croplrx<lrx_opt[0])
      lrx_opt[0]=croplrx;
    if(cropuly>croplry&&cropuly<uly_opt[0])
      uly_opt[0]=cropuly;
    if(croplry<cropuly&&croplry>lry_opt[0])
      lry_opt[0]=croplry;
    if(cut_opt.size()||eoption_opt.size())
      extentReader.open(extent_opt[0]);
  }
  else if(cx_opt.size()&&cy_opt.size()&&nx_opt.size()&&ny_opt.size()){
    ulx_opt[0]=cx_opt[0]-nx_opt[0]/2.0;
    uly_opt[0]=(isGeoRef) ? cy_opt[0]+ny_opt[0]/2.0 : cy_opt[0]-ny_opt[0]/2.0;
    lrx_opt[0]=cx_opt[0]+nx_opt[0]/2.0;
    lry_opt[0]=(isGeoRef) ? cy_opt[0]-ny_opt[0]/2.0 : cy_opt[0]+ny_opt[0]/2.0;
  }
  else if(cx_opt.size()&&cy_opt.size()&&ns_opt.size()&&nl_opt.size()){
    ulx_opt[0]=cx_opt[0]-ns_opt[0]*dx/2.0;
    uly_opt[0]=(isGeoRef) ? cy_opt[0]+nl_opt[0]*dy/2.0 : cy_opt[0]-nl_opt[0]*dy/2.0;
    lrx_opt[0]=cx_opt[0]+ns_opt[0]*dx/2.0;
    lry_opt[0]=(isGeoRef) ? cy_opt[0]-nl_opt[0]*dy/2.0 : cy_opt[0]+nl_opt[0]*dy/2.0;
  }

  if(verbose_opt[0])
    cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;

  int ncropcol=0;
  int ncroprow=0;

  ImgRaster maskReader;
  if(extent_opt.size()&&(cut_opt[0]||eoption_opt.size())){
    if(mask_opt.size()){
      string errorString="Error: can only either mask or extent extent with cutline, not both";
      throw(errorString);
    }
    try{
      ncropcol=abs(static_cast<unsigned int>(ceil((lrx_opt[0]-ulx_opt[0])/dx)));
      ncroprow=abs(static_cast<unsigned int>(ceil((uly_opt[0]-lry_opt[0])/dy)));
      maskReader.open(ncropcol,ncroprow,1,GDT_Float64);
      double gt[6];
      gt[0]=ulx_opt[0];
      gt[1]=dx;
      gt[2]=0;
      gt[3]=uly_opt[0];
      gt[4]=0;
      gt[5]=-dy;
      maskReader.setGeoTransform(gt);
      if(projection_opt.size())
	maskReader.setProjectionProj4(projection_opt[0]);
      else if(projectionString.size())
	maskReader.setProjection(projectionString);
      
      vector<double> burnValues(1,1);//burn value is 1 (single band)
      maskReader.rasterizeBuf(extentReader,burnValues,eoption_opt);
    }
    catch(string error){
      cerr << error << std::endl;
      exit(2);
    }
  }
  else if(mask_opt.size()==1){
    try{
      //there is only a single mask
      maskReader.open(mask_opt[0]);
      if(mskband_opt[0]>=maskReader.nrOfBand()){
	string errorString="Error: illegal mask band";
	throw(errorString);
      }
    }
    catch(string error){
      cerr << error << std::endl;
      exit(2);
    }
  }

  //determine number of output bands
  int writeBand=0;//write band

  if(scale_opt.size()){
    while(scale_opt.size()<band_opt.size())
      scale_opt.push_back(scale_opt[0]);
  }
  if(offset_opt.size()){
    while(offset_opt.size()<band_opt.size())
      offset_opt.push_back(offset_opt[0]);
  }
  if(autoscale_opt.size()){
    assert(autoscale_opt.size()%2==0);
  }

  // for(int iimg=0;iimg<input_opt.size();++iimg){
  for(imit=begin();imit!=end();++imit){
  // for(int iimg=0;iimg<imgReader.size();++iimg){
    // if(verbose_opt[0])
    //   cout << "opening image " << input_opt[iimg] << endl;
    // try{
    //   imgReader.open(input_opt[iimg],GA_ReadOnly,memory_opt[0]);
    // }
    // catch(string error){
    //   cerr << error << std::endl;
    //   exit(2);
    // }
    //if output type not set, get type from input image
    if(theType==GDT_Unknown){
      theType=(*imit)->getDataType();
      if(verbose_opt[0])
        cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
    }
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=(*imit)->getInterleave();
      option_opt.push_back(theInterleave);
    }
    // if(verbose_opt[0])
    //   cout << "size of " << input_opt[iimg] << ": " << ncol << " cols, "<< nrow << " rows" << endl;
    double uli,ulj,lri,lrj;//image coordinates
    bool forceEUgrid=false;
    if(projection_opt.size())
      forceEUgrid=(!(projection_opt[0].compare("EPSG:3035"))||!(projection_opt[0].compare("EPSG:3035"))||projection_opt[0].find("ETRS-LAEA")!=string::npos);
    if(ulx_opt[0]>=lrx_opt[0]){//default bounding box: no cropping
      uli=0;
      lri=(*imit)->nrOfCol()-1;
      ulj=0;
      lrj=(*imit)->nrOfRow()-1;
      ncropcol=(*imit)->nrOfCol();
      ncroprow=(*imit)->nrOfRow();
      (*imit)->getBoundingBox(cropulx,cropuly,croplrx,croplry);
      double magicX=1,magicY=1;
      // (*imit)->getMagicPixel(magicX,magicY);
      if(forceEUgrid){
	//force to LAEA grid
	Egcs egcs;
        egcs.setLevel(egcs.res2level(dx));
	egcs.force2grid(cropulx,cropuly,croplrx,croplry);
	(*imit)->geo2image(cropulx+(magicX-1.0)*(*imit)->getDeltaX(),cropuly-(magicY-1.0)*(*imit)->getDeltaY(),uli,ulj);
	(*imit)->geo2image(croplrx+(magicX-2.0)*(*imit)->getDeltaX(),croplry-(magicY-2.0)*(*imit)->getDeltaY(),lri,lrj);
      }
      (*imit)->geo2image(cropulx+(magicX-1.0)*(*imit)->getDeltaX(),cropuly-(magicY-1.0)*(*imit)->getDeltaY(),uli,ulj);
      (*imit)->geo2image(croplrx+(magicX-2.0)*(*imit)->getDeltaX(),croplry-(magicY-2.0)*(*imit)->getDeltaY(),lri,lrj);
      ncropcol=abs(static_cast<unsigned int>(ceil((croplrx-cropulx)/dx)));
      ncroprow=abs(static_cast<unsigned int>(ceil((cropuly-croplry)/dy)));
    }
    else{
      double magicX=1,magicY=1;
      // (*imit)->getMagicPixel(magicX,magicY);
      cropulx=ulx_opt[0];
      cropuly=uly_opt[0];
      croplrx=lrx_opt[0];
      croplry=lry_opt[0];
      if(forceEUgrid){
	//force to LAEA grid
	Egcs egcs;
        egcs.setLevel(egcs.res2level(dx));
	egcs.force2grid(cropulx,cropuly,croplrx,croplry);
      }
      else if(align_opt[0]){
      	if(cropulx>(*imit)->getUlx())
      	  cropulx-=fmod(cropulx-(*imit)->getUlx(),dx);
      	else if(cropulx<(*imit)->getUlx())
      	  cropulx+=fmod((*imit)->getUlx()-cropulx,dx)-dx;
      	if(croplrx<(*imit)->getLrx())
      	  croplrx+=fmod((*imit)->getLrx()-croplrx,dx);
      	else if(croplrx>(*imit)->getLrx())
      	  croplrx-=fmod(croplrx-(*imit)->getLrx(),dx)+dx;
      	if(croplry>(*imit)->getLry())
      	  croplry-=fmod(croplry-(*imit)->getLry(),dy);
      	else if(croplry<(*imit)->getLry())
      	  croplry+=fmod((*imit)->getLry()-croplry,dy)-dy;
      	if(cropuly<(*imit)->getUly())
      	  cropuly+=fmod((*imit)->getUly()-cropuly,dy);
      	else if(cropuly>(*imit)->getUly())
      	  cropuly-=fmod(cropuly-(*imit)->getUly(),dy)+dy;
      }
      (*imit)->geo2image(cropulx+(magicX-1.0)*(*imit)->getDeltaX(),cropuly-(magicY-1.0)*(*imit)->getDeltaY(),uli,ulj);
      (*imit)->geo2image(croplrx+(magicX-2.0)*(*imit)->getDeltaX(),croplry-(magicY-2.0)*(*imit)->getDeltaY(),lri,lrj);

      ncropcol=abs(static_cast<unsigned int>(ceil((croplrx-cropulx)/dx)));
      ncroprow=abs(static_cast<unsigned int>(ceil((cropuly-croplry)/dy)));
      uli=floor(uli);
      ulj=floor(ulj);
      lri=floor(lri);
      lrj=floor(lrj);
    }

    // double deltaX=(*imit)->getDeltaX();
    // double deltaY=(*imit)->getDeltaY();
    if(!imgWriter->nrOfBand()){//not opened yet
      if(verbose_opt[0]){
	cout << "cropulx: " << cropulx << endl;
	cout << "cropuly: " << cropuly << endl;
	cout << "croplrx: " << croplrx << endl;
	cout << "croplry: " << croplry << endl;
	cout << "ncropcol: " << ncropcol << endl;
	cout << "ncroprow: " << ncroprow << endl;
	cout << "cropulx+ncropcol*dx: " << cropulx+ncropcol*dx << endl;
	cout << "cropuly-ncroprow*dy: " << cropuly-ncroprow*dy << endl;
	cout << "upper left column of input image: " << uli << endl;
	cout << "upper left row of input image: " << ulj << endl;
	cout << "lower right column of input image: " << lri << endl;
	cout << "lower right row of input image: " << lrj << endl;
	cout << "new number of cols: " << ncropcol << endl;
	cout << "new number of rows: " << ncroprow << endl;
	cout << "new number of bands: " << ncropband << endl;
      }
      string imageType;//=(*imit)->getImageType();
      if(oformat_opt.size())//default
        imageType=oformat_opt[0];
      try{
        imgWriter->open(ncropcol,ncroprow,ncropband,theType);
        imgWriter->setNoData(nodata_opt);
        // imgWriter->open(output_opt[0],ncropcol,ncroprow,ncropband,theType,imageType,memory_opt[0],option_opt);
	// if(nodata_opt.size()){
	//   imgWriter->setNoData(nodata_opt);
	// }
      }
      catch(string errorstring){
        cout << errorstring << endl;
        exit(4);
      }
      if(description_opt.size())
	imgWriter->setImageDescription(description_opt[0]);
      double gt[6];
      gt[0]=cropulx;
      gt[1]=dx;
      gt[2]=0;
      gt[3]=cropuly;
      gt[4]=0;
      gt[5]=((*imit)->isGeoRef())? -dy : dy;
      imgWriter->setGeoTransform(gt);
      if(projection_opt.size()){
	if(verbose_opt[0])
	  cout << "projection: " << projection_opt[0] << endl;
	imgWriter->setProjectionProj4(projection_opt[0]);
      }
      else
	imgWriter->setProjection((*imit)->getProjection());
      if(imgWriter->getDataType()==GDT_Byte){
	if(colorTable_opt.size()){
	  if(colorTable_opt[0]!="none")
	    imgWriter->setColorTable(colorTable_opt[0]);
	}
	else if ((*imit)->getColorTable()!=NULL)//copy colorTable from input image
	  imgWriter->setColorTable((*imit)->getColorTable());
      }
    }

    double startCol=uli;
    double endCol=lri;
    if(uli<0)
      startCol=0;
    else if(uli>=(*imit)->nrOfCol())
      startCol=(*imit)->nrOfCol()-1;
    if(lri<0)
      endCol=0;
    else if(lri>=(*imit)->nrOfCol())
      endCol=(*imit)->nrOfCol()-1;
    double startRow=ulj;
    double endRow=lrj;
    if(ulj<0)
      startRow=0;
    else if(ulj>=(*imit)->nrOfRow())
      startRow=(*imit)->nrOfRow()-1;
    if(lrj<0)
      endRow=0;
    else if(lrj>=(*imit)->nrOfRow())
      endRow=(*imit)->nrOfRow()-1;

    vector<double> readBuffer;
    unsigned int nband=(band_opt.size())?band_opt.size() : (*imit)->nrOfBand();
    for(unsigned int iband=0;iband<nband;++iband){
      unsigned int readBand=(band_opt.size()>iband)?band_opt[iband]:iband;
      if(verbose_opt[0]){
	cout << "extracting band " << readBand << endl;
	pfnProgress(progress,pszMessage,pProgressArg);
      }
      double theMin=0;
      double theMax=0;
      if(autoscale_opt.size()){
	try{
	  (*imit)->getMinMax(static_cast<unsigned int>(startCol),static_cast<unsigned int>(endCol),static_cast<unsigned int>(startRow),static_cast<unsigned int>(endRow),readBand,theMin,theMax);
	}
	catch(string errorString){
	  cout << errorString << endl;
	}
	if(verbose_opt[0])
	  cout << "minmax: " << theMin << ", " << theMax << endl;
	double theScale=(autoscale_opt[1]-autoscale_opt[0])/(theMax-theMin);
	double theOffset=autoscale_opt[0]-theScale*theMin;
	(*imit)->setScale(theScale,readBand);
	(*imit)->setOffset(theOffset,readBand);
      }	
      else{
	if(scale_opt.size()){
	  if(scale_opt.size()>iband)
	    (*imit)->setScale(scale_opt[iband],readBand);
	  else
	    (*imit)->setScale(scale_opt[0],readBand);
	}
	if(offset_opt.size()){
	  if(offset_opt.size()>iband)
	    (*imit)->setOffset(offset_opt[iband],readBand);
	  else
	    (*imit)->setOffset(offset_opt[0],readBand);
	}
      }

      double readRow=0;
      double readCol=0;
      double lowerCol=0;
      double upperCol=0;
      for(int irow=0;irow<imgWriter->nrOfRow();++irow){
	vector<double> lineMask;
	double x=0;
	double y=0;
	//convert irow to geo
	imgWriter->image2geo(0,irow,x,y);
	//lookup corresponding row for irow in this file
	(*imit)->geo2image(x,y,readCol,readRow);
	vector<double> writeBuffer;
	if(readRow<0||readRow>=(*imit)->nrOfRow()){
	  for(int icol=0;icol<imgWriter->nrOfCol();++icol)
	    writeBuffer.push_back(nodataValue);
	}
	else{
	  try{
            if(endCol<(*imit)->nrOfCol()-1){
              (*imit)->readData(readBuffer,startCol,endCol+1,readRow,readBand,theResample);
            }
            else{
              (*imit)->readData(readBuffer,startCol,endCol,readRow,readBand,theResample);
            }
	    double oldRowMask=-1;//keep track of row mask to optimize number of line readings
	    for(int icol=0;icol<imgWriter->nrOfCol();++icol){
	      imgWriter->image2geo(icol,irow,x,y);
	      //lookup corresponding row for irow in this file
	      (*imit)->geo2image(x,y,readCol,readRow);
	      if(readCol<0||readCol>=(*imit)->nrOfCol()){
		writeBuffer.push_back(nodataValue);
	      }
	      else{
                bool valid=true;
		double geox=0;
		double geoy=0;
                if(maskReader.isInit()){
		  //read mask
		  double colMask=0;
		  double rowMask=0;

		  imgWriter->image2geo(icol,irow,geox,geoy);
		  maskReader.geo2image(geox,geoy,colMask,rowMask);
		  colMask=static_cast<unsigned int>(colMask);
		  rowMask=static_cast<unsigned int>(rowMask);
		  if(rowMask>=0&&rowMask<maskReader.nrOfRow()&&colMask>=0&&colMask<maskReader.nrOfCol()){
		    if(static_cast<unsigned int>(rowMask)!=static_cast<unsigned int>(oldRowMask)){

		      try{
			maskReader.readData(lineMask,static_cast<unsigned int>(rowMask),mskband_opt[0]);
		      }
		      catch(string errorstring){
			cerr << errorstring << endl;
			exit(1);
		      }
		      catch(...){
			cerr << "error caught" << std::endl;
			exit(3);
		      }
		      oldRowMask=rowMask;
		    }
                    for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
                      if(lineMask[colMask]==msknodata_opt[ivalue]){
                        valid=false;
			if(nodata_opt.size()>ivalue)
			  nodataValue=nodata_opt[ivalue];
                      }
                    }
		  }
		}
                if(!valid)
                  writeBuffer.push_back(nodataValue);
                else{
                  switch(theResample){
                  case(BILINEAR):
                    lowerCol=readCol-0.5;
                    lowerCol=static_cast<unsigned int>(lowerCol);
                    upperCol=readCol+0.5;
                    upperCol=static_cast<unsigned int>(upperCol);
                    if(lowerCol<0)
                      lowerCol=0;
                    if(upperCol>=(*imit)->nrOfCol())
                      upperCol=(*imit)->nrOfCol()-1;
                    writeBuffer.push_back((readCol-0.5-lowerCol)*readBuffer[upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[lowerCol-startCol]);
                    break;
                  default:
                    readCol=static_cast<unsigned int>(readCol);
                    readCol-=startCol;//we only start reading from startCol
                    writeBuffer.push_back(readBuffer[readCol]);
                    break;
		  }
		}
	      }
	    }
	  }
	  catch(string errorstring){
	    cout << errorstring << endl;
	    exit(2);
	  }
	}
	if(writeBuffer.size()!=imgWriter->nrOfCol())
	  cout << "writeBuffer.size()=" << writeBuffer.size() << ", imgWriter->nrOfCol()=" << imgWriter->nrOfCol() << endl;

	assert(writeBuffer.size()==imgWriter->nrOfCol());
	try{
	  imgWriter->writeData(writeBuffer,irow,writeBand);
	}
	catch(string errorstring){
	  cout << errorstring << endl;
	  exit(3);
	}
	if(verbose_opt[0]){
	  progress=(1.0+irow);
	  progress/=imgWriter->nrOfRow();
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
	else{
	  progress=(1.0+irow);
	  progress+=(imgWriter->nrOfRow()*writeBand);
	  progress/=imgWriter->nrOfBand()*imgWriter->nrOfRow();
	  assert(progress>=0);
	  assert(progress<=1);
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
      }
      ++writeBand;
    }
    // (*imit)->close();
  }
  if(extent_opt.size()&&(cut_opt[0]||eoption_opt.size())){
    extentReader.close();
  }
  if(maskReader.isInit())
    maskReader.close();
  return(0);
}

/**
 * @param imgRaster output raster composite dataset
 * @return output image
 **/
shared_ptr<ImgRaster> ImgCollection::composite(const AppFactory& app){
  std::vector<std::shared_ptr<ImgRaster> >::const_iterator imit=begin();
  //test
  cout << "typeid of first element in collection: " << typeid(*imit).name() << endl;
  shared_ptr<ImgRaster> imgWriter=(*imit)->clone();//create clone to first object, allowing for polymorphism in case of derived ImgRaster objects
  //test
  cout << "typeid of imgWriter: " << typeid(imgWriter).name() << endl;
  composite(imgWriter, app);
  return(imgWriter);
}

/**
 * @param imgRaster output raster composite dataset
 * @return 0 if successful
 **/
int ImgCollection::composite(shared_ptr<ImgRaster> imgWriter, const AppFactory& app){
  Optionpk<unsigned int>  band_opt("b", "band", "band index(es) to crop (leave empty if all bands must be retained)");
  Optionpk<double>  dx_opt("dx", "dx", "Output resolution in x (in meter) (empty: keep original resolution)");
  Optionpk<double>  dy_opt("dy", "dy", "Output resolution in y (in meter) (empty: keep original resolution)");
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file");
  Optionpk<bool> cut_opt("cut", "crop_to_cutline", "Crop the extent of the target dataset to the extent of the cutline.",false);
  Optionpk<string> eoption_opt("eo","eo", "special extent options controlling rasterization: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG, e.g., -eo ATTRIBUTE=fieldname");
  Optionpk<string> mask_opt("m", "mask", "Use the first band of the specified file as a validity mask (0 is nodata).");
  Optionpk<unsigned int> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed)", 0);
  Optionpk<float> msknodata_opt("msknodata", "msknodata", "Mask value not to consider for composite.", 0);
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box", 0.0);
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box", 0.0);
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box", 0.0);
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box", 0.0);
  Optionpk<string> crule_opt("cr", "crule", "Composite rule (overwrite, maxndvi, maxband, minband, mean, mode (only for byte images), median, sum, maxallbands, minallbands, stdev", "overwrite");
  Optionpk<unsigned int> ruleBand_opt("cb", "cband", "band index used for the composite rule (e.g., for ndvi, use --cband=0 --cband=1 with 0 and 1 indices for red and nir band respectively", 0);
  Optionpk<double> srcnodata_opt("srcnodata", "srcnodata", "invalid value(s) for input raster dataset");
  Optionpk<unsigned int> bndnodata_opt("bndnodata", "bndnodata", "Band(s) in input image to check if pixel is valid (used for srcnodata, min and max options)", 0);
  Optionpk<double> minValue_opt("min", "min", "flag values smaller or equal to this value as invalid.");
  Optionpk<double> maxValue_opt("max", "max", "flag values larger or equal to this value as invalid.");
  Optionpk<double>  dstnodata_opt("dstnodata", "dstnodata", "nodata value to put in output raster dataset if not valid or out of bounds.", 0);
  Optionpk<string>  resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation).", "near");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the spatial reference for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<short> file_opt("file", "file", "write number of observations (1) or sequence nr of selected file (2) for each pixels as additional layer in composite", 0);
  Optionpk<short> weight_opt("w", "weight", "Weights (type: short) for the composite, use one weight for each input file in same order as input files are provided). Use value 1 for equal weights.", 1);
  Optionpk<short> class_opt("c", "class", "classes for multi-band output image: each band represents the number of observations for one specific class. Use value 0 for no multi-band output image.", 0);
  Optionpk<string>  colorTable_opt("ct", "ct", "color table file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<string>  description_opt("d", "description", "Set image description");
  Optionpk<bool>  align_opt("align", "align", "Align output bounding box to input image",false);
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);
  Optionpk<short>  verbose_opt("v", "verbose", "verbose", 0,2);

  extent_opt.setHide(1);
  cut_opt.setHide(1);
  eoption_opt.setHide(1);
  mask_opt.setHide(1);
  msknodata_opt.setHide(1);
  mskband_opt.setHide(1);
  file_opt.setHide(1);
  weight_opt.setHide(1);
  class_opt.setHide(1);
  colorTable_opt.setHide(1);
  description_opt.setHide(1);
  memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=band_opt.retrieveOption(app.getArgc(),app.getArgv());
    dx_opt.retrieveOption(app.getArgc(),app.getArgv());
    dy_opt.retrieveOption(app.getArgc(),app.getArgv());
    extent_opt.retrieveOption(app.getArgc(),app.getArgv());
    cut_opt.retrieveOption(app.getArgc(),app.getArgv());
    eoption_opt.retrieveOption(app.getArgc(),app.getArgv());
    mask_opt.retrieveOption(app.getArgc(),app.getArgv());
    msknodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    mskband_opt.retrieveOption(app.getArgc(),app.getArgv());
    ulx_opt.retrieveOption(app.getArgc(),app.getArgv());
    uly_opt.retrieveOption(app.getArgc(),app.getArgv());
    lrx_opt.retrieveOption(app.getArgc(),app.getArgv());
    lry_opt.retrieveOption(app.getArgc(),app.getArgv());
    crule_opt.retrieveOption(app.getArgc(),app.getArgv());
    ruleBand_opt.retrieveOption(app.getArgc(),app.getArgv());
    srcnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    bndnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    minValue_opt.retrieveOption(app.getArgc(),app.getArgv());
    maxValue_opt.retrieveOption(app.getArgc(),app.getArgv());
    dstnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    resample_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    projection_opt.retrieveOption(app.getArgc(),app.getArgv());
    file_opt.retrieveOption(app.getArgc(),app.getArgv());
    weight_opt.retrieveOption(app.getArgc(),app.getArgv());
    class_opt.retrieveOption(app.getArgc(),app.getArgv());
    colorTable_opt.retrieveOption(app.getArgc(),app.getArgv());
    description_opt.retrieveOption(app.getArgc(),app.getArgv());
    align_opt.retrieveOption(app.getArgc(),app.getArgv());
    memory_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    std::ostringstream helpStream;
    helpStream << "Usage: pkcomposite -i input [-i input]* -o output" << std::endl;
    helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    throw(helpStream.str());//help was invoked, stop processing
  }
  if(empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(0);
  }
  // ///band index(es) to crop (leave empty if all bands must be retained)
  // vector<int>  band_opt;
  // ///Output resolution in x (in meter) (empty: keep original resolution)
  // vector<double>  dx_opt;
  // ///Output resolution in y (in meter) (empty: keep original resolution)
  // vector<double>  dy_opt;
  // ///get boundary from extent from polygons in vector file
  // vector<string>  extent_opt;
  // ///Crop the extent of the target dataset to the extent of the cutline.
  // vector<bool> cut_opt;
  // ///special extent options controlling rasterization: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG, e.g., -eo ATTRIBUTE=fieldname
  // vector<string> eoption_opt;
  // ///Use the first band of the specified file as a validity mask (0 is nodata).
  // vector<string> mask_opt;
  // ///Mask value not to consider for composite.
  // vector<float> msknodata_opt;
  // ///Mask band to read (0 indexed)
  // vector<short> mskband_opt;
  // ///Upper left x value bounding box
  // vector<double>  ulx_opt;
  // ///Upper left y value bounding box
  // vector<double>  uly_opt;
  // ///Lower right x value bounding box
  // vector<double>  lrx_opt;
  // ///Lower right y value bounding box
  // vector<double>  lry_opt;
  // ///Composite rule (overwrite, maxndvi, maxband, minband, mean, mode (only for byte images), median, sum, maxallbands, minallbands, stdev
  // vector<string> crule_opt;
  // ///band index used for the composite rule (e.g., for ndvi, use --cband=0 --cband=1 with 0 and 1 indices for red and nir band respectively
  // vector<int> ruleBand_opt;
  // ///invalid value(s) for input raster dataset
  // vector<double> srcnodata_opt;
  // ///Band(s) in input image to check if pixel is valid (used for srcnodata, min and max options)
  // vector<int> bndnodata_opt;
  // ///flag values smaller or equal to this value as invalid.
  // vector<double> minValue_opt;
  // ///flag values larger or equal to this value as invalid.
  // vector<double> maxValue_opt;
  // ///nodata value to put in output raster dataset if not valid or out of bounds.
  // vector<double>  dstnodata_opt;
  // ///Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation).
  // vector<string>  resample_opt("near",1);
  // ///Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image
  // vector<string>  otype_opt("",1);
  // ///Creation option for output file. Multiple options can be specified.
  // vector<string> option_opt;
  // ///Override the spatial reference for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid
  // vector<string>  projection_opt;
  // ///write number of observations (1) or sequence nr of selected file (2) for each pixels as additional layer in composite
  // vector<short> file_opt;
  // ///Weights (type: short) for the composite, use one weight for each input file in same order as input files are provided). Use value 1 for equal weights.
  // vector<short> weight_opt;
  // ///classes for multi-band output image: each band represents the number of observations for one specific class. Use value 0 for no multi-band output image.
  // vector<short> class_opt;
  // ///color table file with 5 columns: id R G B ALFA (0: transparent, 255: solid)
  // vector<string>  colorTable_opt;
  // ///Set image description
  // vector<string>  description_opt;
  // ///Align output bounding box to input image
  // vector<bool>  align_opt;
  // ///verbose
  // vector<short> verbose_opt(0,1);

  std::map<std::string, CRULE_TYPE> cruleMap;
  // //initialize cruleMap
  // enum CRULE_TYPE {overwrite=0, maxndvi=1, maxband=2, minband=3, validband=4, mean=5, mode=6, median=7,sum=8};
  
  cruleMap["overwrite"]=overwrite;
  cruleMap["maxndvi"]=maxndvi;
  cruleMap["maxband"]=maxband;
  cruleMap["minband"]=minband;
  cruleMap["validband"]=validband;
  cruleMap["mean"]=mean;
  cruleMap["mode"]=mode;
  cruleMap["median"]=median;
  cruleMap["sum"]=sum;
  cruleMap["maxallbands"]=maxallbands;
  cruleMap["minallbands"]=minallbands;
  cruleMap["stdev"]=stdev;

  if(srcnodata_opt.size()){
    while(srcnodata_opt.size()<bndnodata_opt.size())
      srcnodata_opt.push_back(srcnodata_opt[0]);
  }
  while(bndnodata_opt.size()<srcnodata_opt.size())
    bndnodata_opt.push_back(bndnodata_opt[0]);
  if(minValue_opt.size()){
    while(minValue_opt.size()<bndnodata_opt.size())
      minValue_opt.push_back(minValue_opt[0]);
    while(bndnodata_opt.size()<minValue_opt.size())
      bndnodata_opt.push_back(bndnodata_opt[0]);
  }
  if(maxValue_opt.size()){
    while(maxValue_opt.size()<bndnodata_opt.size())
      maxValue_opt.push_back(maxValue_opt[0]);
    while(bndnodata_opt.size()<maxValue_opt.size())
      bndnodata_opt.push_back(bndnodata_opt[0]);
  }

  RESAMPLE theResample;
  if(resample_opt[0]=="near"){
    theResample=NEAR;
    if(verbose_opt[0])
      cout << "resampling: nearest neighbor" << endl;
  }
  else if(resample_opt[0]=="bilinear"){
    theResample=BILINEAR;
    if(verbose_opt[0])
      cout << "resampling: bilinear interpolation" << endl;
  }
  else{
    std::cout << "Error: resampling method " << resample_opt[0] << " not supported" << std::endl;
    exit(1);
  }

  int nband=0;
  int nwriteBand=0;
  vector<short> bands;

  //get bounding box
  double maxLRX=lrx_opt[0];
  double maxULY=uly_opt[0];
  double minULX=ulx_opt[0];
  double minLRY=lry_opt[0];
  double magic_x=1,magic_y=1;//magic pixel for GDAL map info

  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0])
    cout << "possible output data types: ";
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(verbose_opt[0]){
    cout << endl;
    if(theType==GDT_Unknown)
      cout << "Unknown output pixel type: " << otype_opt[0] << endl;
    else
      cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
  }

  double dx=0;
  double dy=0;
  //get bounding box from extentReader if defined
  ImgReaderOgr extentReader;
  if(extent_opt.size()){
    double e_ulx;
    double e_uly;
    double e_lrx;
    double e_lry;
    for(int iextent=0;iextent<extent_opt.size();++iextent){
      extentReader.open(extent_opt[iextent]);
      if(!(extentReader.getExtent(e_ulx,e_uly,e_lrx,e_lry))){
        cerr << "Error: could not get extent from " << extent_opt[0] << endl;
        exit(1);
      }
      if(!iextent){
	ulx_opt[0]=e_ulx;
	uly_opt[0]=e_uly;
	lrx_opt[0]=e_lrx;
	lry_opt[0]=e_lry;
      }
      else{
	if(e_ulx<ulx_opt[0])
	  ulx_opt[0]=e_ulx;
	if(e_uly>uly_opt[0])
	  uly_opt[0]=e_uly;
	if(e_lrx>lrx_opt[0])
	  lrx_opt[0]=e_lrx;
	if(e_lry<lry_opt[0])
	  lry_opt[0]=e_lry;
      }
      extentReader.close();
    }
    if(maxLRX>minULX&&minULX>ulx_opt[0])
      ulx_opt[0]=minULX;
    if(maxLRX>minULX&&maxLRX<lrx_opt[0])
      lrx_opt[0]=maxLRX;
    if(maxULY>minLRY&&maxULY<uly_opt[0])
      uly_opt[0]=maxULY;
    if(minLRY<maxULY&&minLRY>lry_opt[0])
      lry_opt[0]=minLRY;
    if(cut_opt.size()||eoption_opt.size())
      extentReader.open(extent_opt[0]);
  }

  if(verbose_opt[0])
    cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;


  string theProjection="";
  // GDALColorTable* theColorTable=NULL;
  bool init=false;

  std::vector<std::shared_ptr<ImgRaster> >::const_iterator imit=begin();
  for(imit=begin();imit!=end();++imit){
  // for(int ifile=0;ifile<(*imit).size();++ifile){
    //todo: must be in init part only?
    if(colorTable_opt.empty())
      // if((*imit)->getColorTable())
      //   theColorTable=((*imit)->getColorTable()->Clone());
    if(projection_opt.empty())
      theProjection=(*imit)->getProjection();
    // if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    //   string theInterleave="INTERLEAVE=";
    //   theInterleave+=(*imit)->getInterleave();
    //   option_opt.push_back(theInterleave);
    // }

    if((ulx_opt[0]||uly_opt[0]||lrx_opt[0]||lry_opt[0])&&(!(*imit)->covers(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
      if(verbose_opt[0])
	cout << "Input not within bounding box, skipping..." << endl;
      continue;
    }
    double theULX, theULY, theLRX, theLRY;
    (*imit)->getBoundingBox(theULX,theULY,theLRX,theLRY);
    if(theLRY>theULY){
      cerr << "Error: input is not georeferenced, only referenced images are supported for pkcomposite " << endl;
      exit(1);
    }
    if(verbose_opt[0])
      cout << "Bounding Box (ULX ULY LRX LRY): " << fixed << setprecision(6) << theULX << " " << theULY << " " << theLRX << " " << theLRY << endl;
    if(!init){
      if(verbose_opt[0]){
        switch(cruleMap[crule_opt[0]]){
        default:
        case(overwrite):
          cout << "Composite rule: overwrite" << endl;
          break;
        case(maxndvi):
          cout << "Composite rule: max ndvi" << endl;
          break;
        case(maxband):
          cout << "Composite rule: max band" << endl;
          break;
        case(minband):
          cout << "Composite rule: min band" << endl;
          break;
        case(validband):
          cout << "Composite rule: valid band" << endl;
          break;
        case(mean):
          cout << "Composite rule: mean value" << endl;
          break;
        case(mode):
          cout << "Composite rule: max voting (only for byte images)" << endl;
          break;
        case(median):
          cout << "Composite rule: median" << endl;
          break;
        case(stdev):
          cout << "Composite rule: stdev" << endl;
          break;
        case(sum):
          cout << "Composite rule: sum" << endl;
          break;
        case(minallbands):
          cout << "Composite rule: minallbands" << endl;
          break;
        case(maxallbands):
          cout << "Composite rule: maxallbands" << endl;
          break;
        }
      }
      if(band_opt.size()){
	nband=band_opt.size();
        bands.resize(band_opt.size());
        for(unsigned int iband=0;iband<band_opt.size();++iband){
          bands[iband]=band_opt[iband];
          assert(bands[iband]<(*imit)->nrOfBand());
        }
      }
      else{
	nband=(*imit)->nrOfBand();
        bands.resize(nband);
        for(unsigned int iband=0;iband<nband;++iband)
          bands[iband]=iband;
      }
      for(unsigned int iband=0;iband<bndnodata_opt.size();++iband){
        assert(bndnodata_opt[iband]>=0&&bndnodata_opt[iband]<nband);
      }
      //if output type not set, get type from input image
      if(theType==GDT_Unknown){
        theType=(*imit)->getDataType();
        if(verbose_opt[0])
          cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
      }

      if(verbose_opt[0]){
        cout << "input data type: " << theType << endl;
        cout << "nband: " << nband << endl;
      }
      
      maxLRX=theLRX;
      maxULY=theULY;
      minULX=theULX;
      minLRY=theLRY;
      if(dx_opt.size())
	dx=dx_opt[0];
      else
        dx=(*imit)->getDeltaX();
      if(dy_opt.size())
	dy=dy_opt[0];
      else
        dy=(*imit)->getDeltaY();
      init=true;
    }
    else{
      maxLRX=(theLRX>maxLRX)?theLRX:maxLRX;
      maxULY=(theULY>maxULY)?theULY:maxULY;
      minULX=(theULX<minULX)?theULX:minULX;
      minLRY=(theLRY<minLRY)?theLRY:minLRY;
    }
    // (*imit).close();
  }
  if(verbose_opt[0])
    cout << "bounding box input images (ULX ULY LRX LRY): " << fixed << setprecision(6) << minULX << " " << maxULY << " " << maxLRX << " " << minLRY << endl;
  if(ulx_opt[0]||uly_opt[0]||lrx_opt[0]||lry_opt[0]){
    maxLRX=lrx_opt[0];
    maxULY=uly_opt[0];
    minULX=ulx_opt[0];
    minLRY=lry_opt[0];
  }
  
  bool forceEUgrid=false;
  if(projection_opt.size())
    forceEUgrid=(!(projection_opt[0].compare("EPSG:3035"))||!(projection_opt[0].compare("EPSG:3035"))||projection_opt[0].find("ETRS-LAEA")!=string::npos);
  if(forceEUgrid){
    //force to LAEA grid
    minULX=floor(minULX);
    minULX-=static_cast<unsigned int>(minULX)%(static_cast<unsigned int>(dx));
    maxULY=ceil(maxULY);
    if(static_cast<unsigned int>(maxULY)%static_cast<unsigned int>(dy))
      maxULY+=dy;
    maxULY-=static_cast<unsigned int>(maxULY)%(static_cast<unsigned int>(dy));
    maxLRX=ceil(maxLRX);
    if(static_cast<unsigned int>(maxLRX)%static_cast<unsigned int>(dx))
      maxLRX+=dx;
    maxLRX-=static_cast<unsigned int>(maxLRX)%(static_cast<unsigned int>(dx));
    minLRY=floor(minLRY);
    minLRY-=static_cast<unsigned int>(minLRY)%(static_cast<unsigned int>(dy));
  }
  else if(align_opt[0]){
    if(minULX>front()->getUlx())
      minULX-=fmod(minULX-front()->getUlx(),dx);
    else if(minULX<front()->getUlx())
      minULX+=fmod(front()->getUlx()-minULX,dx)-dx;
    if(maxLRX<front()->getLrx())
      maxLRX+=fmod(front()->getLrx()-maxLRX,dx);
    else if(maxLRX>front()->getLrx())
      maxLRX-=fmod(maxLRX-front()->getLrx(),dx)+dx;
    if(minLRY>front()->getLry())
      minLRY-=fmod(minLRY-front()->getLry(),dy);
    else if(minLRY<front()->getLry())
      minLRY+=fmod(front()->getLry()-minLRY,dy)-dy;
    if(maxULY<front()->getUly())
      maxULY+=fmod(front()->getUly()-maxULY,dy);
    else if(maxULY>front()->getUly())
      maxULY-=fmod(maxULY-front()->getUly(),dy)+dy;
  }

  if(verbose_opt[0])
    cout << "bounding box composite image (ULX ULY LRX LRY): " << fixed << setprecision(6) << minULX << " " << maxULY << " " << maxLRX << " " << minLRY << endl;
  //initialize image
  if(verbose_opt[0])
    cout << "initializing composite image..." << endl;
//   double dcol=(maxLRX-minULX+dx-1)/dx;
//   double drow=(maxULY-minLRY+dy-1)/dy;
//   int ncol=static_cast<unsigned int>(dcol);
//   int nrow=static_cast<unsigned int>(drow);

  int ncol=ceil((maxLRX-minULX)/dx);
  int nrow=ceil((maxULY-minLRY)/dy);

  if(verbose_opt[0])
    cout << "composite image dim (nrow x ncol): " << nrow << " x " << ncol << endl;
  while(weight_opt.size()<size())
    weight_opt.push_back(weight_opt[0]);
  if(verbose_opt[0]){
    std::cout << weight_opt << std::endl;
  }
  if(cruleMap[crule_opt[0]]==mode){
    nwriteBand=(file_opt[0])? class_opt.size()+1:class_opt.size();
  }
  else
    nwriteBand=(file_opt[0])? bands.size()+1:bands.size();

  try{
    imgWriter->open(ncol,nrow,nwriteBand,theType);
    imgWriter->setNoData(dstnodata_opt);
  }
  catch(string error){
    cout << error << endl;
  }
  double gt[6];
  gt[0]=minULX;
  gt[1]=dx;
  gt[2]=0;
  gt[3]=maxULY;
  gt[4]=0;
  gt[5]=-dy;
  imgWriter->setGeoTransform(gt);

  if(projection_opt.size()){
    if(verbose_opt[0])
      cout << "projection: " << projection_opt[0] << endl;
    imgWriter->setProjectionProj4(projection_opt[0]);
  }
  else if(theProjection!=""){
    if(verbose_opt[0])
      cout << "projection: " << theProjection << endl;
    imgWriter->setProjection(theProjection);
  }

  //test
  std::cout << "creating mask" << std::endl;
  ImgRaster maskReader;
  //test
  std::cout << "end creating mask" << std::endl;
  if(extent_opt.size()&&(cut_opt[0]||eoption_opt.size())){
    if(mask_opt.size()){
      string errorString="Error: can only either mask or extent extent with cutline, not both";
      throw(errorString);
    }
    try{
      maskReader.open(ncol,nrow,1,GDT_Float64);
      double gt[6];
      gt[0]=minULX;
      gt[1]=dx;
      gt[2]=0;
      gt[3]=maxULY;
      gt[4]=0;
      gt[5]=-dy;
      maskReader.setGeoTransform(gt);
      if(projection_opt.size())
	maskReader.setProjectionProj4(projection_opt[0]);
      else if(theProjection!=""){
	if(verbose_opt[0])
	  cout << "projection: " << theProjection << endl;
	maskReader.setProjection(theProjection);
      }
      vector<double> burnValues(1,1);//burn value is 1 (single band)
      maskReader.rasterizeOgr(extentReader,burnValues,eoption_opt);
    }
    catch(string error){
      cerr << error << std::endl;
      exit(2);
    }
    catch(...){
      cerr << "error caught" << std::endl;
      exit(1);
    }
    //todo: support multiple masks
  }
  else if(mask_opt.size()==1){
    try{
      //there is only a single mask
      maskReader.open(mask_opt[0]);
      if(mskband_opt[0]>=maskReader.nrOfBand()){
	string errorString="Error: illegal mask band";
	throw(errorString);
      }
    }
    catch(string error){
      cerr << error << std::endl;
      exit(2);
    }
  }

  //create composite image
  if(verbose_opt[0])
     cout << "creating composite image" << endl;
  Vector2d<double> writeBuffer(nband,imgWriter->nrOfCol());
  vector<short> fileBuffer(ncol);//holds the number of used files
  Vector2d<short> maxBuffer;//buffer used for maximum voting
  // Vector2d<double> readBuffer(nband);
  vector<Vector2d<unsigned short> > readBuffer(size());
  for(int ifile=0;ifile<size();++ifile)
    readBuffer[ifile].resize((this->at(ifile))->nrOfBand());

  statfactory::StatFactory stat;
  if(cruleMap[crule_opt[0]]==maxndvi)//ndvi
    assert(ruleBand_opt.size()==2);
  if(cruleMap[crule_opt[0]]==mode){//max voting
    maxBuffer.resize(imgWriter->nrOfCol(),256);//use only byte images for max voting
    for(int iclass=0;iclass<class_opt.size();++iclass)
      assert(class_opt[iclass]<maxBuffer.size());
  }
  double readRow=0;
  double readCol=0;
  double lowerCol=0;
  double upperCol=0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(unsigned int irow=0;irow<imgWriter->nrOfRow();++irow){
    vector<float> lineMask;
    Vector2d< vector<double> > storeBuffer;
    vector<bool> writeValid(ncol);

    //convert irow to geo
    double x=0;
    double y=0;
    imgWriter->image2geo(0,irow,x,y);


    if(cruleMap[crule_opt[0]]==mean ||
       cruleMap[crule_opt[0]]==median ||
       cruleMap[crule_opt[0]]==sum ||
       cruleMap[crule_opt[0]]==minallbands ||
       cruleMap[crule_opt[0]]==maxallbands ||
       cruleMap[crule_opt[0]]==stdev)
      storeBuffer.resize(nband,ncol);
    for(unsigned int icol=0;icol<imgWriter->nrOfCol();++icol){
      writeValid[icol]=false;
      fileBuffer[icol]=0;
      if(cruleMap[crule_opt[0]]==mode){//max voting
        for(int iclass=0;iclass<256;++iclass)
          maxBuffer[icol][iclass]=0;
      }
      else{
        for(unsigned int iband=0;iband<nband;++iband)
          writeBuffer[iband][icol]=dstnodata_opt[0];
      }
    }

    double oldRowMask=-1;//keep track of row mask to optimize number of line readings

    for(unsigned int ifile=0;ifile<size();++ifile){
      //(*imit) already open...
      // try{
      //   (*imit).open((*imit));
      // }
      // catch(string error){
      //   cout << error << endl;
      // }
      // assert((*imit).getDataType()==theType);
      assert((this->at(ifile))->nrOfBand()>=nband);
      if(!(this->at(ifile))->covers(minULX,maxULY,maxLRX,minLRY)){
        // (this->at(ifile)).close();
        continue;
      }
      double uli,ulj,lri,lrj;
      (this->at(ifile))->geo2image(minULX+(magic_x-1.0)*(this->at(ifile))->getDeltaX(),maxULY-(magic_y-1.0)*(this->at(ifile))->getDeltaY(),uli,ulj);
      (this->at(ifile))->geo2image(maxLRX+(magic_x-2.0)*(this->at(ifile))->getDeltaX(),minLRY-(magic_y-2.0)*(this->at(ifile))->getDeltaY(),lri,lrj);
      uli=floor(uli);
      ulj=floor(ulj);
      lri=floor(lri);
      lrj=floor(lrj);
        
      double startCol=uli;
      double endCol=lri;
      if(uli<0)
        startCol=0;
      else if(uli>=(this->at(ifile))->nrOfCol())
        startCol=(this->at(ifile))->nrOfCol()-1;
      if(lri<0)
        endCol=0;
      else if(lri>=(this->at(ifile))->nrOfCol())
        endCol=(this->at(ifile))->nrOfCol()-1;

      //lookup corresponding row for irow in this file
      (this->at(ifile))->geo2image(x,y,readCol,readRow);
      if(readRow<0||readRow>=(this->at(ifile))->nrOfRow()){
        // (this->at(ifile)).close();
        continue;
      }
      // for(int iband=0;iband<(this->at(ifile)).nrOfBand();++iband){
      for(unsigned int iband=0;iband<nband;++iband){
	unsigned int readBand=(band_opt.size()>iband)? band_opt[iband] : iband;
        // readBuffer[iband].resize(readncol);
	try{

          (this->at(ifile))->readData(readBuffer[ifile][iband],startCol,endCol,readRow,readBand,theResample);
	  // if(readRow==0&&iband==0){
	  //   for(unsigned int icol=0;icol<10;++icol)
	  //     cout << readBuffer[0][0][icol] << " ";
	  //   cout << endl;
	  // }
	}
	catch(string error){
	  cerr << "error reading image " << ifile << ": " << endl;
	  throw;
	}
      }
      for(int ib=0;ib<ncol;++ib){
        imgWriter->image2geo(ib,irow,x,y);
	//check mask first
	bool valid=true;
	if(maskReader.isInit()){
	  //read mask
	  double colMask=0;
	  double rowMask=0;

	  maskReader.geo2image(x,y,colMask,rowMask);
	  colMask=static_cast<unsigned int>(colMask);
	  rowMask=static_cast<unsigned int>(rowMask);
	  if(rowMask>=0&&rowMask<maskReader.nrOfRow()&&colMask>=0&&colMask<maskReader.nrOfCol()){
	    if(static_cast<unsigned int>(rowMask)!=static_cast<unsigned int>(oldRowMask)){

	      try{
		maskReader.readData(lineMask,static_cast<unsigned int>(rowMask),mskband_opt[0]);
	      }
	      catch(string errorstring){
		cerr << errorstring << endl;
		exit(1);
	      }
	      catch(...){
		cerr << "error caught" << std::endl;
		exit(3);
	      }
	      oldRowMask=rowMask;
	    }
            for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
              if(lineMask[colMask]==msknodata_opt[ivalue]){
                valid=false;
              }
            }
	  }
	}
	if(!valid)
	  continue;

        //lookup corresponding row for irow in this file
        (this->at(ifile))->geo2image(x,y,readCol,readRow);
        if(readCol<0||readCol>=(this->at(ifile))->nrOfCol())
          continue;
        double val_current=0;
        double val_new=0;
        bool readValid=true;
        switch(theResample){
        case(BILINEAR):
          lowerCol=readCol-0.5;
          lowerCol=static_cast<unsigned int>(lowerCol);
          upperCol=readCol+0.5;
          upperCol=static_cast<unsigned int>(upperCol);
          if(lowerCol<0)
            lowerCol=0;
          if(upperCol>=(this->at(ifile))->nrOfCol())
            upperCol=(this->at(ifile))->nrOfCol()-1;
          for(int vband=0;vband<bndnodata_opt.size();++vband){
            val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][bndnodata_opt[vband]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][bndnodata_opt[vband]][lowerCol-startCol];
	    if(minValue_opt.size()>vband){
	      if(val_new<=minValue_opt[vband]){
		readValid=false;
		break;
	      }
	    }
	    if(maxValue_opt.size()>vband){
	      if(val_new>=maxValue_opt[vband]){
		readValid=false;
		break;
	      }
	    }
	    if(srcnodata_opt.size()>vband){
	      if(val_new==srcnodata_opt[vband]){
		readValid=false;
		break;
	      }
	    }
	  }
          break;
        default:
          readCol=static_cast<unsigned int>(readCol);
          for(int vband=0;vband<bndnodata_opt.size();++vband){
            val_new=readBuffer[ifile][bndnodata_opt[vband]][readCol-startCol];
	    if(minValue_opt.size()>vband){
	      if(val_new<=minValue_opt[vband]){
		readValid=false;
		break;
	      }
	    }
	    if(maxValue_opt.size()>vband){
	      if(val_new>=maxValue_opt[vband]){
		readValid=false;
		break;
	      }
	    }
	    if(srcnodata_opt.size()>vband){
	      if(val_new==srcnodata_opt[vband]){
		readValid=false;
		break;
	      }
	    }
	  }
          break;
	}
	if(readValid){
	  if(file_opt[0]==1)
	    ++fileBuffer[ib];
          if(writeValid[ib]){
            unsigned int iband=0;
	    switch(cruleMap[crule_opt[0]]){
	    case(maxndvi):{//max ndvi
              double red_current=writeBuffer[ruleBand_opt[0]][ib];
              double nir_current=writeBuffer[ruleBand_opt[1]][ib];
	      double ndvi_current=0;
              if(red_current+nir_current>0&&red_current>=0&&nir_current>=0)
                ndvi_current=(nir_current-red_current)/(nir_current+red_current);
	      double ndvi_new=0;
              double red_new=0;
              double nir_new=0;
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                red_new=(readCol-0.5-lowerCol)*readBuffer[ifile][ruleBand_opt[0]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][ruleBand_opt[0]][lowerCol-startCol];
                nir_new=(readCol-0.5-lowerCol)*readBuffer[ifile][ruleBand_opt[1]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][ruleBand_opt[1]][lowerCol-startCol];
                if(red_new+nir_new>0&&red_new>=0&&nir_new>=0)
                  ndvi_new=(nir_new-red_new)/(nir_new+red_new);
                if(ndvi_new>=ndvi_current){
                  for(iband=0;iband<nband;++iband){
                    val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
                    writeBuffer[iband][ib]=val_new;
                  }
		  if(file_opt[0]>1)
		    fileBuffer[ib]=ifile;
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                red_new=readBuffer[ifile][ruleBand_opt[0]][readCol-startCol];
                nir_new=readBuffer[ifile][ruleBand_opt[1]][readCol-startCol];
                if(red_new+nir_new>0&&red_new>=0&&nir_new>=0)
                  ndvi_new=(nir_new-red_new)/(nir_new+red_new);
                if(ndvi_new>=ndvi_current){
                  for(iband=0;iband<nband;++iband){
                    val_new=readBuffer[ifile][iband][readCol-startCol];
                    writeBuffer[iband][ib]=val_new;
                  }
		  if(file_opt[0]>1)
		    fileBuffer[ib]=ifile;
                }
                break;
              }
	      break;
            }
	    case(maxband):
            case(minband):
            case(validband)://max,min,valid band
              val_current=writeBuffer[ruleBand_opt[0]][ib];
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][ruleBand_opt[0]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][ruleBand_opt[0]][lowerCol-startCol];
                val_new*=weight_opt[ifile];
                if((cruleMap[crule_opt[0]]==maxband&&val_new>val_current)||(cruleMap[crule_opt[0]]==minband&&val_new<val_current)||(cruleMap[crule_opt[0]]==validband)){//&&val_new>minValue_opt[0]&&val_new<maxValue_opt[0])){
                  for(iband=0;iband<nband;++iband){
                    val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
                    val_new*=weight_opt[ifile];
                    writeBuffer[iband][ib]=val_new;
                  }
		  if(file_opt[0]>1)
		    fileBuffer[ib]=ifile;
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                val_new=readBuffer[ifile][ruleBand_opt[0]][readCol-startCol];
                val_new*=weight_opt[ifile];
                if((cruleMap[crule_opt[0]]==maxband&&val_new>val_current)||(cruleMap[crule_opt[0]]==minband&&val_new<val_current)||(cruleMap[crule_opt[0]]==validband)){//&&val_new>minValue_opt[0]&&val_new<maxValue_opt[0])){
                  for(iband=0;iband<nband;++iband){
                    val_new=readBuffer[ifile][iband][readCol-startCol];
                    val_new*=weight_opt[ifile];
                    writeBuffer[iband][ib]=val_new;
                  }
		  if(file_opt[0]>1)
		    fileBuffer[ib]=ifile;
                }
                break;
              }
	      break;
            case(mode)://max voting (only for Byte images)
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
                  // ++(maxBuffer[ib][val_new]);
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[ifile][iband][readCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
		}
                break;
	      }
              break;
            case(mean)://mean value
	    case(median)://median value
	    case(sum)://sum value
	    case(minallbands)://minimum for each and every band
	    case(maxallbands)://maximum for each and every band
	    case(stdev)://maximum for each and every band
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[ifile][iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                  assert(ifile>0);
                  // assert(weight_opt[ifile]>=0);
                  // assert(storeBuffer[iband][ib].back()>=0);
                }
                break;
              }
	    if(file_opt[0]>1)
	      fileBuffer[ib]=ifile;
	      break;
	    case(overwrite):
	    default:
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[ifile][iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              }
	    if(file_opt[0]>1)
	      fileBuffer[ib]=ifile;
            break;
	    }
	  }
	  else{
            writeValid[ib]=true;//readValid was true
            unsigned int iband=0;
	    switch(cruleMap[crule_opt[0]]){
            case(mean):
            case(median):
            case(sum):
            case(minallbands):
            case(maxallbands):
            case(stdev):
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[ifile][iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                }
                break;
              }
	    if(file_opt[0]>1)
	      fileBuffer[ib]=ifile;
	    break;
            case(mode):
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
                  // ++(maxBuffer[ib][val_new]);
		}
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                for(iband=0;iband<nband;++iband){
		  val_new=readBuffer[ifile][iband][readCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
		}
                  // ++(maxBuffer[ib][val_new]);
                break;
              }
              break;
            default:
              switch(theResample){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<unsigned int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<unsigned int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=(this->at(ifile))->nrOfCol())
                  upperCol=(this->at(ifile))->nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[ifile][iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ifile][iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              default:
                readCol=static_cast<unsigned int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[ifile][iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              }
	      if(file_opt[0]>1)
		fileBuffer[ib]=ifile;
              break;
            }
          }
        }
      }
      // (this->at(ifile)).close();
    }
    if(cruleMap[crule_opt[0]]==mode){
      vector<short> classBuffer(imgWriter->nrOfCol());
      if(class_opt.size()>1){
        for(int iclass=0;iclass<class_opt.size();++iclass){
          for(unsigned int icol=0;icol<imgWriter->nrOfCol();++icol)
            classBuffer[icol]=maxBuffer[icol][class_opt[iclass]];
          try{
            imgWriter->writeData(classBuffer,irow,iclass);
          }
          catch(string error){
            cerr << "error writing image file: " << error << endl;
            throw;
          }
        }
      }
      else{
        for(unsigned int icol=0;icol<imgWriter->nrOfCol();++icol){
          vector<short>::iterator maxit=maxBuffer[icol].begin();
          maxit=stat.mymax(maxBuffer[icol],maxBuffer[icol].begin(),maxBuffer[icol].end());
          writeBuffer[0][icol]=distance(maxBuffer[icol].begin(),maxit);
	  if(file_opt[0]>1)
	    fileBuffer[icol]=*(maxit);
        }
        try{
          imgWriter->writeData(writeBuffer[0],irow,0);
          if(file_opt[0])
            imgWriter->writeData(fileBuffer,irow,1);
        }
        catch(string error){
          cerr << "error writing image file: " << error << endl;
          throw;
        }
      }
    }
    else{
      for(unsigned int iband=0;iband<bands.size();++iband){
        // assert(writeBuffer[bands[iband]].size()==imgWriter->nrOfCol());
        assert(writeBuffer[iband].size()==imgWriter->nrOfCol());
        for(unsigned int icol=0;icol<imgWriter->nrOfCol();++icol){
	  try{
	    switch(cruleMap[crule_opt[0]]){
	    case(mean):
	      // writeBuffer[iband][icol]=stat.mean(storeBuffer[bands[iband]][icol]);
	      writeBuffer[iband][icol]=stat.mean(storeBuffer[iband][icol]);
	      break;
	    case(median):
	      // writeBuffer[iband][icol]=stat.median(storeBuffer[bands[iband]][icol]);
	      writeBuffer[iband][icol]=stat.median(storeBuffer[iband][icol]);
	      break;
	    case(sum):
	      // writeBuffer[iband][icol]=stat.sum(storeBuffer[bands[iband]][icol]);
	      writeBuffer[iband][icol]=stat.sum(storeBuffer[iband][icol]);
	      break;
	    case(minallbands):
	      // writeBuffer[iband][icol]=stat.mymin(storeBuffer[bands[iband]][icol]);
	      writeBuffer[iband][icol]=stat.mymin(storeBuffer[iband][icol]);
	      break;
	    case(maxallbands):
	      // writeBuffer[iband][icol]=stat.mymax(storeBuffer[bands[iband]][icol]);
	      writeBuffer[iband][icol]=stat.mymax(storeBuffer[iband][icol]);
	      break;
	    case(stdev):
	      // writeBuffer[iband][icol]=sqrt(stat.var(storeBuffer[bands[iband]][icol]));
	      writeBuffer[iband][icol]=sqrt(stat.var(storeBuffer[iband][icol]));
	      break;
	    default:
	      break;
	    }
	  }
	  catch(string error){
	    if(verbose_opt[0])
	      cerr << error << endl;
	    writeBuffer[iband][icol]=dstnodata_opt[0];
	    continue;
	  }
        }
        try{
          imgWriter->writeData(writeBuffer[iband],irow,iband);
        }
        catch(string error){
          cerr << error << endl;
          throw;
        }
      }
      if(file_opt[0]){
        try{
          imgWriter->writeData(fileBuffer,irow,bands.size());
        }
        catch(string error){
          cerr << error << endl;
          throw;
        }
      }
    }
    progress=static_cast<float>(irow+1.0)/imgWriter->nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  if(extent_opt.size()&&(cut_opt[0]||eoption_opt.size())){
    extentReader.close();
  }
  if(maskReader.isInit())
    maskReader.close();
  return(0);
}

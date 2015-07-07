/**********************************************************************
pkcrop.cc: perform raster data operations on image such as crop, extract and stack bands
Copyright (C) 2008-2014 Pieter Kempeneers

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
#include <assert.h>
#include <cstdlib>
#include <string>
#include <list>
#include <iostream>
#include <algorithm>
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Optionpk.h"
#include "algorithms/Egcs.h"

/******************************************************************************/
/*! \page pkcrop pkcrop
 perform raster data operations on image such as crop, extract and stack bands
## SYNOPSIS

<code>
  Usage: pkcrop -i input -o output
</code>

<code>

  Options: [-of out_format] [-ot {Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}] [-b band]* [-ulx ULX -uly ULY -lrx LRX -lry LRY] [-dx xres] [-dy yres] [-r resampling_method] [-a_srs epsg:number] [-nodata value] 

  Advanced options:
  	   [-e vector [-cut]] [-sband band -eband band]* [-co NAME=VALUE]* [-x center_x -y center_y] [-nx size_x -ny size_y] [-ns nsample -nl nlines] [-as min -as max] [-scale value]* [-off offset]* [-ct colortable] [-d description] 
</code>

\section pkcrop_description Description

The utility pkcrop can subset and stack raster images. In the spatial domain it can crop a bounding box from a larger image. The output bounding box is selected by setting the new corner coordinates using the options -ulx -uly -lrx -lry. Alternatively you can set the new image center (-x -y) and size. This can be done either in projected coordinates (using the options -nx -ny) or in image coordinates (using the options -ns -nl). You can also use a vector file to set the new bounding box (option -e). In the spectral domain, pkcrop allows you to select individual bands from one or more input image(s). Bands are stored in the same order as provided on the command line, using the option -b. Band numbers start with index 0 (indicating the first band). The default is to select all input bands. If more input images are provided, the bands are stacked into a multi-band image. If the bounding boxes or spatial resolution are not identical for all input images, you should explicitly set them via the options. The pkcrop utility is not suitable to mosaic or composite images. Consider the utility pkcomposite instead.\section pkcrop_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image file(s). If input contains multiple images, a multi-band output is created | 
 | o      | output               | std::string |       |Output image file | 
 | a_srs  | a_srs                | std::string |       |Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid | 
 | ulx    | ulx                  | double | 0     |Upper left x value bounding box | 
 | uly    | uly                  | double | 0     |Upper left y value bounding box | 
 | lrx    | lrx                  | double | 0     |Lower right x value bounding box | 
 | lry    | lry                  | double | 0     |Lower right y value bounding box | 
 | b      | band                 | unsigned short |       |band index to crop (leave empty to retain all bands) | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 | as     | autoscale            | double |       |scale output to min and max, e.g., --autoscale 0 --autoscale 255 | 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string |       |Output image format (see also gdal_translate). Empty string: inherit from input image | 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 
 | dx     | dx                   | double |       |Output resolution in x (in meter) (empty: keep original resolution) | 
 | dy     | dy                   | double |       |Output resolution in y (in meter) (empty: keep original resolution) | 
 | r      | resampling-method    | std::string | near  |Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation). | 
 | e      | extent               | std::string |       |get boundary from extent from polygons in vector file | 
 | cut      | crop_to_cutline    | bool | false |Crop the extent of the target dataset to the extent of the cutline | 
 | m      | mask                 | std::string |       |Use the first band of the specified file as a validity mask (0 is nodata) | 
 | msknodata | msknodata            | float | 0     |Mask value not to consider for crop
 | msknodata | msknodata            | float | 0     |Mask value not to consider for crop
 | mskband | mskband              | short | 0     |Mask band to read (0 indexed). Provide band for each mask. | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | x      | x                    | double |       |x-coordinate of image center to crop (in meter) | 
 | y      | y                    | double |       |y-coordinate of image center to crop (in meter) | 
 | nx     | nx                   | double |       |image size in x to crop (in meter) | 
 | ny     | ny                   | double |       |image size in y to crop (in meter) | 
 | ns     | ns                   | int  |       |number of samples  to crop (in pixels) | 
 | nl     | nl                   | int  |       |number of lines to crop (in pixels) | 
 | scale  | scale                | double |       |output=scale*input+offset | 
 | off    | offset               | double |       |output=scale*input+offset | 
 | nodata | nodata               | float |       |Nodata value to put in image if out of bounds. | 
 | d      | description          | std::string |       |Set image description | 

Usage: pkcrop -i input -o output


Examples
========
Some examples how to use pkcrop can be found \ref examples_pkcrop "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string>  input_opt("i", "input", "Input image file(s). If input contains multiple images, a multi-band output is created");
  Optionpk<string>  output_opt("o", "output", "Output image file");
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  //todo: support layer names
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file");
  Optionpk<bool> cut_opt("cut", "crop_to_cutline", "Crop the extent of the target dataset to the extent of the cutline.",false);
  Optionpk<string> mask_opt("m", "mask", "Use the first band of the specified file as a validity mask (0 is nodata).");
  Optionpk<float> msknodata_opt("msknodata", "msknodata", "Mask value not to consider for crop.", 0);
  Optionpk<short> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed). Provide band for each mask.", 0);
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
  Optionpk<int> ns_opt("ns", "ns", "number of samples  to crop (in pixels)");
  Optionpk<int> nl_opt("nl", "nl", "number of lines to crop (in pixels)");
  Optionpk<unsigned short>  band_opt("b", "band", "band index to crop (leave empty to retain all bands)");
  Optionpk<unsigned short> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned short> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<double> autoscale_opt("as", "autoscale", "scale output to min and max, e.g., --autoscale 0 --autoscale 255");
  Optionpk<double> scale_opt("scale", "scale", "output=scale*input+offset");
  Optionpk<double> offset_opt("offset", "offset", "output=scale*input+offset");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<float>  nodata_opt("nodata", "nodata", "Nodata value to put in image if out of bounds.");
  Optionpk<string>  resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation).", "near");
  Optionpk<string>  description_opt("d", "description", "Set image description");
  Optionpk<short>  verbose_opt("v", "verbose", "verbose", 0,2);

  extent_opt.setHide(1);
  cut_opt.setHide(1);
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
  
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    bstart_opt.retrieveOption(argc,argv);
    bend_opt.retrieveOption(argc,argv);
    autoscale_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
    extent_opt.retrieveOption(argc,argv);
    cut_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    mskband_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    cx_opt.retrieveOption(argc,argv);
    cy_opt.retrieveOption(argc,argv);
    nx_opt.retrieveOption(argc,argv);
    ny_opt.retrieveOption(argc,argv);
    ns_opt.retrieveOption(argc,argv);
    nl_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(verbose_opt[0])
    cout << setprecision(12) << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;

  if(!doProcess){
    cout << endl;
    cout << "Usage: pkcrop -i input -o output" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }
  if(input_opt.empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(0);
  }
  if(output_opt.empty()){
    std::cerr << "No output file provided (use option -o). Use --help for help information" << std::endl;
    exit(0);
  }

  float nodataValue=nodata_opt.size()? nodata_opt[0] : 0;
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
  ImgReaderGdal imgReader;
  ImgWriterGdal imgWriter;
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
	for(int iband=bstart_opt[ipair];iband<=bend_opt[ipair];++iband)
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
  for(int iimg=0;iimg<input_opt.size();++iimg){
    imgReader.open(input_opt[iimg]);
    if(!isGeoRef)
      isGeoRef=imgReader.isGeoRef();
    if(imgReader.isGeoRef()&&projection_opt.empty())
      projectionString=imgReader.getProjection();
    if(dx_opt.empty()){
      if(!iimg||imgReader.getDeltaX()<dx)
        dx=imgReader.getDeltaX();
    }
    if(dy_opt.empty()){
      if(!iimg||imgReader.getDeltaY()<dy)
        dy=imgReader.getDeltaY();
    }
    if(band_opt.size())
      ncropband+=band_opt.size();
    else
      ncropband+=imgReader.nrOfBand();
    imgReader.close();
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
    if(cut_opt.size())
      extentReader.open(extent_opt[0]);
  }
  else if(cx_opt.size()&&cy_opt.size()&&nx_opt.size()&&ny_opt.size()){
    ulx_opt[0]=cx_opt[0]-nx_opt[0]/2.0;
    uly_opt[0]=(isGeoRef) ? cy_opt[0]+ny_opt[0]/2.0 : cy_opt[0]-ny_opt[0]/2.0;
    lrx_opt[0]=cx_opt[0]+nx_opt[0]/2.0;
    lry_opt[0]=(isGeoRef) ? cy_opt[0]-ny_opt[0]/2.0 : cy_opt[0]+ny_opt[0]/2.0;
    // if(cropulx<ulx_opt[0])
    //   cropulx=ulx_opt[0];
    // if(cropuly>uly_opt[0])
    //   cropuly=uly_opt[0];
    // if(croplrx>lrx_opt[0])
    //   croplrx=lrx_opt[0];
    // if(croplry<lry_opt[0])
    //   croplry=lry_opt[0];
  }
  else if(cx_opt.size()&&cy_opt.size()&&ns_opt.size()&&nl_opt.size()){
    ulx_opt[0]=cx_opt[0]-ns_opt[0]*dx/2.0;
    uly_opt[0]=(isGeoRef) ? cy_opt[0]+nl_opt[0]*dy/2.0 : cy_opt[0]-nl_opt[0]*dy/2.0;
    lrx_opt[0]=cx_opt[0]+ns_opt[0]*dx/2.0;
    lry_opt[0]=(isGeoRef) ? cy_opt[0]-nl_opt[0]*dy/2.0 : cy_opt[0]+nl_opt[0]*dy/2.0;
    // if(cropulx<ulx_opt[0])
    //   cropulx=ulx_opt[0];
    // if(cropuly>uly_opt[0])
    //   cropuly=uly_opt[0];
    // if(croplrx>lrx_opt[0])
    //   croplrx=lrx_opt[0];
    // if(croplry<lry_opt[0])
    //   croplry=lry_opt[0];
  }

  if(verbose_opt[0])
    cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;

  int ncropcol=0;
  int ncroprow=0;

  ImgWriterGdal maskWriter;
  if(extent_opt.size()&&cut_opt[0]){
    try{
      ncropcol=abs(static_cast<int>(ceil((lrx_opt[0]-ulx_opt[0])/dx)));
      ncroprow=abs(static_cast<int>(ceil((uly_opt[0]-lry_opt[0])/dy)));
      maskWriter.open("/vsimem/mask.tif",ncropcol,ncroprow,1,GDT_Float32,"GTiff",option_opt);
      double gt[6];
      gt[0]=ulx_opt[0];
      gt[1]=dx;
      gt[2]=0;
      gt[3]=uly_opt[0];
      gt[4]=0;
      gt[5]=-dy;
      maskWriter.setGeoTransform(gt);
      if(projection_opt.size())
	maskWriter.setProjectionProj4(projection_opt[0]);
      else if(projectionString.size())
	maskWriter.setProjection(projectionString);
	
      //todo: handle multiple extent options
      vector<double> burnValues(1,1);//burn value is 1 (single band)
      maskWriter.rasterizeOgr(extentReader,burnValues);
      maskWriter.close();
    }
    catch(string error){
      cerr << error << std::endl;
      exit(2);
    }
    catch(...){
      cerr << "error catched" << std::endl;
      exit(1);
    }
    //todo: support multiple masks
    mask_opt.clear();
    mask_opt.push_back("/vsimem/mask.tif");
  }
  ImgReaderGdal maskReader;
  if(mask_opt.size()){
    try{
      if(verbose_opt[0]>=1)
	std::cout << "opening mask image file " << mask_opt[0] << std::endl;
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
    catch(...){
      cerr << "error catched" << std::endl;
      exit(1);
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
    // while(autoscale_opt.size()<band_opt.size()*2){
    //   autoscale_opt.push_back(autoscale_opt[0]);
    //   autoscale_opt.push_back(autoscale_opt[1]);
    // }
  }

  for(int iimg=0;iimg<input_opt.size();++iimg){
    if(verbose_opt[0])
      cout << "opening image " << input_opt[iimg] << endl;
    imgReader.open(input_opt[iimg]);
    //if output type not set, get type from input image
    if(theType==GDT_Unknown){
      theType=imgReader.getDataType();
      if(verbose_opt[0])
        cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
    }
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=imgReader.getInterleave();
      option_opt.push_back(theInterleave);
    }
    int nrow=imgReader.nrOfRow();
    int ncol=imgReader.nrOfCol();
    // if(!dx||!dy){
    //   dx=imgReader.getDeltaX();
    //   dy=imgReader.getDeltaY();
    // }
    if(verbose_opt[0])
      cout << "size of " << input_opt[iimg] << ": " << ncol << " cols, "<< nrow << " rows" << endl;
    double uli,ulj,lri,lrj;//image coordinates
    bool forceEUgrid=false;
    if(projection_opt.size())
      forceEUgrid=(!(projection_opt[0].compare("EPSG:3035"))||!(projection_opt[0].compare("EPSG:3035"))||projection_opt[0].find("ETRS-LAEA")!=string::npos);
    if(ulx_opt[0]>=lrx_opt[0]){//default bounding box: no cropping
      uli=0;
      lri=imgReader.nrOfCol()-1;
      ulj=0;
      lrj=imgReader.nrOfRow()-1;
      ncropcol=imgReader.nrOfCol();
      ncroprow=imgReader.nrOfRow();
      imgReader.getBoundingBox(cropulx,cropuly,croplrx,croplry);
      double magicX=1,magicY=1;
      // imgReader.getMagicPixel(magicX,magicY);
      if(forceEUgrid){
	//force to LAEA grid
	Egcs egcs;
        egcs.setLevel(egcs.res2level(dx));
	egcs.force2grid(cropulx,cropuly,croplrx,croplry);
	imgReader.geo2image(cropulx+(magicX-1.0)*imgReader.getDeltaX(),cropuly-(magicY-1.0)*imgReader.getDeltaY(),uli,ulj);
	imgReader.geo2image(croplrx+(magicX-2.0)*imgReader.getDeltaX(),croplry-(magicY-2.0)*imgReader.getDeltaY(),lri,lrj);
      }
      imgReader.geo2image(cropulx+(magicX-1.0)*imgReader.getDeltaX(),cropuly-(magicY-1.0)*imgReader.getDeltaY(),uli,ulj);
      imgReader.geo2image(croplrx+(magicX-2.0)*imgReader.getDeltaX(),croplry-(magicY-2.0)*imgReader.getDeltaY(),lri,lrj);
      //test
      ncropcol=abs(static_cast<int>(ceil((croplrx-cropulx)/dx)));
      ncroprow=abs(static_cast<int>(ceil((cropuly-croplry)/dy)));
    }
    else{
      double magicX=1,magicY=1;
      // imgReader.getMagicPixel(magicX,magicY);
      cropulx=ulx_opt[0];
      cropuly=uly_opt[0];
      croplrx=lrx_opt[0];
      croplry=lry_opt[0];
      if(forceEUgrid){
	//force to LAEA grid
	Egcs egcs;
        egcs.setLevel(egcs.res2level(dx));
	egcs.force2grid(cropulx,cropuly,croplrx,croplry);
	imgReader.geo2image(cropulx+(magicX-1.0)*imgReader.getDeltaX(),cropuly-(magicY-1.0)*imgReader.getDeltaY(),uli,ulj);
	imgReader.geo2image(croplrx+(magicX-2.0)*imgReader.getDeltaX(),croplry-(magicY-2.0)*imgReader.getDeltaY(),lri,lrj);
      }
      imgReader.geo2image(cropulx+(magicX-1.0)*imgReader.getDeltaX(),cropuly-(magicY-1.0)*imgReader.getDeltaY(),uli,ulj);
      imgReader.geo2image(croplrx+(magicX-2.0)*imgReader.getDeltaX(),croplry-(magicY-2.0)*imgReader.getDeltaY(),lri,lrj);

      ncropcol=abs(static_cast<int>(ceil((croplrx-cropulx)/dx)));
      ncroprow=abs(static_cast<int>(ceil((cropuly-croplry)/dy)));
      uli=floor(uli);
      ulj=floor(ulj);
      lri=floor(lri);
      lrj=floor(lrj);
    }

    double dcropcol=0;
    double dcroprow=0;
    double deltaX=imgReader.getDeltaX();
    double deltaY=imgReader.getDeltaY();
    dcropcol=(lri-uli+1)/(dx/deltaX);
    dcroprow=(lrj-ulj+1)/(dy/deltaY);
    if(!imgWriter.nrOfBand()){//not opened yet
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
      // string theCompression;
      // if(compress_opt[0]!="")//default
      //   theCompression=compress_opt[0];
      // else
      //   theCompression=imgReader.getCompression();
      // string theInterleave;
      // if(interleave_opt[0]!="")//default
      //   theInterleave=interleave_opt[0];
      // else
      //   theInterleave=imgReader.getInterleave();
      string imageType=imgReader.getImageType();
      if(oformat_opt.size())//default
        imageType=oformat_opt[0];
      try{
        imgWriter.open(output_opt[0],ncropcol,ncroprow,ncropband,theType,imageType,option_opt);
	if(nodata_opt.size()){
	  for(int iband=0;iband<ncropband;++iband)
	    imgWriter.GDALSetNoDataValue(nodata_opt[0],iband);
	}
      }
      catch(string errorstring){
        cout << errorstring << endl;
        exit(4);
      }
      if(description_opt.size())
	imgWriter.setImageDescription(description_opt[0]);
      double gt[6];
      gt[0]=cropulx;
      gt[1]=dx;
      gt[2]=0;
      gt[3]=cropuly;
      gt[4]=0;
      gt[5]=(imgReader.isGeoRef())? -dy : dy;
      imgWriter.setGeoTransform(gt);
      if(projection_opt.size()){
	if(verbose_opt[0])
	  cout << "projection: " << projection_opt[0] << endl;
	imgWriter.setProjectionProj4(projection_opt[0]);
      }
      else
	imgWriter.setProjection(imgReader.getProjection());
      if(imgWriter.getDataType()==GDT_Byte){
	if(colorTable_opt.size()){
	  if(colorTable_opt[0]!="none")
	    imgWriter.setColorTable(colorTable_opt[0]);
	}
	else if (imgReader.getColorTable()!=NULL)//copy colorTable from input image
	  imgWriter.setColorTable(imgReader.getColorTable());
      }
    }

    double startCol=uli;
    double endCol=lri;
    if(uli<0)
      startCol=0;
    else if(uli>=imgReader.nrOfCol())
      startCol=imgReader.nrOfCol()-1;
    if(lri<0)
      endCol=0;
    else if(lri>=imgReader.nrOfCol())
      endCol=imgReader.nrOfCol()-1;
    double startRow=ulj;
    double endRow=lrj;
    if(ulj<0)
      startRow=0;
    else if(ulj>=imgReader.nrOfRow())
      startRow=imgReader.nrOfRow()-1;
    if(lrj<0)
      endRow=0;
    else if(lrj>=imgReader.nrOfRow())
      endRow=imgReader.nrOfRow()-1;



    int readncol=endCol-startCol+1;
    vector<double> readBuffer(readncol+1);
    int nband=(band_opt.size())?band_opt.size() : imgReader.nrOfBand();
    for(int iband=0;iband<nband;++iband){
      int readBand=(band_opt.size()>iband)?band_opt[iband]:iband;
      if(verbose_opt[0]){
	cout << "extracting band " << readBand << endl;
	pfnProgress(progress,pszMessage,pProgressArg);
      }
      double theMin=0;
      double theMax=0;
      if(autoscale_opt.size()){
	try{
	  imgReader.getMinMax(static_cast<int>(startCol),static_cast<int>(endCol),static_cast<int>(startRow),static_cast<int>(endRow),readBand,theMin,theMax);
	}
	catch(string errorString){
	  cout << errorString << endl;
	}
	if(verbose_opt[0])
	  cout << "minmax: " << theMin << ", " << theMax << endl;
	double theScale=(autoscale_opt[1]-autoscale_opt[0])/(theMax-theMin);
	double theOffset=autoscale_opt[0]-theScale*theMin;
	imgReader.setScale(theScale,readBand);
	imgReader.setOffset(theOffset,readBand);
      }	
      else{
	if(scale_opt.size()){
	  if(scale_opt.size()>iband)
	    imgReader.setScale(scale_opt[iband],readBand);
	  else
	    imgReader.setScale(scale_opt[0],readBand);
	}
	if(offset_opt.size()){
	  if(offset_opt.size()>iband)
	    imgReader.setOffset(offset_opt[iband],readBand);
	  else
	    imgReader.setOffset(offset_opt[0],readBand);
	}
      }

      double readRow=0;
      double readCol=0;
      double lowerCol=0;
      double upperCol=0;
      for(int irow=0;irow<imgWriter.nrOfRow();++irow){
	vector<float> lineMask;
	double x=0;
	double y=0;
	//convert irow to geo
	imgWriter.image2geo(0,irow,x,y);
	//lookup corresponding row for irow in this file
	imgReader.geo2image(x,y,readCol,readRow);
	vector<double> writeBuffer;
	if(readRow<0||readRow>=imgReader.nrOfRow()){
	  //if(readRow<0)
	  //readRow=0;
	  //else if(readRow>=imgReader.nrOfRow())
	  //readRow=imgReader.nrOfRow()-1;
	  for(int icol=0;icol<imgWriter.nrOfCol();++icol)
	    writeBuffer.push_back(nodataValue);
	}
	else{
	  if(verbose_opt[0]>1)
	    cout << "reading row: " << readRow << endl;
	  try{
            if(endCol<imgReader.nrOfCol()-1)
              imgReader.readData(readBuffer,GDT_Float64,startCol,endCol+1,readRow,readBand,theResample);
            else
              imgReader.readData(readBuffer,GDT_Float64,startCol,endCol,readRow,readBand,theResample);
	    // for(int icol=0;icol<ncropcol;++icol){
	    double oldRowMask=-1;//keep track of row mask to optimize number of line readings
	    for(int icol=0;icol<imgWriter.nrOfCol();++icol){
	      imgWriter.image2geo(icol,irow,x,y);
	      //lookup corresponding row for irow in this file
	      imgReader.geo2image(x,y,readCol,readRow);
	      if(readCol<0||readCol>=imgReader.nrOfCol()){
	      // if(readCol<0||readCol>=imgReader.nrOfCol()){
		//               if(readCol<0)
		//                 readCol=0;
		//               else if(readCol>=imgReader.nrOfCol())
		//                 readCol=imgReader.nrOfCol()-1;
		writeBuffer.push_back(nodataValue);
	      }
	      else{
                bool valid=true;
		double geox=0;
		double geoy=0;
                if(mask_opt.size()){
		  //read mask
		  double colMask=0;
		  double rowMask=0;

		  imgWriter.image2geo(icol,irow,geox,geoy);
		  maskReader.geo2image(geox,geoy,colMask,rowMask);
		  colMask=static_cast<int>(colMask);
		  rowMask=static_cast<int>(rowMask);
		  if(rowMask>=0&&rowMask<maskReader.nrOfRow()&&colMask>=0&&colMask<maskReader.nrOfCol()){
		    if(static_cast<int>(rowMask)!=static_cast<int>(oldRowMask)){

		      assert(rowMask>=0&&rowMask<maskReader.nrOfRow());
		      try{
			maskReader.readData(lineMask,GDT_Float32,static_cast<int>(rowMask),mskband_opt[0]);
		      }
		      catch(string errorstring){
			cerr << errorstring << endl;
			exit(1);
		      }
		      catch(...){
			cerr << "error catched" << std::endl;
			exit(3);
		      }
		      oldRowMask=rowMask;
		    }
		    if(lineMask[colMask]==msknodata_opt[0])
		      valid=false;
		  }
		}

                if(!valid)
                  writeBuffer.push_back(nodataValue);
                else{
                  switch(theResample){
                  case(BILINEAR):
                    lowerCol=readCol-0.5;
                    lowerCol=static_cast<int>(lowerCol);
                    upperCol=readCol+0.5;
                    upperCol=static_cast<int>(upperCol);
                    if(lowerCol<0)
                      lowerCol=0;
                    if(upperCol>=imgReader.nrOfCol())
                      upperCol=imgReader.nrOfCol()-1;
                    // writeBuffer.push_back((readCol-0.5-lowerCol)*(readBuffer[upperCol-startCol]*theScale+theOffset)+(1-readCol+0.5+lowerCol)*(readBuffer[lowerCol-startCol]*theScale+theOffset));
                    writeBuffer.push_back((readCol-0.5-lowerCol)*readBuffer[upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[lowerCol-startCol]);
                    break;
                  default:
                    readCol=static_cast<int>(readCol);
                    readCol-=startCol;//we only start reading from startCol
                    // writeBuffer.push_back(readBuffer[readCol]*theScale+theOffset);
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
	if(writeBuffer.size()!=imgWriter.nrOfCol())
	  cout << "writeBuffer.size()=" << writeBuffer.size() << ", imgWriter.nrOfCol()=" << imgWriter.nrOfCol() << endl;
	assert(writeBuffer.size()==imgWriter.nrOfCol());
	try{
	  imgWriter.writeData(writeBuffer,GDT_Float64,irow,writeBand);
	}
	catch(string errorstring){
	  cout << errorstring << endl;
	  exit(3);
	}
	if(verbose_opt[0]){
	  progress=(1.0+irow);
	  progress/=imgWriter.nrOfRow();
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
	else{
	  progress=(1.0+irow);
	  progress+=(imgWriter.nrOfRow()*writeBand);
	  progress/=imgWriter.nrOfBand()*imgWriter.nrOfRow();
	  assert(progress>=0);
	  assert(progress<=1);
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
      }
      ++writeBand;
    }
    imgReader.close();
  }
  if(extent_opt.size()&&cut_opt.size()){
    extentReader.close();
  }
  if(mask_opt.size())
    maskReader.close();
  imgWriter.close();
}

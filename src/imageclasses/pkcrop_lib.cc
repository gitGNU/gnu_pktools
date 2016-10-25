/**********************************************************************
pkcrop_lib.cc: perform raster data operations on image such as crop, extract and stack bands
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
#include <assert.h>
#include <string>
#include <iostream>
#include <algorithm>
#include <memory>
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgCollection.h"
#include "base/Optionpk.h"
#include "algorithms/Egcs.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRasterGdal> ImgCollection::crop(const AppFactory& app){
  shared_ptr<ImgRasterGdal> imgWriter=ImgRasterGdal::createImg();
  crop(*imgWriter, app);
  return(imgWriter);
}

/**
 * @param imgWriter output raster crop dataset
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgCollection::crop(ImgRasterGdal& imgWriter, const AppFactory& app){
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  //todo: support layer names
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file");
  Optionpk<bool> cut_opt("cut", "crop_to_cutline", "Crop the extent of the target dataset to the extent of the cutline.",false);
  Optionpk<string> eoption_opt("eo","eo", "special extent options controlling rasterization: ATTRIBUTE|CHUNKYSIZE|ALL_TOUCHED|BURN_VALUE_FROM|MERGE_ALG, e.g., -eo ATTRIBUTE=fieldname");
  Optionpk<string> mask_opt("m", "mask", "Use the the specified file as a validity mask (0 is nodata).");
  Optionpk<unsigned int> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed)", 0);
  Optionpk<double> msknodata_opt("msknodata", "msknodata", "Mask value not to consider for crop.", 0);
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
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  // Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<double>  nodata_opt("nodata", "nodata", "Nodata value to put in image if out of bounds.");
  Optionpk<string>  resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbor, bilinear: bi-linear interpolation).", "near");
  Optionpk<string>  description_opt("d", "description", "Set image description");
  Optionpk<bool>  align_opt("align", "align", "Align output bounding box to input image",false);
  Optionpk<short>  verbose_opt("v", "verbose", "verbose", 0,2);

  extent_opt.setHide(1);
  cut_opt.setHide(1);
  eoption_opt.setHide(1);
  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  mask_opt.setHide(1);
  msknodata_opt.setHide(1);
  mskband_opt.setHide(1);
  // option_opt.setHide(1);
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
    mskband_opt.retrieveOption(app.getArgc(),app.getArgv());
    msknodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    // option_opt.retrieveOption(app.getArgc(),app.getArgv());
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
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
  if(!doProcess){
    cout << endl;
    std::ostringstream helpStream;
    helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    //test
    helpStream << "throwing exception" << std::endl;
    throw(helpStream.str());//help was invoked, stop processing
  }
  // if(input_opt.empty()){
  if(empty()){
    std::cerr << "Input collection is empty. Use --help for more help information" << std::endl;
    return(CE_Failure);
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
    return(CE_Failure);
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
    return(CE_Failure);
  }

  bool isGeoRef=false;
  string projectionString;
  // for(int iimg=0;iimg<input_opt.size();++iimg){

  std::vector<std::shared_ptr<ImgRasterGdal> >::const_iterator imit=begin();

  for(imit=begin();imit!=end();++imit){
    // while((imgReader=getNextImage())){
    if(verbose_opt[0])
      cout << "m_index: " << m_index << endl;
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
  if(otype_opt.size()){
    theType=getGDALDataType(otype_opt[0]);
    if(theType==GDT_Unknown)
      std::cout << "Warning: unknown output pixel type: " << otype_opt[0] << ", using input type as default" << std::endl;
  }
  if(verbose_opt[0])
    cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

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
        return(CE_Failure);
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

  ImgRasterGdal maskReader;
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
      return(CE_Failure);
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
      return(CE_Failure);
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
    if(theType==GDT_Unknown){
      theType=(*imit)->getDataType();
      if(verbose_opt[0])
        cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
    }
    // if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    //   string theInterleave="INTERLEAVE=";
    //   theInterleave+=(*imit)->getInterleave();
    //   option_opt.push_back(theInterleave);
    // }
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
      string imageType;//=(*imit)->getImageType();
      if(oformat_opt.size())//default
        imageType=oformat_opt[0];
      try{
        imgWriter.open(ncropcol,ncroprow,ncropband,theType);
        if(nodata_opt.size()){
          imgWriter.setNoData(nodata_opt);
          for(int iband=0;iband<ncropband;++iband)
            imgWriter.GDALSetNoDataValue(nodata_opt[0],iband);
        }
      }
      catch(string errorstring){
        cout << errorstring << endl;
        return(CE_Failure);
      }
      if(description_opt.size())
        imgWriter.setImageDescription(description_opt[0]);
      double gt[6];
      gt[0]=cropulx;
      gt[1]=dx;
      gt[2]=0;
      gt[3]=cropuly;
      gt[4]=0;
      gt[5]=((*imit)->isGeoRef())? -dy : dy;
      imgWriter.setGeoTransform(gt);
      if(projection_opt.size()){
        if(verbose_opt[0])
          cout << "projection: " << projection_opt[0] << endl;
        imgWriter.setProjectionProj4(projection_opt[0]);
      }
      else
        imgWriter.setProjection((*imit)->getProjection());
      if(imgWriter.getDataType()==GDT_Byte){
        if(colorTable_opt.size()){
          if(colorTable_opt[0]!="none")
            imgWriter.setColorTable(colorTable_opt[0]);
        }
        else if ((*imit)->getColorTable()!=NULL)//copy colorTable from input image
          imgWriter.setColorTable((*imit)->getColorTable());
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
      for(int irow=0;irow<imgWriter.nrOfRow();++irow){
        vector<double> lineMask;
        double x=0;
        double y=0;
        //convert irow to geo
        imgWriter.image2geo(0,irow,x,y);
        //lookup corresponding row for irow in this file
        (*imit)->geo2image(x,y,readCol,readRow);
        vector<double> writeBuffer;
        if(readRow<0||readRow>=(*imit)->nrOfRow()){
          for(int icol=0;icol<imgWriter.nrOfCol();++icol)
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
            for(int icol=0;icol<imgWriter.nrOfCol();++icol){
              imgWriter.image2geo(icol,irow,x,y);
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

                  imgWriter.image2geo(icol,irow,geox,geoy);
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
                        return(CE_Failure);
                      }
                      catch(...){
                        cerr << "error caught" << std::endl;
                        return(CE_Failure);
                      }
                      oldRowMask=rowMask;
                    }
                    for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
                      if(lineMask[colMask]==msknodata_opt[ivalue]){
                        if(nodata_opt.size()>ivalue)
                          nodataValue=nodata_opt[ivalue];
                        valid=false;
                        break;
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
            return(CE_Failure);
          }
        }
        if(writeBuffer.size()!=imgWriter.nrOfCol())
          cout << "writeBuffer.size()=" << writeBuffer.size() << ", imgWriter.nrOfCol()=" << imgWriter.nrOfCol() << endl;

        assert(writeBuffer.size()==imgWriter.nrOfCol());
        try{
          imgWriter.writeData(writeBuffer,irow,writeBand);
        }
        catch(string errorstring){
          cout << errorstring << endl;
          return(CE_Failure);
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
    // (*imit)->close();
  }
  if(extent_opt.size()&&(cut_opt[0]||eoption_opt.size())){
    extentReader.close();
  }
  if(maskReader.isInit())
    maskReader.close();
  return(CE_None);
}

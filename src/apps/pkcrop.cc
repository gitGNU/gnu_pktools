/**********************************************************************
pkcrop.cc: perform raster data operations on image such as crop, extract and stack bands
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
#include <assert.h>
#include <string>
#include <list>
#include <iostream>
#include <algorithm>
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Optionpk.h"
#include "algorithms/Egcs.h"

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string>  input_opt("i", "input", "Input image file(s). If input contains multiple images, a multi-band output is created");
  Optionpk<string>  output_opt("o", "output", "Output image file");
  Optionpk<string>  projection_opt("a_srs", "a_srs", "Override the projection for the output file (leave blank to copy from input file, use epsg:3035 to use European projection and force to European grid");
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file");
  Optionpk<bool> mask_opt("m","mask","mask values out of polygon in extent file to flag option (tip: for better performance, use gdal_rasterize -i -burn 0 -l extent extent.shp output (with output the result of pkcrop)",false);
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
  Optionpk<int>  band_opt("b", "band", "band index to crop (leave empty to retain all bands)");
  Optionpk<double> autoscale_opt("as", "autoscale", "scale output to min and max, e.g., --autoscale 0 --autoscale 255");
  Optionpk<double> scale_opt("s", "scale", "output=scale*input+offset");
  Optionpk<double> offset_opt("off", "offset", "output=scale*input+offset");
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<double>  nodata_opt("nodata", "nodata", "Nodata value to put in image if out of bounds.");
  Optionpk<string>  resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbour, bilinear: bi-linear interpolation).", "near");
  Optionpk<string>  description_opt("d", "description", "Set image description");
  Optionpk<bool>  verbose_opt("v", "verbose", "verbose", false);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    extent_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    autoscale_opt.retrieveOption(argc,argv);
    scale_opt.retrieveOption(argc,argv);
    offset_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    cx_opt.retrieveOption(argc,argv);
    cy_opt.retrieveOption(argc,argv);
    nx_opt.retrieveOption(argc,argv);
    ny_opt.retrieveOption(argc,argv);
    ns_opt.retrieveOption(argc,argv);
    nl_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
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

  double nodataValue=nodata_opt.size()? nodata_opt[0] : 0;
  RESAMPLE theResample;
  if(resample_opt[0]=="near"){
    theResample=NEAR;
    if(verbose_opt[0])
      cout << "resampling: nearest neighbour" << endl;
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
  double dx=(dx_opt.size())? dx_opt[0]:0;
  double dy=(dy_opt.size())? dy_opt[0]:0;
  for(int iimg=0;iimg<input_opt.size();++iimg){
    imgReader.open(input_opt[iimg]);
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
    for(int iextent=0;iextent<extent_opt.size();++iextent){
      extentReader.open(extent_opt[iextent]);
      if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
        cerr << "Error: could not get extent from " << extent_opt[0] << endl;
        exit(1);
      }
      extentReader.close();
    }
    if(mask_opt[0])
      extentReader.open(extent_opt[0]);
  }
  else if(cx_opt.size()&&cy_opt.size()&&nx_opt.size()&&ny_opt.size()){
    ulx_opt[0]=cx_opt[0]-nx_opt[0]/2.0;
    uly_opt[0]=cy_opt[0]+ny_opt[0]/2.0;
    lrx_opt[0]=cx_opt[0]+nx_opt[0]/2.0;
    lry_opt[0]=cy_opt[0]-ny_opt[0]/2.0;
  }
  else if(cx_opt.size()&&cy_opt.size()&&ns_opt.size()&&nl_opt.size()){
    ulx_opt[0]=cx_opt[0]-ns_opt[0]*dx/2.0;
    uly_opt[0]=cy_opt[0]+nl_opt[0]*dy/2.0;
    lrx_opt[0]=cx_opt[0]+ns_opt[0]*dx/2.0;
    lry_opt[0]=cy_opt[0]-nl_opt[0]*dy/2.0;
  }
  if(ulx_opt[0]<cropulx)
    cropulx=ulx_opt[0];
  if(uly_opt[0]>cropuly)
    cropuly=uly_opt[0];
  if(lry_opt[0]<croplry)
    croplry=lry_opt[0];
  if(lrx_opt[0]>croplrx)
    croplrx=lrx_opt[0];
  if(verbose_opt[0])
    cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;
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
    int ncropcol=0;
    int ncroprow=0;
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
      // ncropcol=ceil((croplrx-cropulx)/dx);
      // ncroprow=ceil((cropuly-croplry)/dy);
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

      ncropcol=ceil((croplrx-cropulx)/dx);
      ncroprow=ceil((cropuly-croplry)/dy);
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
      gt[5]=-dy;
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
	imgReader.getMinMax(static_cast<int>(startCol),static_cast<int>(endCol),static_cast<int>(startRow),static_cast<int>(endRow),readBand,theMin,theMax);
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
	double x=0;
	double y=0;
	//convert irow to geo
	imgWriter.image2geo(0,irow,x,y);
	//lookup corresponding row for irow in this file
	  //test
	  // cout << "x: " << x << ", y: " << y << endl;
	imgReader.geo2image(x,y,readCol,readRow);
	//test
	// cout << "readRow: " << readRow << ", readCol: " << readCol << flush << endl;
	// double lowerCol=0;
	// double upperCol=0;
	vector<double> writeBuffer;
	if(readRow<0||readRow>=imgReader.nrOfRow()){
	  //if(readRow<0)
	  //readRow=0;
	  //else if(readRow>=imgReader.nrOfRow())
	  //readRow=imgReader.nrOfRow()-1;
	  for(int ib=0;ib<imgWriter.nrOfCol();++ib)
	    writeBuffer.push_back(nodataValue);
	}
	else{
	  if(verbose_opt[0])
	    cout << "reading row: " << readRow << endl;
	  try{
            if(endCol<imgReader.nrOfCol()-1)
              imgReader.readData(readBuffer,GDT_Float64,startCol,endCol+1,readRow,readBand,theResample);
            else
              imgReader.readData(readBuffer,GDT_Float64,startCol,endCol,readRow,readBand,theResample);
	    // for(int ib=0;ib<ncropcol;++ib){
	    for(int ib=0;ib<imgWriter.nrOfCol();++ib){
	      assert(imgWriter.image2geo(ib,irow,x,y));
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
                if(mask_opt[0]&&extent_opt.size()){
                  valid=false;
                  OGRPoint thePoint;
                  thePoint.setX(x);
                  thePoint.setY(y);
                  OGRLayer  *readLayer;
                  readLayer = extentReader.getDataSource()->GetLayer(0);
                  readLayer->ResetReading();
                  OGRFeature *readFeature;
                  while( (readFeature = readLayer->GetNextFeature()) != NULL ){
                    OGRGeometry *poGeometry;
                    poGeometry = readFeature->GetGeometryRef();
                    assert(poGeometry!=NULL);
                    //check if point is on surface
                    OGRPolygon readPolygon = *((OGRPolygon *) poGeometry);
                    readPolygon.closeRings();
                    if(readPolygon.Contains(&thePoint)){
                      valid=true;
                      break;
                    }
                    else
                      continue;
                  }
                }
                if(!valid)
                  writeBuffer.push_back(nodataValue);
                else{
                  // double theScale=1;
                  // double theOffset=0;
		  // if(autoscale_opt.size()){
		  //   theScale=(autoscale_opt[1]-autoscale_opt[0])/(theMax-theMin);
		  //   theOffset=autoscale_opt[0]-theScale*theMin;
		  // }
		  // else{
		  //   theScale=(scale_opt.size()>1)?scale_opt[iband]:scale_opt[0];
		  //   theOffset=(offset_opt.size()>1)?offset_opt[iband]:offset_opt[0];
		  // }
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
  if(extent_opt.size()&&mask_opt[0]){
    extentReader.close();
  }
  imgWriter.close();
}

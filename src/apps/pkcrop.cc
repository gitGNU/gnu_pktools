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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

int main(int argc, char *argv[])
{
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string>  input_opt("i", "input", "Input image file(s). If input contains multiple images, a multi-band output is created", "");
  Optionpk<string>  output_opt("o", "output", "Output image file", "");
  Optionpk<string>  projection_opt("p", "projection", "projection in EPSG format (leave blank to copy from input file, use EPSG:3035 to use European projection and to force to European grid", "");
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file", "");
  Optionpk<bool> mask_opt("m","mask","mask values out of polygon in extent file to flag option (tip: for better performance, use gdal_rasterize -i -burn 0 -l extent extent.shp output (with output the result of pkcrop)",false);
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  dx_opt("dx", "dx", "Output resolution in x (in meter) (0.0: keep original resolution)", 0.0);
  Optionpk<double>  dy_opt("dy", "dy", "Output resolution in y (in meter) (0.0: keep original resolution)", 0.0);
  Optionpk<int>  band_opt("b", "band", "band index to crop (-1: crop all bands)", -1);
  Optionpk<double> scale_opt("s", "scale", "output=scale*input+offset", 1);
  Optionpk<double> offset_opt("off", "offset", "output=scale*input+offset", 0);
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]");
  Optionpk<short>  flag_opt("f", "flag", "Flag value to put in image if out of bounds.", 0);
  Optionpk<string>  resample_opt("r", "resampling-method", "Resampling method (near: nearest neighbour, bilinear: bi-linear interpolation).", "near");
  Optionpk<string>  description_opt("d", "description", "Set image description", "");
  Optionpk<bool>  verbose_opt("v", "verbose", "verbose", false);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  projection_opt.retrieveOption(argc,argv);
  extent_opt.retrieveOption(argc,argv);
  mask_opt.retrieveOption(argc,argv);
  ulx_opt.retrieveOption(argc,argv);
  uly_opt.retrieveOption(argc,argv);
  lrx_opt.retrieveOption(argc,argv);
  lry_opt.retrieveOption(argc,argv);
  band_opt.retrieveOption(argc,argv);
  scale_opt.retrieveOption(argc,argv);
  offset_opt.retrieveOption(argc,argv);
  otype_opt.retrieveOption(argc,argv);
  oformat_opt.retrieveOption(argc,argv);
  colorTable_opt.retrieveOption(argc,argv);
  dx_opt.retrieveOption(argc,argv);
  dy_opt.retrieveOption(argc,argv);
  option_opt.retrieveOption(argc,argv);
  flag_opt.retrieveOption(argc,argv);
  resample_opt.retrieveOption(argc,argv);
  description_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

  if(version_opt[0]||todo_opt[0]){
    cout << version_opt.getHelp() << endl;
    cout << "todo: " << todo_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }
  if(help_opt[0]){
    cout << "usage: pkcrop -i inputimage -o outputimage [OPTIONS]" << endl;
    exit(0);
  }

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
  double dx=dx_opt[0];
  double dy=dy_opt[0];
  //bounding box of cropped image
  double cropulx=ulx_opt[0];
  double cropuly=uly_opt[0];
  double croplrx=lrx_opt[0];
  double croplry=lry_opt[0];
  //get bounding box from extentReader if defined
  ImgReaderOgr extentReader;
  if(extent_opt[0]!=""){
    for(int iextent=0;iextent<extent_opt.size();++iextent){
      extentReader.open(extent_opt[iextent]);
      if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
        cerr << "Error: could not get extent from " << extent_opt[0] << endl;
        exit(1);
      }
      if(ulx_opt[0]<cropulx)
        cropulx=ulx_opt[0];
      if(uly_opt[0]>cropuly)
        cropuly=uly_opt[0];
      if(lry_opt[0]<croplry)
        croplry=lry_opt[0];
      if(lrx_opt[0]>croplrx)
        croplrx=lrx_opt[0];
      extentReader.close();
    }
    if(mask_opt[0])
      extentReader.open(extent_opt[0]);
  }
  if(verbose_opt[0])
    cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;
  //determine number of output bands
  int ncropband=0;//total number of bands to write
  int writeBand=0;//write band

  while(scale_opt.size()<band_opt.size())
    scale_opt.push_back(scale_opt[0]);
  while(offset_opt.size()<band_opt.size())
    offset_opt.push_back(offset_opt[0]);

  for(int iimg=0;iimg<input_opt.size();++iimg){
    imgReader.open(input_opt[iimg]);
    if(band_opt[0]>=0)
      ncropband+=band_opt.size();
    else
      ncropband+=imgReader.nrOfBand();
    imgReader.close();
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
    if(!dx||!dy){
      dx=imgReader.getDeltaX();
      dy=imgReader.getDeltaY();
    }
    if(verbose_opt[0])
      cout << "size of " << input_opt[iimg] << ": " << ncol << " cols, "<< nrow << " rows" << endl;
    double uli,ulj,lri,lrj;//image coordinates
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
      if(!imgReader.getProjection().compare("ETRS-LAEA")||!projection_opt[0].compare("EPSG:3035")||!projection_opt[0].compare("epsg:3035")){
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
    }
    else{
      double magicX=1,magicY=1;
      // imgReader.getMagicPixel(magicX,magicY);
      cropulx=ulx_opt[0];
      cropuly=uly_opt[0];
      croplrx=lrx_opt[0];
      croplry=lry_opt[0];
      if(!imgReader.getProjection().compare("ETRS-LAEA")||!projection_opt[0].compare("EPSG:3035")||!projection_opt[0].compare("epsg:3035")){
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
      }
      catch(string errorstring){
        cout << errorstring << endl;
        exit(4);
      }
      if(description_opt[0]!="")
	imgWriter.setImageDescription(description_opt[0]);
      imgWriter.setGeoTransform(cropulx,cropuly,dx,dy,0,0);
      if(projection_opt[0]!=""){
	if(verbose_opt[0])
	  cout << "projection: " << projection_opt[0] << endl;
	imgWriter.setProjectionProj4(projection_opt[0]);
      }
      else if(imgReader.isGeoRef())
	imgWriter.setProjection(imgReader.getProjection());
      if(colorTable_opt[0]!=""){
        if(verbose_opt[0])
          cout << "set colortable " << colorTable_opt[0] << endl;
        assert(imgWriter.getDataType()==GDT_Byte);
        imgWriter.setColorTable(colorTable_opt[0]);
      }
      else if(imgReader.getColorTable()!=NULL){
        if(verbose_opt[0])
          cout << "set colortable from input image" << endl;
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
    int readncol=endCol-startCol+1;
    vector<double> readBuffer(readncol+1);
    int nband=(band_opt[0]<0)?imgReader.nrOfBand():band_opt.size();
    for(int iband=0;iband<nband;++iband){
      int readBand=(band_opt[0]<0)?iband:band_opt[iband];
      if(verbose_opt[0]){
	cout << "extracting band " << readBand << endl;
	pfnProgress(progress,pszMessage,pProgressArg);
      }
      double readRow=0;
      double readCol=0;
      double lowerCol=0;
      double upperCol=0;
      for(int irow=0;irow<ncroprow;++irow){
	double x=0;
	double y=0;
	//convert irow to geo
	imgWriter.image2geo(0,irow,x,y);
	//lookup corresponding row for irow in this file
	imgReader.geo2image(x,y,readCol,readRow);
	// double lowerCol=0;
	// double upperCol=0;
	vector<double> writeBuffer;
	if(readRow<0||readRow>=imgReader.nrOfRow()){
	  //if(readRow<0)
	  //readRow=0;
	  //else if(readRow>=imgReader.nrOfRow())
	  //readRow=imgReader.nrOfRow()-1;
	  for(int ib=0;ib<ncropcol;++ib)
	    writeBuffer.push_back(flag_opt[0]);
	}
	else{
	  try{
            if(endCol<imgReader.nrOfCol()-1)
              imgReader.readData(readBuffer,GDT_Float64,startCol,endCol+1,readRow,readBand,theResample);
            else
              imgReader.readData(readBuffer,GDT_Float64,startCol,endCol,readRow,readBand,theResample);
	    for(int ib=0;ib<ncropcol;++ib){
	      assert(imgWriter.image2geo(ib,irow,x,y));
	      //lookup corresponding row for irow in this file
	      imgReader.geo2image(x,y,readCol,readRow);
	      if(readCol<0||readCol>=imgReader.nrOfCol()){
	      // if(readCol<0||readCol>=imgReader.nrOfCol()){
		//               if(readCol<0)
		//                 readCol=0;
		//               else if(readCol>=imgReader.nrOfCol())
		//                 readCol=imgReader.nrOfCol()-1;
		writeBuffer.push_back(flag_opt[0]);
	      }
	      else{
                bool valid=true;
                if(mask_opt[0]&&extent_opt[0]!=""){
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
                  writeBuffer.push_back(flag_opt[0]);
                else{
                  double theScale=(scale_opt.size()>1)?scale_opt[iband]:scale_opt[0];
                  double theOffset=(offset_opt.size()>1)?offset_opt[iband]:offset_opt[0];
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
                    writeBuffer.push_back((readCol-0.5-lowerCol)*(readBuffer[upperCol-startCol]*theScale+theOffset)+(1-readCol+0.5+lowerCol)*(readBuffer[lowerCol-startCol]*theScale+theOffset));
                    break;
                  default:
                    readCol=static_cast<int>(readCol);
                    readCol-=startCol;//we only start reading from startCol
                    writeBuffer.push_back(readBuffer[readCol]*theScale+theOffset);
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
	assert(writeBuffer.size()==ncropcol);
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
  if(extent_opt[0]!=""&&mask_opt[0]){
    extentReader.close();
  }
  imgWriter.close();
}

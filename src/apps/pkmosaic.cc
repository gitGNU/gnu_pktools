/**********************************************************************
pkmosaic.cc: program to create mosaic geo-referenced images
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
#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"

namespace mrule{
  enum MRULE_TYPE {overwrite=0, maxndvi=1, maxband=2, minband=3, validband=4, mean=5, maxvote=6, median=7,sum=8};
}

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string>  input_opt("i", "input", "Input image file(s). If input contains multiple images, a multi-band output is created", "");
  Optionpk<string>  output_opt("o", "output", "Output image file", "");
  Optionpk<string>  projection_opt("p", "projection", "projection in EPSG format (leave blank to copy from input file, use EPSG:3035 to use European projection and to force to European grid", "");
  Optionpk<string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file", "");
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box (in geocoordinates if georef is true)", 0.0);
  Optionpk<double>  dx_opt("dx", "dx", "Output resolution in x (in meter) (0.0: keep original resolution)", 0.0);
  Optionpk<double>  dy_opt("dy", "dy", "Output resolution in y (in meter) (0.0: keep original resolution)", 0.0);
  Optionpk<int>  band_opt("b", "band", "band index to crop (-1: crop all bands)", -1);
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image", "");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]", "INTERLEAVE=BAND");
  Optionpk<short>  flag_opt("f", "flag", "Flag value to put in image if out of bounds.", 0);
  Optionpk<unsigned short>  resample_opt("r", "resample", "Resampling method (0: nearest neighbour, 1: bi-linear interpolation).", 0);
  Optionpk<string>  description_opt("\0", "description", "Set image description", "");
  Optionpk<string> mrule_opt("m", "mrule", "Mosaic rule for mosaic (overwrite, maxndvi, maxband, minband, validband, mean, maxvote (only for byte images), median, sum", "overwrite");
  Optionpk<int> ruleBand_opt("rb", "rband", "band index used for the rule (for ndvi, use --ruleBand=redBand --ruleBand=nirBand", 0);
  Optionpk<int> validBand_opt("vb", "validBand", "valid band index(es)", 0);
  Optionpk<double> invalid_opt("t", "invalid", "invalid value for valid band", 0);
  Optionpk<double> minValue_opt("min", "min", "flag values smaller or equal to this value as invalid.", -99999999);
  Optionpk<double> maxValue_opt("max", "max", "flag values larger or equal to this value as invalid.", 99999999);
  Optionpk<bool> file_opt("file", "file", "write number of observations for each pixels as additional layer in mosaic", false);
  Optionpk<short> weight_opt("w", "weight", "Weights (type: short) for the mosaic, use one weight for each input file in same order as input files are provided). Use value 1 for equal weights.", 1);
  Optionpk<short> class_opt("c", "class", "classes for multi-band output image: each band represents the number of observations for one specific class. Use value 0 for no multi-band output image).", 0);
  Optionpk<bool>  verbose_opt("v", "verbose", "verbose", false);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    extent_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    flag_opt.retrieveOption(argc,argv);
    resample_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
    mrule_opt.retrieveOption(argc,argv);
    ruleBand_opt.retrieveOption(argc,argv);
    validBand_opt.retrieveOption(argc,argv);
    invalid_opt.retrieveOption(argc,argv);
    minValue_opt.retrieveOption(argc,argv);
    maxValue_opt.retrieveOption(argc,argv);
    file_opt.retrieveOption(argc,argv);
    weight_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
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

  std::map<std::string, mrule::MRULE_TYPE> mruleMap;
  //initialize mruleMap
  enum MRULE_TYPE {overwrite=0, maxndvi=1, maxband=2, minband=3, validband=4, mean=5, maxvote=6, median=7,sum=8};

  mruleMap["overwrite"]=mrule::overwrite;
  mruleMap["maxndvi"]=mrule::maxndvi;
  mruleMap["maxband"]=mrule::maxband;
  mruleMap["minband"]=mrule::minband;
  mruleMap["validband"]=mrule::validband;
  mruleMap["mean"]=mrule::mean;
  mruleMap["maxvote"]=mrule::maxvote;
  mruleMap["median"]=mrule::median;
  mruleMap["sum"]=mrule::sum;

  while(invalid_opt.size()<validBand_opt.size())
    invalid_opt.push_back(invalid_opt[0]);
  RESAMPLE theResample;
  switch(resample_opt[0]){
  case(BILINEAR):
    theResample=BILINEAR;
    if(verbose_opt[0])
      cout << "resampling: bilinear interpolation" << endl;
    break;
  default:
    theResample=NEAR;
    if(verbose_opt[0])
      cout << "resampling: nearest neighbour" << endl;
    break;
  }
  
  int nband=0;
  int nwriteBand=0;
  int writeBand=0;
  vector<short> bands;
  
  //get bounding box
  double maxLRX=0;
  double maxULY=0;
  double minULX=0;
  double minLRY=0;
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

  double dx=dx_opt[0];
  double dy=dy_opt[0];
  //get bounding box from extentReader if defined
  ImgReaderOgr extentReader;
  if(extent_opt[0]!=""){
    extentReader.open(extent_opt[0]);
    if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
       cerr << "Error: could not get extent from " << extent_opt[0] << endl;
       exit(1);
      }
    else if(verbose_opt[0])
      cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << endl;
  }

  ImgReaderGdal imgReader;
  string theProjection="";
  GDALColorTable* theColorTable=NULL;
  string imageType;
  bool init=false;
  for(int ifile=0;ifile<input_opt.size();++ifile){
    try{
      imgReader.open(input_opt[ifile]);
    }
    catch(string errorstring){
      cerr << errorstring << " " << input_opt[ifile] << endl;
    }
    if(colorTable_opt[0]=="")
      if(imgReader.getColorTable())
        theColorTable=(imgReader.getColorTable()->Clone());
    if(projection_opt[0]=="")
      theProjection=imgReader.getProjection();
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=imgReader.getInterleave();
      option_opt.push_back(theInterleave);
    }

    if((ulx_opt[0]||uly_opt[0]||lrx_opt[0]||lry_opt[0])&&(!imgReader.covers(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
      if(verbose_opt[0])
	cout << input_opt[ifile] << " not within bounding box, skipping..." << endl;
      imgReader.close();
      continue;
    }
    double theULX, theULY, theLRX, theLRY;
    imgReader.getBoundingBox(theULX,theULY,theLRX,theLRY);
    if(verbose_opt[0])
      cout << "Bounding Box (ULX ULY LRX LRY): " << fixed << setprecision(6) << theULX << " " << theULY << " " << theLRX << " " << theLRY << endl;
    if(!init){
      if(verbose_opt[0]){
        switch(mruleMap[mrule_opt[0]]){
        default:
        case(mrule::overwrite):
          cout << "Mosaic rule: overwrite" << endl;
          break;
        case(mrule::maxndvi):
          cout << "Mosaic rule: max ndvi" << endl;
          break;
        case(mrule::maxband):
          cout << "Mosaic rule: max band" << endl;
          break;
        case(mrule::minband):
          cout << "Mosaic rule: min band" << endl;
          break;
        case(mrule::validband):
          cout << "Mosaic rule: valid band" << endl;
          break;
        case(mrule::mean):
          cout << "Mosaic rule: mean value" << endl;
          break;
        case(mrule::maxvote):
          cout << "Mosaic rule: max voting (only for byte images)" << endl;
          break;
        case(mrule::median):
          cout << "Mosaic rule: median" << endl;
          break;
        case(mrule::sum):
          cout << "Mosaic rule: sum" << endl;
          break;
        }
      }
      if(band_opt[0]>=0){
	nband=band_opt.size();
        bands.resize(band_opt.size());
        for(int iband=0;iband<band_opt.size();++iband){
          bands[iband]=band_opt[iband];
          assert(bands[iband]<imgReader.nrOfBand());
        }
      }
      else{
	nband=imgReader.nrOfBand();
        bands.resize(nband);
        for(int iband=0;iband<nband;++iband)
          bands[iband]=iband;
      }
      assert(validBand_opt.size()==minValue_opt.size());
      assert(validBand_opt.size()==maxValue_opt.size());
      for(int iband=0;iband<validBand_opt.size();++iband){
        assert(validBand_opt[iband]>=0&&validBand_opt[iband]<nband);
        if(verbose_opt[0]){
          cout << "band " << validBand_opt[iband] << " is valid in ] " << minValue_opt[iband] << " , " << maxValue_opt[iband] << " [" << endl;
        }
      }
      //if output type not set, get type from input image
      if(theType==GDT_Unknown){
        theType=imgReader.getDataType();
        if(verbose_opt[0])
          cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
      }

      if(oformat_opt[0]!="")//default
        imageType=oformat_opt[0];
      else
        imageType=imgReader.getImageType();

      // dataType=imgReader.getDataType(0);
      if(verbose_opt[0]){
        cout << "type of data for " << input_opt[ifile] << ": " << theType << endl;
        cout << "nband: " << nband << endl;
      }
      
      maxLRX=theLRX;
      maxULY=theULY;
      minULX=theULX;
      minLRY=theLRY;
      if(!dx||!dy){
        dx=imgReader.getDeltaX();
        dy=imgReader.getDeltaY();
      }
      // imgReader.getMagicPixel(magic_x,magic_y);
      init=true;
    }
    else{
      //convert bounding box to magic coordinates
      //check uniformity magic pixel
      // double test_x,test_y;
      // imgReader.getMagicPixel(test_x,test_y);
      // if(verbose_opt[0]){
      //   cout << "magic_x, magic_y: " << magic_x << ", " << magic_y << endl;
      // }
      // assert(magic_x==test_x);
      // assert(magic_y==test_y);
      maxLRX=(theLRX>maxLRX)?theLRX:maxLRX;
      maxULY=(theULY>maxULY)?theULY:maxULY;
      minULX=(theULX<minULX)?theULX:minULX;
      minLRY=(theLRY<minLRY)?theLRY:minLRY;
    }
    imgReader.close();
  }
  if(verbose_opt[0])
    cout << "bounding box input images (ULX ULY LRX LRY): " << fixed << setprecision(6) << minULX << " " << maxULY << " " << maxLRX << " " << minLRY << endl;
  if(ulx_opt[0]||uly_opt[0]||lrx_opt[0]||lry_opt[0]){
    maxLRX=lrx_opt[0];
    maxULY=uly_opt[0];
    minULX=ulx_opt[0];
    minLRY=lry_opt[0];
  }
  
  if(projection_opt[0].find("ETRS-LAEA")!=string::npos||projection_opt[0].find("EPSG:3035")!=string::npos||projection_opt[0].find("epsg:3035")!=string::npos){
    //force to LAEA grid
    minULX=floor(minULX);
    minULX-=static_cast<int>(minULX)%(static_cast<int>(dx));
    maxULY=ceil(maxULY);
    if(static_cast<int>(maxULY)%static_cast<int>(dy))
      maxULY+=dy;
    maxULY-=static_cast<int>(maxULY)%(static_cast<int>(dy));
    maxLRX=ceil(maxLRX);
    if(static_cast<int>(maxLRX)%static_cast<int>(dx))
      maxLRX+=dx;
    maxLRX-=static_cast<int>(maxLRX)%(static_cast<int>(dx));
    minLRY=floor(minLRY);
    minLRY-=static_cast<int>(minLRY)%(static_cast<int>(dy));
  }

  if(verbose_opt[0])
    cout << "bounding box mosaic image (ULX ULY LRX LRY): " << fixed << setprecision(6) << minULX << " " << maxULY << " " << maxLRX << " " << minLRY << endl;
  //initialize image
  if(verbose_opt[0])
    cout << "initializing mosaic image..." << endl;
//   double dcol=(maxLRX-minULX+dx-1)/dx;
//   double drow=(maxULY-minLRY+dy-1)/dy;
//   int ncol=static_cast<int>(dcol);
//   int nrow=static_cast<int>(drow);

  int ncol=ceil((maxLRX-minULX)/dx);
  int nrow=ceil((maxULY-minLRY)/dy);

  if(verbose_opt[0])
    cout << "mosaic image dim (nrow x ncol): " << nrow << " x " << ncol << endl;
  ImgWriterGdal imgWriter;
  while(weight_opt.size()<input_opt.size())
    weight_opt.push_back(weight_opt[0]);
  if(verbose_opt[0]){
    std::cout << weight_opt << std::endl;
  }
  if(mruleMap[mrule_opt[0]]==mrule::maxvote){
    nwriteBand=(file_opt[0])? class_opt.size()+1:class_opt.size();
  }
  else
    nwriteBand=(file_opt[0])? bands.size()+1:bands.size();
  if(verbose_opt[0])
    cout << "open output image " << output_opt[0] << " with " << nwriteBand << " bands" << endl << flush;
  try{
    imgWriter.open(output_opt[0],ncol,nrow,nwriteBand,theType,imageType,option_opt);
  }
  catch(string error){
    cout << error << endl;
  }
  if(description_opt[0]!="")
    imgWriter.setImageDescription(description_opt[0]);
  imgWriter.setGeoTransform(minULX,maxULY,dx,dy,0,0);
  if(projection_opt[0]!=""){
    if(verbose_opt[0])
      cout << "projection: " << projection_opt[0] << endl;
    imgWriter.setProjectionProj4(projection_opt[0]);
  }
  else if(theProjection!=""){
    if(verbose_opt[0])
      cout << "projection: " << theProjection << endl;
    imgWriter.setProjection(theProjection);
  }
    
  if(colorTable_opt[0]!=""){
    assert(imgWriter.getDataType()==GDT_Byte);
    imgWriter.setColorTable(colorTable_opt[0]);
  }
  else if(theColorTable)
    imgWriter.setColorTable(theColorTable);

  //create mosaic image
  if(verbose_opt[0])
     cout << "creating mosaic image" << endl;
  Vector2d<double> writeBuffer(nband,imgWriter.nrOfCol());
  vector<short> fileBuffer(ncol);//holds the number of used files
  Vector2d<short> maxBuffer;//buffer used for maximum voting
  Vector2d<double> readBuffer(nband);
  statfactory::StatFactory stat;
  if(mruleMap[mrule_opt[0]]==maxndvi)//ndvi
    assert(ruleBand_opt.size()==2);
  if(mruleMap[mrule_opt[0]]==mrule::maxvote){//max voting
    maxBuffer.resize(imgWriter.nrOfCol(),256);//use only byte images for max voting
    for(int iclass=0;iclass<class_opt.size();++iclass)
      assert(class_opt[iclass]<maxBuffer.size());
  }
  int jb=0;
  double readRow=0;
  double readCol=0;
  double lowerCol=0;
  double upperCol=0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int irow=0;irow<imgWriter.nrOfRow();++irow){
    Vector2d< vector<double> > storeBuffer;
    vector<bool> writeValid(ncol);

    if(mruleMap[mrule_opt[0]]==mrule::mean||mruleMap[mrule_opt[0]]==mrule::median||mruleMap[mrule_opt[0]]==mrule::sum)//mean, median or (weighted) sum value
      storeBuffer.resize(nband,ncol);
    for(int icol=0;icol<imgWriter.nrOfCol();++icol){
      writeValid[icol]=false;
      fileBuffer[icol]=0;
      if(mruleMap[mrule_opt[0]]==mrule::maxvote){//max voting
        for(int iclass=0;iclass<256;++iclass)
          maxBuffer[icol][iclass]=0;
      }
      else{
        for(int iband=0;iband<nband;++iband)
          writeBuffer[iband][icol]=flag_opt[0];
      }
    }
    for(int ifile=0;ifile<input_opt.size();++ifile){
      try{
        imgReader.open(input_opt[ifile]);
      }
      catch(string error){
        cout << error << endl;
      }
      // assert(imgReader.getDataType()==theType);
      assert(imgReader.nrOfBand()>=nband);
      if(!imgReader.covers(minULX,maxULY,maxLRX,minLRY)){
        imgReader.close();
        continue;
      }
      double uli,ulj,lri,lrj;
      imgReader.geo2image(minULX+(magic_x-1.0)*imgReader.getDeltaX(),maxULY-(magic_y-1.0)*imgReader.getDeltaY(),uli,ulj);
      imgReader.geo2image(maxLRX+(magic_x-2.0)*imgReader.getDeltaX(),minLRY-(magic_y-2.0)*imgReader.getDeltaY(),lri,lrj);
      uli=floor(uli);
      ulj=floor(ulj);
      lri=floor(lri);
      lrj=floor(lrj);
        
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

      //convert irow to geo
      double x=0;
      double y=0;
      imgWriter.image2geo(0,irow,x,y);
      //lookup corresponding row for irow in this file
      imgReader.geo2image(x,y,readCol,readRow);
      if(readRow<0||readRow>=imgReader.nrOfRow()){
        imgReader.close();
        continue;
      }
      // for(int iband=0;iband<imgReader.nrOfBand();++iband){
      for(int iband=0;iband<nband;++iband){
	int readBand=(band_opt[0]<0)?iband:band_opt[iband];
        readBuffer[iband].resize(readncol);
	try{
          imgReader.readData(readBuffer[iband],GDT_Float64,startCol,endCol,readRow,readBand,theResample);
	}
	catch(string error){
	  cerr << "error reading image " << input_opt[ifile] << ": " << endl;
	  throw;
	}
      }
        
      for(int ib=0;ib<ncol;++ib){
        assert(imgWriter.image2geo(ib,irow,x,y));
        //lookup corresponding row for irow in this file
        imgReader.geo2image(x,y,readCol,readRow);
        if(readCol<0||readCol>=imgReader.nrOfCol())
          continue;
        double val_current=0;
        double val_new=0;
        bool readValid=true;
        switch(resample_opt[0]){
        case(BILINEAR):
          lowerCol=readCol-0.5;
          lowerCol=static_cast<int>(lowerCol);
          upperCol=readCol+0.5;
          upperCol=static_cast<int>(upperCol);
          if(lowerCol<0)
            lowerCol=0;
          if(upperCol>=imgReader.nrOfCol())
            upperCol=imgReader.nrOfCol()-1;
          for(int vband=0;vband<validBand_opt.size();++vband){
            val_new=(readCol-0.5-lowerCol)*readBuffer[validBand_opt[vband]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[validBand_opt[vband]][lowerCol-startCol];
            if(val_new<=minValue_opt[vband]||val_new>=maxValue_opt[vband]||val_new==invalid_opt[vband]){
              readValid=false;
              break;
            }
          }
          break;
        default:
          readCol=static_cast<int>(readCol);
          for(int vband=0;vband<validBand_opt.size();++vband){
            val_new=readBuffer[validBand_opt[vband]][readCol-startCol];
            if(val_new<=minValue_opt[vband]||val_new>=maxValue_opt[vband]||val_new==invalid_opt[vband]){
              readValid=false;
              break;
            }
          }
          break;
        }
	if(readValid){
          if(writeValid[ib]){
            int iband=0;
	    switch(mruleMap[mrule_opt[0]]){
	    case(mrule::maxndvi):{//max ndvi
              double red_current=writeBuffer[ruleBand_opt[0]][ib];
              double nir_current=writeBuffer[ruleBand_opt[1]][ib];
	      double ndvi_current=0;
              if(red_current+nir_current>0&&red_current>=0&&nir_current>=0)
                ndvi_current=(nir_current-red_current)/(nir_current+red_current);
	      double ndvi_new=0;
              double red_new=0;
              double nir_new=0;
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                red_new=(readCol-0.5-lowerCol)*readBuffer[ruleBand_opt[0]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ruleBand_opt[0]][lowerCol-startCol];
                nir_new=(readCol-0.5-lowerCol)*readBuffer[ruleBand_opt[1]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ruleBand_opt[1]][lowerCol-startCol];
                if(red_new+nir_new>0&&red_new>=0&&nir_new>=0)
                  ndvi_new=(nir_new-red_new)/(nir_new+red_new);
                if(ndvi_new>=ndvi_current){
                  for(iband=0;iband<nband;++iband){
                    val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
                    writeBuffer[iband][ib]=val_new;
                  }
                  // fileBuffer[ib]=ifile;
                  ++fileBuffer[ib];
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                red_new=readBuffer[ruleBand_opt[0]][readCol-startCol];
                nir_new=readBuffer[ruleBand_opt[1]][readCol-startCol];
                if(red_new+nir_new>0&&red_new>=0&&nir_new>=0)
                  ndvi_new=(nir_new-red_new)/(nir_new+red_new);
                if(ndvi_new>=ndvi_current){
                  for(iband=0;iband<nband;++iband){
                    val_new=readBuffer[iband][readCol-startCol];
                    writeBuffer[iband][ib]=val_new;
                  }
                  ++fileBuffer[ib];
                  // fileBuffer[ib]=ifile;
                }
                break;
              }
	      break;
            }
	    case(mrule::maxband):
            case(mrule::minband):
            case(mrule::validband)://max,min,valid band
              val_current=writeBuffer[ruleBand_opt[0]][ib];
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                val_new=(readCol-0.5-lowerCol)*readBuffer[ruleBand_opt[0]][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[ruleBand_opt[0]][lowerCol-startCol];
                val_new*=weight_opt[ifile];
                if((mruleMap[mrule_opt[0]]==mrule::maxband&&val_new>val_current)||(mruleMap[mrule_opt[0]]==mrule::minband&&val_new<val_current)||(mruleMap[mrule_opt[0]]==mrule::validband)){//&&val_new>minValue_opt[0]&&val_new<maxValue_opt[0])){
                  for(iband=0;iband<nband;++iband){
                    val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
                    val_new*=weight_opt[ifile];
                    writeBuffer[iband][ib]=val_new;
                  }
                  // fileBuffer[ib]=ifile;
                  ++fileBuffer[ib];
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                val_new=readBuffer[ruleBand_opt[0]][readCol-startCol];
                val_new*=weight_opt[ifile];
                if((mruleMap[mrule_opt[0]]==mrule::maxband&&val_new>val_current)||(mruleMap[mrule_opt[0]]==mrule::minband&&val_new<val_current)||(mruleMap[mrule_opt[0]]==mrule::validband)){//&&val_new>minValue_opt[0]&&val_new<maxValue_opt[0])){
                  for(iband=0;iband<nband;++iband){
                    val_new=readBuffer[iband][readCol-startCol];
                    val_new*=weight_opt[ifile];
                    writeBuffer[iband][ib]=val_new;
                  }
                  // fileBuffer[ib]=ifile;
                  ++fileBuffer[ib];
                }
                break;
              }
	      break;
            case(mrule::maxvote)://max voting (only for Byte images)
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
                  // ++(maxBuffer[ib][val_new]);
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[iband][readCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
		}
                break;
	      }
              break;
            case(mrule::mean)://mean value
	    case(mrule::median)://median value
	    case(mrule::sum)://sum value
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                  assert(ifile>0);
                  // assert(weight_opt[ifile]>=0);
                  // assert(storeBuffer[iband][ib].back()>=0);
                }
                break;
              }
              ++fileBuffer[ib];
	      break;
	    case(mrule::overwrite):
	    default:
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              }
            // fileBuffer[ib]=ifile;
            ++fileBuffer[ib];
            break;
	    }
          }
	  else{
            writeValid[ib]=true;//readValid was true
            int iband=0;
	    switch(mruleMap[mrule_opt[0]]){
            case(mrule::mean):
            case(mrule::median):
            case(mrule::sum):
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  storeBuffer[iband][ib].push_back(val_new);
                }
                break;
              }
              ++fileBuffer[ib];
              break;
            case(mrule::maxvote):
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
                  // ++(maxBuffer[ib][val_new]);
		}
                break;
              default:
                readCol=static_cast<int>(readCol);
                for(iband=0;iband<nband;++iband){
		  val_new=readBuffer[iband][readCol-startCol];
		  maxBuffer[ib][val_new]=maxBuffer[ib][val_new]+weight_opt[ifile];
		}
                  // ++(maxBuffer[ib][val_new]);
                break;
              }
              break;
            default:
              switch(resample_opt[0]){
              case(BILINEAR):
                lowerCol=readCol-0.5;
                lowerCol=static_cast<int>(lowerCol);
                upperCol=readCol+0.5;
                upperCol=static_cast<int>(upperCol);
                if(lowerCol<0)
                  lowerCol=0;
                if(upperCol>=imgReader.nrOfCol())
                  upperCol=imgReader.nrOfCol()-1;
                for(iband=0;iband<nband;++iband){
                  val_new=(readCol-0.5-lowerCol)*readBuffer[iband][upperCol-startCol]+(1-readCol+0.5+lowerCol)*readBuffer[iband][lowerCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              default:
                readCol=static_cast<int>(readCol);
                for(iband=0;iband<nband;++iband){
                  val_new=readBuffer[iband][readCol-startCol];
                  val_new*=weight_opt[ifile];
                  writeBuffer[iband][ib]=val_new;
                }
                break;
              }
              // fileBuffer[ib]=ifile;
              ++fileBuffer[ib];
              break;
            }
          }
        }
      }
      imgReader.close();
    }
    if(mruleMap[mrule_opt[0]]==mrule::maxvote){
      vector<short> classBuffer(imgWriter.nrOfCol());
      if(class_opt.size()>1){
        for(int iclass=0;iclass<class_opt.size();++iclass){
          for(int icol=0;icol<imgWriter.nrOfCol();++icol)
            classBuffer[icol]=maxBuffer[icol][class_opt[iclass]];
          try{
            imgWriter.writeData(classBuffer,GDT_Int16,irow,iclass);
          }
          catch(string error){
            cerr << "error writing image file " << output_opt[0] << ": " << error << endl;
            throw;
          }
        }
      }
      else{
        for(int icol=0;icol<imgWriter.nrOfCol();++icol){
          vector<short>::iterator maxit=maxBuffer[icol].begin();
          maxit=stat.max(maxBuffer[icol],maxBuffer[icol].begin(),maxBuffer[icol].end());
          writeBuffer[0][icol]=distance(maxBuffer[icol].begin(),maxit);
          fileBuffer[icol]=*(maxit);
        }
        try{
          imgWriter.writeData(writeBuffer[0],GDT_Float64,irow,0);
          if(file_opt[0])
            imgWriter.writeData(fileBuffer,GDT_Int16,irow,1);
        }
        catch(string error){
          cerr << "error writing image file " << output_opt[0] << ": " << error << endl;
          throw;
        }
      }
    }
    else{
      for(int iband=0;iband<bands.size();++iband){
        // assert(writeBuffer[bands[iband]].size()==imgWriter.nrOfCol());
        assert(writeBuffer[iband].size()==imgWriter.nrOfCol());
        for(int icol=0;icol<imgWriter.nrOfCol();++icol){
          switch(mruleMap[mrule_opt[0]]){
          case(mrule::mean):
            assert(storeBuffer[bands[iband]][icol].size()==fileBuffer[icol]);
            if(storeBuffer[bands[iband]][icol].size())
              writeBuffer[iband][icol]=stat.mean(storeBuffer[bands[iband]][icol]);
            break;
          case(mrule::median):
            assert(storeBuffer[bands[iband]][icol].size()==fileBuffer[icol]);
            if(storeBuffer[bands[iband]][icol].size())
              writeBuffer[iband][icol]=stat.median(storeBuffer[bands[iband]][icol]);
            break;
          case(mrule::sum)://sum
            assert(storeBuffer[bands[iband]][icol].size()==fileBuffer[icol]);
            if(storeBuffer[bands[iband]][icol].size())
              writeBuffer[iband][icol]=stat.sum(storeBuffer[bands[iband]][icol]);
            break;
          default:
            break;
          }
        }
        try{
          imgWriter.writeData(writeBuffer[iband],GDT_Float64,irow,iband);
        }
        catch(string error){
          cerr << "error writing image file " << output_opt[0] << ": " << error << endl;
          throw;
        }
      }
      if(file_opt[0]){
        try{
          imgWriter.writeData(fileBuffer,GDT_Int16,irow,bands.size());
        }
        catch(string error){
          cerr << "error writing image file " << output_opt[0] << ": " << error << endl;
          throw;
        }
      }
    }
    progress=static_cast<float>(irow+1.0)/imgWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  imgWriter.close();
}
  

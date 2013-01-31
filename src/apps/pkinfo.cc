/**********************************************************************
pkinfo.cc: program to retrieve information from raster images
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
#include <sstream>
#include <list>
#include "base/Optionpk.h"
#include "algorithms/Egcs.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgReaderOgr.h"

int main(int argc, char *argv[])
{
  Optionpk<std::string> input_opt("i","input","Input image file");
  Optionpk<bool>  bbox_opt("bb", "bbox", "Shows bounding box ", false,1);
  Optionpk<bool>  bbox_te_opt("te", "te", "Shows bounding box in GDAL format: xmin ymin xmax ymax ", false,1);
  Optionpk<bool>  centre_opt("c", "centre", "Image centre in projected X,Y coordinates ", false,1);
  Optionpk<bool>  colorTable_opt("ct", "colourtable", "Shows colour table ", false,1);
  Optionpk<bool>  samples_opt("s", "samples", "Number of samples in image ", false,1);
  Optionpk<bool>  lines_opt("l", "lines", "Number of lines in image ", false,1);
  Optionpk<bool>  nband_opt("nb", "nband", "Show number of bands in image", false,1);
  Optionpk<short>  band_opt("b", "band", "Band specific information", 0,1);
  Optionpk<bool>  dx_opt("dx", "dx", "Gets resolution in x (in m)", false,1);
  Optionpk<bool>  dy_opt("dy", "dy", "Gets resolution in y (in m)", false,1);
  Optionpk<bool>  minmax_opt("mm", "minmax", "Shows min and max value of the image ", false,1);
  Optionpk<bool>  stat_opt("stat", "stat", "Shows statistics (min,max, mean and stdDev of the image)", false,1);
  Optionpk<double>  min_opt("min", "min", "Sets minimum for histogram");
  Optionpk<double>  max_opt("max", "max", "Sets maximum for histogram");
  Optionpk<bool>  relative_opt("rel", "rel", "Calculates relative histogram in percentage", false,1);
  Optionpk<bool>  projection_opt("p", "projection", "Shows projection of the image ", false,1);
  Optionpk<bool>  geo_opt("geo", "geo", "Gets geotransform  ", false,1);
  Optionpk<bool>  interleave_opt("il", "interleave", "Shows interleave ", false,1);
  Optionpk<bool>  filename_opt("f", "filename", "Shows image filename ", false,1);
  Optionpk<bool>  cover_opt("cov", "cover", "Image covers bounding box (or x and y pos) if printed to std out ", false,1);
  Optionpk<double>  x_opt("x", "xpos", "x pos");
  Optionpk<double>  y_opt("y", "ypos", "y pos");
  Optionpk<bool>  read_opt("r", "read", "Reads row y (in projected coordinates if geo option is set, otherwise in image coordinates, 0 based)",false,1);
  Optionpk<bool>  refpixel_opt("ref", "ref", "Gets reference pixel (lower left corner of centre of gravity pixel)", false,1);
  Optionpk<bool>  driver_opt("of", "oformat", "Gets driver description ", false,1);
  Optionpk<std::string>  extent_opt("e", "extent", "Gets boundary from extent from polygons in vector file");
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box");
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box");
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box");
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box");
  Optionpk<bool>  hist_opt("hist", "hist", "Calculates histogram. Use --rel for a relative histogram output. ", false,1);
  Optionpk<short>  nbin_opt("nbin", "nbin", "Number of bins used in histogram. Use 0 for all input values as integers", 0,1);
  Optionpk<bool>  type_opt("ot", "otype", "Returns data type", false,1);
  Optionpk<bool>  description_opt("d", "description", "Returns image description", false,1);
  Optionpk<bool>  metadata_opt("meta", "meta", "Shows meta data ", false,1);
  Optionpk<double> nodata_opt("nodata", "nodata", "Sets no data value(s) for calculations (flags in input image)");

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    bbox_opt.retrieveOption(argc,argv);
    bbox_te_opt.retrieveOption(argc,argv);
    centre_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    samples_opt.retrieveOption(argc,argv);
    lines_opt.retrieveOption(argc,argv);
    nband_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    minmax_opt.retrieveOption(argc,argv);
    stat_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    max_opt.retrieveOption(argc,argv);
    relative_opt.retrieveOption(argc,argv);
    projection_opt.retrieveOption(argc,argv);
    geo_opt.retrieveOption(argc,argv);
    interleave_opt.retrieveOption(argc,argv);
    filename_opt.retrieveOption(argc,argv);
    cover_opt.retrieveOption(argc,argv);
    x_opt.retrieveOption(argc,argv);
    y_opt.retrieveOption(argc,argv);
    read_opt.retrieveOption(argc,argv);
    refpixel_opt.retrieveOption(argc,argv);
    driver_opt.retrieveOption(argc,argv);
    extent_opt.retrieveOption(argc,argv);
    ulx_opt.retrieveOption(argc,argv);
    uly_opt.retrieveOption(argc,argv);
    lrx_opt.retrieveOption(argc,argv);
    lry_opt.retrieveOption(argc,argv);
    hist_opt.retrieveOption(argc,argv);
    nbin_opt.retrieveOption(argc,argv);
    type_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
    metadata_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }
  //for union
  double maxLRX=0;
  double maxULY=0;
  double minULX=0;
  double minLRY=0;

  //for intersect
  double minLRX=0;
  double minULY=0;
  double maxULX=0;
  double maxLRY=0;
  
  ImgReaderGdal imgReader;
  for(int ifile=0;ifile<input_opt.size();++ifile){
    imgReader.open(input_opt[ifile]);
    for(int inodata=0;inodata<nodata_opt.size();++inodata){
      if(!inodata)
        imgReader.GDALSetNoDataValue(nodata_opt[0],band_opt[0]);//only single no data can be set in GDALRasterBand (used for ComputeStatistics)
      imgReader.pushNoDataValue(nodata_opt[inodata]);
    }
    if(filename_opt[0])
      std::cout << " --input " << input_opt[ifile] << " ";
    if(centre_opt[0]){
      double theX, theY;
      imgReader.getCentrePos(theX,theY);
      std::cout << setprecision(12) << " -x " << theX << " -y " << theY << " ";
    }
    if(refpixel_opt[0]){
      assert(band_opt[0]<imgReader.nrOfBand());
      Egcs egcs;
      double refX,refY;
      //get centre of reference (centre of gravity) pixel in image
      imgReader.getRefPix(refX,refY,band_opt[0]);
      cout << setprecision(12) << "-rx " << refX << " -ry " << refY << endl;
      egcs.setLevel(egcs.res2level(imgReader.getDeltaX()));
      // unsigned short theLevel=egcs.getLevel(imgReader.getDeltaX());
      // egcs.setLevel(theLevel);
      //cout << "cell code at level " << egcs.getLevel() << " (resolution is " << egcs.getResolution() << "): " << egcs.geo2cell(refX,refY) << endl;
    }
    if(bbox_opt[0]||bbox_te_opt[0]){
      double theULX, theULY, theLRX, theLRY;
      imgReader.getBoundingBox(theULX,theULY,theLRX,theLRY);
      if(bbox_te_opt[0])
        std::cout << setprecision(12) << "-te " << theULX << " " << theLRY << " " << theLRX << " " << theULY;
      else
        std::cout << setprecision(12) << "--ulx=" << theULX << " --uly=" << theULY << " --lrx=" << theLRX << " --lry=" << theLRY << " ";
      if(!ifile){
	maxLRX=theLRX;
	maxULY=theULY;
	minULX=theULX;
	minLRY=theLRY;

	minLRX=theLRX;
	minULY=theULY;
	maxULX=theULX;
	maxLRY=theLRY;
      }
      else{
	maxLRX=(theLRX>maxLRX)?theLRX:maxLRX;
	maxULY=(theULY>maxULY)?theULY:maxULY;
	minULX=(theULX<minULX)?theULX:minULX;
	minLRY=(theLRY<minLRY)?theLRY:minLRY;

	minLRX=(theLRX<minLRX)?theLRX:minLRX;
	minULY=(theULY<minULY)?theULY:minULY;
	maxULX=(theULX>maxULX)?theULX:maxULX;
	maxLRY=(theLRY>maxLRY)?theLRY:maxLRY;
      }
    }
    if(dx_opt[0])
      std::cout << "--dx " << imgReader.getDeltaX() << " ";
    if(dy_opt[0])
      std::cout << "--dy " << imgReader.getDeltaY() << " ";
    if(cover_opt[0]){
      //get bounding box from extentReader if defined
      ImgReaderOgr extentReader;
      if(extent_opt.size()){
        extentReader.open(extent_opt[0]);
        if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
          cerr << "Error: could not get extent from " << extent_opt[0] << std::endl;
          exit(1);
        }
        // std::cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << std::endl;
      }
      double theULX, theULY, theLRX, theLRY;
      imgReader.getBoundingBox(theULX,theULY,theLRX,theLRY);
      if((ulx_opt.size()||uly_opt.size()||lrx_opt.size()||lry_opt.size())&&(imgReader.covers(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0])))
	std::cout << " -i " << input_opt[ifile] << " ";
      else if(imgReader.covers(x_opt[0],y_opt[0]))
	std::cout << " -i " << input_opt[ifile] << " ";

    }
    else if(ulx_opt.size()||uly_opt.size()||lrx_opt.size()||lry_opt.size()){
      double ulx,uly,lrx,lry;
      imgReader.getBoundingBox(ulx,uly,lrx,lry);
      if(ulx_opt.size())
        std::cout << " --ulx=" << fixed << ulx << " ";
      if(uly_opt.size())
        std::cout << " --uly=" << fixed << uly << " ";
      if(lrx_opt.size())
        std::cout << " --lrx=" << fixed << lrx << " ";
      if(lry_opt.size())
        std::cout << " --lry=" << fixed << lry << " ";
    }
    if(colorTable_opt[0]){
      GDALColorTable* colorTable=imgReader.getColorTable();
      if(colorTable!=NULL){
        for(int index=0;index<colorTable->GetColorEntryCount();++index){
          GDALColorEntry sEntry=*(colorTable->GetColorEntry(index));
          std::cout << index << " " << sEntry.c1 << " " << sEntry.c2 << " " << sEntry.c3 << " " << sEntry.c4 << std::endl;
        }
      }
      else
        std::cout << "--ct none ";
    }
    if(samples_opt[0])
      std::cout << "--samples " << imgReader.nrOfCol() << " ";
    if(lines_opt[0])
      std::cout << "--rows " << imgReader.nrOfRow() << " ";
    if(nband_opt[0])
      std::cout << "--nband " << imgReader.nrOfBand() << " ";
    double minValue=0;
    double maxValue=0;
    double meanValue=0;
    double stdDev=0;
    if(stat_opt[0]){
      assert(band_opt[0]<imgReader.nrOfBand());
      GDALProgressFunc pfnProgress;
      void* pProgressData;
      GDALRasterBand* rasterBand;
      rasterBand=imgReader.getRasterBand(band_opt[0]);
      rasterBand->ComputeStatistics(0,&minValue,&maxValue,&meanValue,&stdDev,pfnProgress,pProgressData);
      std::cout << "--min " << minValue << " --max " << maxValue << " --mean " << meanValue << " --stdDev " << stdDev << " ";
    }

    if(minmax_opt[0]){
      assert(band_opt[0]<imgReader.nrOfBand());
      imgReader.getMinMax(minValue,maxValue,band_opt[0],true);
      std::cout << "--min " << minValue << " --max " << maxValue << " ";
    }
    if(hist_opt[0]){
      assert(band_opt[0]<imgReader.nrOfBand());
      imgReader.getMinMax(minValue,maxValue,band_opt[0]);
      if(min_opt.size())
        minValue=min_opt[0];
      if(max_opt.size())
        maxValue=max_opt[0];
      int nbin=nbin_opt[0];
      if(nbin_opt[0]==0)
	nbin=maxValue-minValue+1;
      assert(nbin>0);
      std::vector<unsigned long int> output(nbin);
      unsigned long int nsample=0;
      unsigned long int ninvalid=0;
      std::vector<double> lineBuffer(imgReader.nrOfCol());
      for(int i=0;i<nbin;output[i++]=0);
      for(int irow=0;irow<imgReader.nrOfRow();++irow){
	imgReader.readData(lineBuffer,GDT_Float64,irow,band_opt[0]);
	for(int icol=0;icol<imgReader.nrOfCol();++icol){
          if(imgReader.isNoData(lineBuffer[icol]))
            ++ninvalid;
          else if(lineBuffer[icol]>maxValue)
            ++ninvalid;
          else if(lineBuffer[icol]<minValue)
            ++ninvalid;
	  else if(lineBuffer[icol]==maxValue)
	    ++output[nbin-1];
	  // else if(static_cast<double>(lineBuffer[icol]-minValue)/(maxValue-minValue)*nbin>=nbin){
          //   //test
          //   std::cout << "..." << lineBuffer[icol] << std::endl;
	  //   ++output[nbin-1];
          // }
	  else
	    ++output[static_cast<int>(static_cast<double>(lineBuffer[icol]-minValue)/(maxValue-minValue)*nbin)];
	}
      }
      nsample=imgReader.nrOfCol()*imgReader.nrOfRow()-ninvalid;
      std::cout.precision(10);
      for(int bin=0;bin<nbin;++bin){
	nsample+=output[bin];
        if(output[bin]>0){
          std::cout << (maxValue-minValue)*bin/(nbin-1)+minValue << " ";
          if(relative_opt[0])
            std::cout << 100.0*static_cast<double>(output[bin])/static_cast<double>(nsample) << std::endl;
          else
            std::cout << static_cast<double>(output[bin])  << std::endl;
        }
      }
    }
    else{
      int minCol,minRow;
      if(min_opt.size()){
        assert(band_opt[0]<imgReader.nrOfBand());
        std::cout << "--min " << imgReader.getMin(minCol, minRow,band_opt[0]);
      }
      if(max_opt.size()){
        assert(band_opt[0]<imgReader.nrOfBand());
        assert(band_opt[0]<imgReader.nrOfBand());
        std::cout << "--max " << imgReader.getMax(minCol, minRow,band_opt[0]);
      }
    }
    if(projection_opt[0]){
      if(imgReader.isGeoRef())
        std::cout << "--projection " << imgReader.getProjection() << " ";
      else
        std::cout << " --projection none" << " ";
    }
    if(geo_opt[0]&&!read_opt[0]){
      double ulx,uly,deltaX,deltaY,rot1,rot2;
      imgReader.getGeoTransform(ulx,uly,deltaX,deltaY,rot1,rot2);
      std::cout << " --geo " << setprecision(12) << ulx << " " << uly << " " << deltaX << " " << deltaY << " " << rot1 << " " << rot2 << " ";
    }
    if(interleave_opt[0]){
      std::cout << " --interleave " << imgReader.getInterleave() << " ";
    }
    if(type_opt[0]){
      std::cout << "--otype " << GDALGetDataTypeName(imgReader.getDataType(band_opt[0])) << " ";
      // std::cout << " -ot " << GDALGetDataTypeName(imgReader.getDataType(band_opt[0])) << " (" << static_cast<short>(imgReader.getDataType(band_opt[0])) << ")" << std::endl;
    }
    if(description_opt[0]){
      // try{
      // 	std::cout << "image description: " << imgReader.getImageDescription() << std::endl;
      // }
      // catch(...){
      // 	std::cout << "catched" << std::endl;
      // }
      list<std::string> metaData;
      imgReader.getMetadata(metaData);
      list<std::string>::const_iterator lit=metaData.begin();
      std::cout << " --description ";
      while(lit!=metaData.end())
      	std::cout << *(lit++) << " ";
    }
    if(metadata_opt[0]){
      std::cout << "Metadata: " << std::endl;
      list<std::string> lmeta;
      imgReader.getMetadata(lmeta);
      list<std::string>::const_iterator lit=lmeta.begin();
      while(lit!=lmeta.end()){
        std::cout << *lit << std::endl;
        ++lit;
      }
//       char** cmetadata=imgReader.getMetadata();
//       while(*cmetadata!=NULL){
//         std::cout << *(cmetadata) << std::endl;
//         ++cmetadata;
//       }
    }
    if(read_opt[0]){
      int nband=band_opt.size();
      if(band_opt[0]<0)
        nband=imgReader.nrOfBand();
      std::cout.precision(12);
      for(int iband=0;iband<nband;++iband){
        unsigned short theBand=(band_opt[0]<0)? iband : band_opt[iband];
        std::vector<float> rowBuffer;//buffer will be resized in readdata
        for(int iy=0;iy<y_opt.size();++iy){
	  double theRow=y_opt[iy];
	  int ncol=(x_opt.size())? x_opt.size() : imgReader.nrOfCol();
          for(int ix=0;ix<ncol;++ix){
	    double theCol=ix;
	    if(x_opt.size()){
	      if(geo_opt[0])
		imgReader.geo2image(x_opt[ix],y_opt[iy],theCol,theRow);
	      else
		theCol=x_opt[ix];
	    }
            assert(theRow>=0);
            assert(theRow<imgReader.nrOfRow());
            imgReader.readData(rowBuffer,GDT_Float32, static_cast<int>(theRow), theBand);
	    assert(theCol<rowBuffer.size());
	    std::cout << rowBuffer[static_cast<int>(theCol)] << " ";
	  }
          std::cout << std::endl;
        }
      }
    }
    if(driver_opt[0])
      std::cout << " --oformat " << imgReader.getDriverDescription() << " ";
    imgReader.close();
  }
  if((bbox_opt[0]||bbox_te_opt[0])&&input_opt.size()>1){
    if(bbox_te_opt[0])
      std::cout << setprecision(12) << "-te " << minULX << " " << minLRY << " " << maxLRX << " " << maxULY;
    else
      std::cout << "union bounding box: " << setprecision(12) << "--ulx=" << minULX << " --uly=" << maxULY << " --lrx=" << maxLRX << " --lry=" << minLRY << std::endl;
    if(maxULX<minLRX&&minULY>maxLRY){
      if(bbox_te_opt[0])
        std::cout << "intersect bounding box: " << setprecision(12) << "-te " << maxULX << " " << maxLRY << " " << minLRX << " --lry=" << minULY << std::endl;
      else
        std::cout << "intersect bounding box: " << setprecision(12) << "--ulx=" << maxULX << " --uly=" << minULY << " --lrx=" << minLRX << " --lry=" << maxLRY << std::endl;
    }
    else
      std::cout << "no intersect" << std::endl;
  }
  if(!input_opt.size())
    std::cerr << "No input file provided (use option -i). Use pkinfo --help for help information" << std::endl;
  else if(!read_opt[0]&&!hist_opt[0])
    std::cout << std::endl;
}

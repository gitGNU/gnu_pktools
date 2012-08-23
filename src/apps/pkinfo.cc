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
  Optionpk<std::string> input_opt("i","input","input image file","");
  Optionpk<bool>  bbox_opt("bb", "bbox", "show bounding box ", false);
  Optionpk<bool>  bbox_te_opt("te", "te", "show bounding box in GDAL format: xmin ymin xmax ymax ", false);
  Optionpk<bool>  centre_opt("c", "centre", "Image centre in projected X,Y coordinates ", false);
  Optionpk<bool>  colorTable_opt("ct", "colourtable", "show colour table ", false);
  Optionpk<bool>  samples_opt("s", "samples", "Number of samples in image ", false);
  Optionpk<bool>  lines_opt("l", "lines", "Number of lines in image ", false);
  Optionpk<bool>  nband_opt("nb", "nband", "Show number of bands in image", false);
  Optionpk<short>  band_opt("b", "band", "Band specific information", 0);
  Optionpk<bool>  dx_opt("dx", "dx", "get resolution in x (in m)", false);
  Optionpk<bool>  dy_opt("dy", "dy", "get resolution in y (in m)", false);
  Optionpk<bool>  minmax_opt("mm", "minmax", "Show min and max value of the image ", false);
  Optionpk<bool>  stat_opt("stat", "stat", "Show statistics (min,max, mean and stdDev of the image ", false);
  Optionpk<double>  min_opt("m", "min", "Show min value of the image (or set minimum for histogram)", 0);
  Optionpk<double>  max_opt("M", "max", "Show max value of the image (or set maximum for histogram)", 0);
  Optionpk<bool>  relative_opt("rel", "rel", "Calculate relative histogram in percentage", false);
  Optionpk<bool>  projection_opt("p", "projection", "Show projection of the image ", false);
  Optionpk<bool>  geo_opt("geo", "geo", "get geotransform:  ", false);
  Optionpk<bool>  interleave_opt("il", "interleave", "Show interleave ", false);
  Optionpk<bool>  filename_opt("f", "filename", "Image filename ", false);
  Optionpk<bool>  cover_opt("cov", "cover", "Image covers bounding box (or x and y pos) if printed to std out ", false);
  Optionpk<double>  x_opt("x", "xpos", "x pos", -1);
  Optionpk<double>  y_opt("y", "ypos", "y pos", -1);
  Optionpk<bool>  read_opt("r", "read", "read row y (in projected coordinates if geo option is set, otherwise in image coordinates, 0 based)", 0);
  Optionpk<bool>  refpixel_opt("ref", "ref", "get reference pixel (lower left corner of centre of gravity pixel)", false);
  Optionpk<bool>  driver_opt("of", "oformat", "get driver description ", false);
  Optionpk<std::string>  extent_opt("e", "extent", "get boundary from extent from polygons in vector file", "");
  Optionpk<double>  ulx_opt("ulx", "ulx", "Upper left x value bounding box (0)", 0.0);
  Optionpk<double>  uly_opt("uly", "uly", "Upper left y value bounding box (0)", 0.0);
  Optionpk<double>  lrx_opt("lrx", "lrx", "Lower right x value bounding box (0)", 0.0);
  Optionpk<double>  lry_opt("lry", "lry", "Lower right y value bounding box (0)", 0.0);
  Optionpk<bool>  hist_opt("hist", "hist", "Calculate histogram. Use --rel for a relative histogram output. ", false);
  Optionpk<short>  nbin_opt("nbin", "nbin", "Number of bins used in histogram. Use 0 for all input values as integers", 0);
  Optionpk<bool>  type_opt("ot", "otype", "Return data type", false);
  Optionpk<bool>  description_opt("d", "description", "Return image description", false);
  Optionpk<bool>  metadata_opt("meta", "meta", "Show meta data ", false);
  Optionpk<double> maskValue_opt("mask", "mask", "mask value(s) for no data to calculate reference pixel in image",0);
  
  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  if(version_opt[0]){
    std::cout << version_opt.getHelp() << std::endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }

  input_opt.retrieveOption(argc,argv);
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
  maskValue_opt.retrieveOption(argc,argv);

  if(help_opt[0]){
    std::cout << "usage: pkinfo -i imagefile [OPTIONS]" << std::endl;
    exit(0);
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
    for(int inodata=0;inodata<maskValue_opt.size();++inodata)
      imgReader.pushNoDataValue(maskValue_opt[inodata],band_opt[0]);

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
      if(extent_opt[0]!=""){
        extentReader.open(extent_opt[0]);
        if(!(extentReader.getExtent(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0]))){
          cerr << "Error: could not get extent from " << extent_opt[0] << std::endl;
          exit(1);
        }
        // std::cout << "--ulx=" << ulx_opt[0] << " --uly=" << uly_opt[0] << " --lrx=" << lrx_opt[0] << " --lry=" << lry_opt[0] << std::endl;
      }
      double theULX, theULY, theLRX, theLRY;
      imgReader.getBoundingBox(theULX,theULY,theLRX,theLRY);
      if((ulx_opt[0]||uly_opt[0]||lrx_opt[0]||lry_opt[0])&&(imgReader.covers(ulx_opt[0],uly_opt[0],lrx_opt[0],lry_opt[0])))
	std::cout << " -i " << input_opt[ifile] << " ";
      else if(imgReader.covers(x_opt[0],y_opt[0]))
	std::cout << " -i " << input_opt[ifile] << " ";

    }
    else if(ulx_opt[0]||uly_opt[0]||lrx_opt[0]||lry_opt[0]){
      double ulx,uly,lrx,lry;
      imgReader.getBoundingBox(ulx,uly,lrx,lry);
      if(ulx_opt[0])
        std::cout << " --ulx=" << fixed << ulx << " ";
      if(uly_opt[0])
        std::cout << " --uly=" << fixed << uly << " ";
      if(lrx_opt[0])
        std::cout << " --lrx=" << fixed << lrx << " ";
      if(lry_opt[0])
        std::cout << " --lry=" << fixed << lry << " ";
    }
    if(colorTable_opt[0]){
      GDALColorTable* colorTable=imgReader.getColorTable();
      for(int index=0;index<colorTable->GetColorEntryCount();++index){
        GDALColorEntry sEntry=*(colorTable->GetColorEntry(index));
        std::cout << index << " " << sEntry.c1 << " " << sEntry.c2 << " " << sEntry.c3 << " " << sEntry.c4 << std::endl;
      }
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
      if(min_opt[0]==max_opt[0])
        imgReader.getMinMax(minValue,maxValue,band_opt[0]);
      else{
        minValue=min_opt[0];
        maxValue=max_opt[0];
      }
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
	  if(lineBuffer[icol]>maxValue)
            ++ninvalid;
          else if(lineBuffer[icol]<minValue)
            ++ninvalid;
	  else if(lineBuffer[icol]==maxValue)
	    ++output[nbin-1];
	  else if(static_cast<double>(lineBuffer[icol]-minValue)/(maxValue-minValue)*nbin>=nbin)
	    ++output[nbin-1];
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
      if(min_opt[0]){
        assert(band_opt[0]<imgReader.nrOfBand());
        std::cout << "--min " << imgReader.getMin(minCol, minRow,band_opt[0]);
      }
      if(max_opt[0]){
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
          for(int ix=0;ix<x_opt.size();++ix){
	    double theCol=x_opt[iy];
	    if(geo_opt[0])
	      imgReader.geo2image(x_opt[ix],y_opt[iy],theCol,theRow);
            assert(theRow>=0);
            assert(theRow<imgReader.nrOfRow());
            imgReader.readData(rowBuffer,GDT_Float32, static_cast<int>(theRow), theBand);
            if(x_opt[ix]>=0){
              assert(theCol<rowBuffer.size());
              std::cout << rowBuffer[static_cast<int>(theCol)] << " ";
            }
            else{
              for(int i=0;i<rowBuffer.size();++i)
                std::cout << rowBuffer[i] << " ";
            }
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
  //    std::cout << "bounding box mosaic (ULX ULY LRX LRY): " << minULX << " " << maxULY << " " << maxLRX << " " << minLRY << std::endl;
  if(!read_opt[0])
    std::cout << std::endl;
}

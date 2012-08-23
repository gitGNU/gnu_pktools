/**********************************************************************
pklas2img.cc: program to create (e.g., DEM) raster image from las files
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
#include <iostream>
#include "Optionpk.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "fileclasses/FileReaderLas.h"
#include "algorithms/Histogram.h"
#include "algorithms/Filter2d.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

int main(int argc,char **argv) {
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<string> input_opt("i", "input", "Input las file", "");
  // Optionpk<string> mask_opt("m", "mask", "mask image file", "");
  // Optionpk<short> invalid_opt("t", "invalid", "Mask value(s) where image is invalid. Use multiple values for a single mask. Default value is 0", 0);
  Optionpk<short> flag_opt("f", "flag", "Flag value(s) to put in image if not valid. Use as many flags as invalid options (Default is 0)", 0);
  Optionpk<string> attribute_opt("n", "name", "names of the attribute to select: [intensity|return|nreturn|z]", "z");
  Optionpk<bool> disc_opt("circ", "circular", "circular disc kernel for dilation and erosion (default is false)", false);
  Optionpk<double> maxSlope_opt("s", "maxSlope", "Maximum slope used for morphological filtering (default is 0)", 0.0);
  Optionpk<double> hThreshold_opt("ht", "maxHeight", "initial and maximum height threshold for progressive morphological filtering (e.g., -ht 0.2 -ht 2.5)", 0.2);
  Optionpk<short> maxIter_opt("\0", "maxIter", "Maximum number of iterations in post filter (default is 100)", 100.0);
  Optionpk<short> nbin_opt("nb", "nbin", "Number of percentile bins for calculating profile (=number of output bands) (default is 10)", 10.0);
  Optionpk<unsigned short> returns_opt("r", "returns", "number(s) of returns to include (use -r -1 for last return only). Default is 0 (include all returns)", 0);
  Optionpk<string> composite_opt("c", "composite", "composite for multiple points in cell (min, max, median, mean, sum, first, last, profile). Default is last (overwrite cells with latest point", "last");
  Optionpk<string> filter_opt("fir", "filter", "filter las points (single,multiple,all). Default is all", "all");
  Optionpk<string> postFilter_opt("pf", "pfilter", "filter las points (etew_min,promorph (progressive morphological filter),bunting (adapted promorph),open,close,none) . Default is none", "none");
  Optionpk<short> dimx_opt("\0", "dimX", "Dimension X of postFilter (default is 3)", 3);
  Optionpk<short> dimy_opt("\0", "dimY", "Dimension Y of postFilter (default is 3)", 3);
  Optionpk<string> output_opt("o", "output", "Output image file", "");
  Optionpk<string> projection_opt("p", "projection", "projection in EPSG code, e.g., EPSG:3035 (Default is no projection)", "");
  Optionpk<double> ulx_opt("\0", "ulx", "Upper left x value bounding box (in geocoordinates if georef is true). Default is 0: read from input file", 0.0);
  Optionpk<double> uly_opt("\0", "uly", "Upper left y value bounding box (in geocoordinates if georef is true). Default is 0: read from input file", 0.0);
  Optionpk<double> lrx_opt("\0", "lrx", "Lower right x value bounding box (in geocoordinates if georef is true). Default is 0: read from input file", 0.0);
  Optionpk<double> lry_opt("\0", "lry", "Lower right y value bounding box (in geocoordinates if georef is true). Default is 0: read from input file", 0.0);
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "Byte");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image", "GTiff");
  Optionpk<string> option_opt("co", "co", "options: NAME=VALUE [-co COMPRESS=LZW] [-co INTERLEAVE=BAND]", "INTERLEAVE=BAND");
  Optionpk<double> dx_opt("dx", "dx", "Output resolution in x (in meter) (default is 0.0: keep original resolution)", 0.0);
  Optionpk<double> dy_opt("dy", "dy", "Output resolution in y (in meter) (default is 0.0: keep original resolution)", 0.0);
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)", "");
  Optionpk<short> verbose_opt("v", "verbose", "verbose (default is 0)", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);

  if(version_opt[0]){
    cout << version_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }

  input_opt.retrieveOption(argc,argv);
  // mask_opt.retrieveOption(argc,argv);
  // invalid_opt.retrieveOption(argc,argv);
  flag_opt.retrieveOption(argc,argv);
  attribute_opt.retrieveOption(argc,argv);
  disc_opt.retrieveOption(argc,argv);
  maxSlope_opt.retrieveOption(argc,argv);
  hThreshold_opt.retrieveOption(argc,argv);
  maxIter_opt.retrieveOption(argc,argv);
  nbin_opt.retrieveOption(argc,argv);
  returns_opt.retrieveOption(argc,argv);
  composite_opt.retrieveOption(argc,argv);
  filter_opt.retrieveOption(argc,argv);
  postFilter_opt.retrieveOption(argc,argv);
  dimx_opt.retrieveOption(argc,argv);
  dimy_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  projection_opt.retrieveOption(argc,argv);
  ulx_opt.retrieveOption(argc,argv);
  uly_opt.retrieveOption(argc,argv);
  lrx_opt.retrieveOption(argc,argv);
  lry_opt.retrieveOption(argc,argv);
  otype_opt.retrieveOption(argc,argv);
  oformat_opt.retrieveOption(argc,argv);
  option_opt.retrieveOption(argc,argv);
  dx_opt.retrieveOption(argc,argv);
  dy_opt.retrieveOption(argc,argv);
  colorTable_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

  if(help_opt[0]){
    cout << version_opt.getHelp() << endl;
    exit(0);
  }


  GDALAllRegister();

  double dfComplete=0.0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;

  Vector2d<vector<float> > inputData;//row,col,point
   
  ImgReaderGdal maskReader;
  ImgWriterGdal outputWriter;
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
    if(theType==GDT_Unknown)
      cout << "Unknown output pixel type: " << otype_opt[0] << endl;
    else
      cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
  }

  double maxLRX=0;
  double maxULY=0;
  double minULX=0;
  double minLRY=0;

  unsigned long int totalPoints=0;
  unsigned long int nPoints=0;
  unsigned long int ipoint=0;
  for(int iinput=0;iinput<input_opt.size();++iinput){
    assert(input_opt[iinput].find(".las")!=string::npos);
    FileReaderLas lasReader;
    try{
      lasReader.open(input_opt[iinput]);
    }
    catch(string errorString){
      cout << errorString << endl;
      exit(1);
    }
    nPoints=lasReader.getPointCount();
    totalPoints+=nPoints;

    if(ulx_opt[0]>=lrx_opt[0]||uly_opt[0]<=lry_opt[0]){
      double ulx,uly,lrx,lry;
      lasReader.getExtent(ulx,uly,lrx,lry);
      lrx+=dx_opt[0];//pixel coordinates are referenced to upper left corner (las coordinates are centres)
      lry-=dy_opt[0];//pixel coordinates are referenced to upper left corner (las coordinates are centres)
      if(ulx>=lrx){
        ulx=ulx-dx_opt[0]/2.0;
        lrx=ulx+dx_opt[0]/2.0;
      }
      if(uly<=lry){
        uly=lry+dy_opt[0]/2.0;
        lry=lry-dy_opt[0]/2.0;
      }
      if(maxLRX>minULX){
        maxLRX=(lrx>maxLRX)?lrx:maxLRX;
        maxULY=(uly>maxULY)?uly:maxULY;
        minULX=(ulx<minULX)?ulx:minULX;
        minLRY=(lry<minLRY)?lry:minLRY;
      }
      else{//initialize
        maxLRX=lrx;
        maxULY=uly;
        minULX=ulx;
        minLRY=lry;
      }        
    }
    else{
      maxLRX=lrx_opt[0];
      maxULY=uly_opt[0];
      minULX=ulx_opt[0];
      minLRY=lry_opt[0];
    }
    lasReader.close();
  }
  if(verbose_opt[0]){
    std::cout << setprecision(12) << "--ulx=" << minULX << " --uly=" << maxULY << " --lrx=" << maxLRX << " --lry=" << minLRY << std::endl;
    std::cout << "total number of points before filtering: " << totalPoints << std::endl;
    std::cout << "filter set to " << filter_opt[0] << std::endl;
    std::cout << "postFilter set to " << postFilter_opt[0] << std::endl;
  }
  int ncol=ceil(maxLRX-minULX)/dx_opt[0];//number of columns in outputGrid
  int nrow=ceil(maxULY-minLRY)/dy_opt[0];//number of rows in outputGrid
  //todo: multiple bands
  int nband=(composite_opt[0]=="profile")? nbin_opt[0] : 1;
  if(output_opt[0]==""){
    cerr << "Error: no output file defined" << endl;
    exit(1);
  }
  if(verbose_opt[0])
    cout << "opening output file " << output_opt[0] << endl;
  outputWriter.open(output_opt[0],ncol,nrow,nband,theType,oformat_opt[0],option_opt);
  //set projection
  outputWriter.setGeoTransform(minULX,maxULY,dx_opt[0],dy_opt[0],0,0);
  if(projection_opt[0]!=""){
    string projectionString=outputWriter.setProjectionProj4(projection_opt[0]);
    if(verbose_opt[0])
      cout << "projection: " << projectionString << endl;
  }
  if(!outputWriter.isGeoRef())
    cout << "Warning: output image " << output_opt[0] << " is not georeferenced!" << endl;
  if(colorTable_opt[0]!="")
    outputWriter.setColorTable(colorTable_opt[0]);

  inputData.clear();
  inputData.resize(nrow,ncol);
  std::cout << "Reading " << input_opt.size() << " las files" << std::endl;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int iinput=0;iinput<input_opt.size();++iinput){
    FileReaderLas lasReader;
    try{
      lasReader.open(input_opt[iinput]);
    }
    catch(string errorString){
      cout << errorString << endl;
      exit(1);
    }
    //set bounding filter
    // lasReader.addBoundsFilter(minULX,maxULY,maxLRX,minLRY);
    //set returns filter
    if(returns_opt[0])
      lasReader.addReturnsFilter(returns_opt);
    lasReader.setFilters();

    if(attribute_opt[0]!="z"){
      vector<boost::uint16_t> returnsVector;
      vector<string>::iterator ait=attribute_opt.begin();
      while(ait!=attribute_opt.end()){
        if(*ait=="intensity"){
          if(verbose_opt[0])
            std::cout << "writing intensity" << std::endl;
          ++ait;
        }
        else if(*ait=="return"){
          if(verbose_opt[0])
            std::cout << "writing return number" << std::endl;
          ++ait;
        }
        else if(*ait=="nreturn"){
          if(verbose_opt[0])
            std::cout << "writing number of returns" << std::endl;
          ++ait;
        }
        else
          attribute_opt.erase(ait);
      }
    }
    liblas::Point thePoint;
    while(lasReader.readNextPoint(thePoint)){
      progress=static_cast<float>(ipoint)/totalPoints;
      pfnProgress(progress,pszMessage,pProgressArg);
      if(verbose_opt[0]>1)
        cout << "reading point " << ipoint << endl;
      if(thePoint.GetX()<minULX||thePoint.GetX()>=maxLRX||thePoint.GetY()>=maxULY||thePoint.GetY()<minLRY)
        continue;
      if((filter_opt[0]=="single")&&(thePoint.GetNumberOfReturns()!=1))
        continue;
      if((filter_opt[0]=="multiple")&&(thePoint.GetNumberOfReturns()<2))
        continue;
      double dcol,drow;
      outputWriter.geo2image(thePoint.GetX(),thePoint.GetY(),dcol,drow);
      int icol=static_cast<int>(dcol);
      int irow=static_cast<int>(drow);
      assert(irow>=0);
      assert(irow<nrow);
      assert(icol>=0);
      assert(icol<ncol);
      if(attribute_opt[0]=="z")
        inputData[irow][icol].push_back(thePoint.GetZ());
      else if(attribute_opt[0]=="intensity")
        inputData[irow][icol].push_back(thePoint.GetIntensity());
      else if(attribute_opt[0]=="return")
        inputData[irow][icol].push_back(thePoint.GetReturnNumber());
      else if(attribute_opt[0]=="nreturn")
        inputData[irow][icol].push_back(thePoint.GetNumberOfReturns());
      else{
        std::string errorString="attribute not supported";
        throw(errorString);
      }
      ++ipoint;
    }
    if(verbose_opt[0])
      std::cout << "number of points: " << ipoint << std::endl;
    lasReader.close();
  }
  progress=1;
  pfnProgress(progress,pszMessage,pProgressArg);

  std::cout << "processing LiDAR points" << std::endl;
  progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  Histogram hist;
  //fill in inputData in outputData
  Vector2d<float> outputData(nrow,ncol);
  if(composite_opt[0]=="profile"){
    assert(postFilter_opt[0]=="none");
    // for(int iband=0;iband<nband;++iband)
      // outputProfile[iband].resize(nrow,ncol);
  }
  for(int irow=0;irow<nrow;++irow){
    Vector2d<float> outputProfile(nband,ncol);
    for(int icol=0;icol<ncol;++icol){
      std::vector<float> profile;
      if(!inputData[irow][icol].size())
        outputData[irow][icol]=(static_cast<float>((flag_opt[0])));
      else{
        Histogram hist;
        if(composite_opt[0]=="min")
          outputData[irow][icol]=hist.min(inputData[irow][icol]);
        else if(composite_opt[0]=="max")
          outputData[irow][icol]=hist.max(inputData[irow][icol]);
        else if(composite_opt[0]=="median")
          outputData[irow][icol]=hist.median(inputData[irow][icol]);
        else if(composite_opt[0]=="mean")
          outputData[irow][icol]=hist.mean(inputData[irow][icol]);
        else if(composite_opt[0]=="sum")
          outputData[irow][icol]=hist.sum(inputData[irow][icol]);
        else if(composite_opt[0]=="first")
          outputData[irow][icol]=inputData[irow][icol][0];
        else if(composite_opt[0]=="last")
          outputData[irow][icol]=inputData[irow][icol].back();
        else if(composite_opt[0]=="profile"){
          if(inputData[irow][icol].size()<2){
            for(int iband=0;iband<nband;++iband)
              outputProfile[iband][icol]=static_cast<float>(flag_opt[0]);
            continue;
          }
          float min=0;
          float max=0;
          hist.minmax(inputData[irow][icol],inputData[irow][icol].begin(),inputData[irow][icol].end(),min,max);
          if(verbose_opt[0])
            std::cout << "min,max: " << min << "," << max << std::endl;
          if(max>min){
            hist.percentiles(inputData[irow][icol],inputData[irow][icol].begin(),inputData[irow][icol].end(),profile,nband,min,max);
            assert(profile.size()==nband);
            for(int iband=0;iband<nband;++iband)
              outputProfile[iband][icol]=profile[iband];
          }
          else{
            for(int iband=0;iband<nband;++iband)
              outputProfile[iband][icol]=max;
          }
        }
        else{
          std::cout << "Error: composite_opt " << composite_opt[0] << " not supported" << std::endl;
          exit(2);
        }
      }
    }
    if(composite_opt[0]=="profile"){
      for(int iband=0;iband<nband;++iband){
        // assert(outputProfile[iband].size()==outputWriter.nrOfRow());
        assert(outputProfile[iband].size()==outputWriter.nrOfCol());
        try{
          outputWriter.writeData(outputProfile[iband],theType,irow,iband);
        }
        catch(std::string errorString){
          cout << errorString << endl;
          exit(1);
        }
      }
    }
    progress=static_cast<float>(irow)/outputWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  progress=1;
  pfnProgress(progress,pszMessage,pProgressArg);
  inputData.clear();//clean up memory
  //apply post filter
  std::cout << "Applying post processing filter: " << postFilter_opt[0] << std::endl;
  if(postFilter_opt[0]=="etew_min"){
    if(composite_opt[0]!="min")
      std::cout << "Warning: composite option is not set to min!" << std::endl;
    //Elevation Threshold with Expand Window (ETEW) Filter (p.73 frmo Airborne LIDAR Data Processing and Analysis Tools ALDPAT 1.0)
    //first iteration is performed assuming only minima are selected using options -fir all -c min
    unsigned long int nchange=1;
    //increase cells and thresholds until no points from the previous iteration are discarded.
    int dimx=dimx_opt[0];
    int dimy=dimy_opt[0];
    Filter2d::Filter2d morphFilter;
    morphFilter.setNoValue(0);
    Vector2d<float> currentOutput=outputData;
    int iteration=1;
    while(nchange&&iteration<maxIter_opt[0]){
      double hThreshold=maxSlope_opt[0]*dimx;
      Vector2d<float> newOutput;
      nchange=morphFilter.morphology(currentOutput,newOutput,Filter2d::ERODE,dimx,dimy,disc_opt[0],hThreshold);
      currentOutput=newOutput;
      dimx+=2;//change from theory: originally double cellCize
      dimy+=2;//change from theory: originally double cellCize
      std::cout << "iteration " << iteration << ": " << nchange << " pixels changed" << std::endl;
      ++iteration;
    }
    outputData=currentOutput;
  }    
  else if(postFilter_opt[0]=="promorph"||postFilter_opt[0]=="bunting"){
    if(composite_opt[0]!="min")
      std::cout << "Warning: composite option is not set to min!" << std::endl;
    assert(hThreshold_opt.size()>1);
    //Progressive morphological filter tgrs2003_zhang vol41 pp 872-882
    //first iteration is performed assuming only minima are selected using options -fir all -c min
    //increase cells and thresholds until no points from the previous iteration are discarded.
    int dimx=dimx_opt[0];
    int dimy=dimy_opt[0];
    Filter2d::Filter2d theFilter;
    theFilter.setNoValue(0);
    Vector2d<float> currentOutput=outputData;
    double hThreshold=hThreshold_opt[0];
    int iteration=1;
    while(iteration<maxIter_opt[0]){
      std::cout << "iteration " << iteration << " with window size " << dimx << " and dh_max: " << hThreshold << std::endl;
      Vector2d<float> newOutput;
      try{
        theFilter.morphology(outputData,currentOutput,Filter2d::ERODE,dimx,dimy,disc_opt[0],maxSlope_opt[0]);
        theFilter.morphology(currentOutput,outputData,Filter2d::DILATE,dimx,dimy,disc_opt[0],maxSlope_opt[0]);
      //   if(postFilter_opt[0]=="bunting"){//todo: implement doit in Filter2d on Vector2d
      //     theFilter.doit(outputData,currentOutput,Filter2d::MEDIAN,dimx,dimy,1,disc_opt[0]);
      // filter2d.doit(input,output,Filter2d::MEDIAN,dimX_opt[0],dimY_opt[0],down_opt[0],disc_opt[0]);

      //     outputData=currentOutput;
      //   }
      }
      catch(std::string errorString){
        cout << errorString << endl;
        exit(1);
      }
      int newdimx=2*(dimx-1)+1;
      int newdimy=2*(dimy-1)+1;//from PE&RS vol 71 pp313-324
      hThreshold=hThreshold_opt[0]+maxSlope_opt[0]*(newdimx-dimx)*dx_opt[0];
      dimx=newdimx;
      dimy=newdimy;
      if(hThreshold>hThreshold_opt[1])
        hThreshold=hThreshold_opt[1];
      ++iteration;
    }
    outputData=currentOutput;
  }    
  else if(postFilter_opt[0]=="open"){
    if(composite_opt[0]!="min")
      std::cout << "Warning: composite option is not set to min!" << std::endl;
    Filter2d::Filter2d morphFilter;
    morphFilter.setNoValue(0);
    Vector2d<float> filterInput=outputData;
    try{
      morphFilter.morphology(outputData,filterInput,Filter2d::ERODE,dimx_opt[0],dimy_opt[0],disc_opt[0],maxSlope_opt[0]);
      morphFilter.morphology(filterInput,outputData,Filter2d::DILATE,dimx_opt[0],dimy_opt[0],disc_opt[0],maxSlope_opt[0]);
    }
    catch(std::string errorString){
      cout << errorString << endl;
      exit(1);
    }
  }
  else if(postFilter_opt[0]=="close"){
    if(composite_opt[0]!="max")
      std::cout << "Warning: composite option is not set to max!" << std::endl;
    Filter2d::Filter2d morphFilter;
    morphFilter.setNoValue(0);
    Vector2d<float> filterInput=outputData;
    try{
      morphFilter.morphology(outputData,filterInput,Filter2d::DILATE,dimx_opt[0],dimy_opt[0],disc_opt[0],maxSlope_opt[0]);
      morphFilter.morphology(filterInput,outputData,Filter2d::ERODE,dimx_opt[0],dimy_opt[0],disc_opt[0],maxSlope_opt[0]);
    }
    catch(std::string errorString){
      cout << errorString << endl;
      exit(1);
    }
  }
  if(composite_opt[0]!="profile"){
    //write output file
    std::cout << "writing output raster file" << std::endl;
    progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    for(int irow=0;irow<nrow;++irow){
      try{
        assert(outputData.size()==outputWriter.nrOfRow());
        assert(outputData[0].size()==outputWriter.nrOfCol());
        outputWriter.writeData(outputData[irow],theType,irow,0);
      }
      catch(std::string errorString){
        cout << errorString << endl;
        exit(1);
      }
      progress=static_cast<float>(irow)/outputWriter.nrOfRow();
      pfnProgress(progress,pszMessage,pProgressArg);
    }
  }
  progress=1;
  pfnProgress(progress,pszMessage,pProgressArg);
  if(verbose_opt[0])
    std::cout << "closing lasReader" << std::endl;
  outputWriter.close();
}

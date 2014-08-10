/**********************************************************************
pkfilterdem.cc: program to post filter raster images created with pklas2img
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
#include <iostream>
#include <string>
#include "base/Optionpk.h"
#include "base/Vector2d.h"
#include "algorithms/Filter2d.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;
/*------------------
  Main procedure
  ----------------*/
int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i","input","input image file");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<std::string> tmpdir_opt("tmp", "tmp", "Temporary directory","/tmp",2);
  Optionpk<bool> disc_opt("circ", "circular", "circular disc kernel for dilation and erosion", false);
  Optionpk<string> postFilter_opt("f", "filter", "post processing filter: vito, etew_min, promorph (progressive morphological filter),open,close).");
  Optionpk<double> dim_opt("dim", "dim", "maximum filter kernel size (optionally you can set both initial and maximum filter kernel size", 3);
  Optionpk<double> maxSlope_opt("st", "st", "slope threshold used for morphological filtering. Use a low values to remove more height objects in flat terrains", 0.0);
  Optionpk<double> hThreshold_opt("ht", "ht", "initial height threshold for progressive morphological filtering. Use low values to remove more height objects. Optionally, a maximum height threshold can be set via a second argument (e.g., -ht 0.2 -ht 2.5 sets an initial threshold at 0.2 m and caps the threshold at 2.5 m).", 0.2);
  Optionpk<short> minChange_opt("minchange", "minchange", "Stop iterations when no more pixels are changed than this threshold.", 0);
  Optionpk<std::string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image","");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate). Empty string: inherit from input image");
  Optionpk<string>  colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid). Use none to ommit color table");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> nodata_opt("nodata", "nodata", "nodata value(s) for smoothnodata filter");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    tmpdir_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    postFilter_opt.retrieveOption(argc,argv);
    dim_opt.retrieveOption(argc,argv);
    maxSlope_opt.retrieveOption(argc,argv);
    hThreshold_opt.retrieveOption(argc,argv);
    minChange_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
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

  ImgReaderGdal input;
  ImgWriterGdal outputWriter;
  if(input_opt.empty()){
    cerr << "Error: no input file selected, use option -i" << endl;
    exit(1);
  }
  if(output_opt.empty()){
    cerr << "Error: no outputWriter file selected, use option -o" << endl;
    exit(1);
  }
  if(postFilter_opt.empty()){
    cerr << "Error: no filter selected, use option -f" << endl;
    exit(1);
  }
  input.open(input_opt[0]);
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
  if(theType==GDT_Unknown)
    theType=input.getDataType();

  if(verbose_opt[0])
    std::cout << std::endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;

  string imageType=input.getImageType();
  if(oformat_opt.size())
    imageType=oformat_opt[0];

  if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
    string theInterleave="INTERLEAVE=";
    theInterleave+=input.getInterleave();
    option_opt.push_back(theInterleave);
  }

  if(verbose_opt[0])
    cout << "opening output file " << output_opt[0] << endl;
  outputWriter.open(output_opt[0],input.nrOfCol(),input.nrOfRow(),1,theType,imageType,option_opt);
  //set projection
  outputWriter.setProjection(input.getProjection());
  outputWriter.copyGeoTransform(input);
  if(colorTable_opt.size())
    outputWriter.setColorTable(colorTable_opt[0]);   

  //set nodata value
  if(nodata_opt.size()){
      for(int iband=0;iband<outputWriter.nrOfBand();++iband)
	outputWriter.GDALSetNoDataValue(nodata_opt[0],iband);
  }

  Vector2d<double> inputData(input.nrOfRow(),input.nrOfCol());
  Vector2d<double> outputData(outputWriter.nrOfRow(),outputWriter.nrOfCol());
  Vector2d<double> tmpData(outputWriter.nrOfRow(),outputWriter.nrOfCol());
  input.readDataBlock(inputData,GDT_Float64,0,inputData.nCols()-1,0,inputData.nRows()-1);

  //apply post filter
  std::cout << "Applying post processing filter: " << postFilter_opt[0] << std::endl;

  // const char* pszMessage;
  // void* pProgressArg=NULL;
  // GDALProgressFunc pfnProgress=GDALTermProgress;
  // double progress=0;
  // pfnProgress(progress,pszMessage,pProgressArg);

  //make sure dim_opt contains initial [0] and maximum [1] kernel sizes in this order
  if(dim_opt.size()<2)
    dim_opt.insert(dim_opt.begin(),3);
  if(dim_opt[0]>dim_opt[1]){
    dim_opt.insert(dim_opt.begin(),dim_opt[1]);
    dim_opt.erase(dim_opt.begin()+2);
  }

  filter2d::Filter2d theFilter;
  if(nodata_opt.size()){
    for(int inodata=0;inodata<nodata_opt.size();++inodata)
      theFilter.pushNoDataValue(nodata_opt[inodata]);
  }

  unsigned long int nchange=1;
  if(postFilter_opt[0]=="vito"){
    //todo: fill empty pixels
    // hThreshold_opt.resize(4);
    // hThreshold_opt[0]=0.7;
    // hThreshold_opt[1]=0.3;
    // hThreshold_opt[2]=0.1;
    // hThreshold_opt[2]=-0.2;
    vector<int> nlimit(4);
    nlimit[0]=2;
    nlimit[1]=3;
    nlimit[2]=4;
    nlimit[2]=2;
    //init finalMask
    for(int irow=0;irow<tmpData.nRows();++irow)
      for(int icol=0;icol<tmpData.nCols();++icol)
	tmpData[irow][icol]=1;
    for(int iheight=0;iheight<hThreshold_opt.size();++iheight){
      if(verbose_opt[0])
	cout << "height: " << hThreshold_opt[iheight] << endl;
      //todo:replace with binary mask (or short) -> adapt template with T1,T2 in Filter2d
      Vector2d<double> tmpMask(input.nrOfRow(),input.nrOfCol());
      for(int irow=0;irow<tmpMask.nRows();++irow)
	for(int icol=0;icol<tmpMask.nCols();++icol)
	  tmpMask[irow][icol]=1;//1=surface, 0=terrain
      if(verbose_opt[0])
	cout << "filtering NWSE" << endl;
      //from here
      // Vector2d<double> tmpDSM(inputData);
      // int dimX=dim_opt[0];
      // int dimY=dim_opt[0];
      // assert(dimX);
      // assert(dimY);
      // statfactory::StatFactory stat;
      // Vector2d<double> inBuffer(dimY,inputData.nCols());
      // if(tmpData.size()!=inputData.nRows())
      // 	tmpData.resize(inputData.nRows());
      // int indexI=0;
      // int indexJ=inputData.nRows()-1;
      // // int indexJ=0;
      // //initialize last half of inBuffer
      // for(int j=-(dimY-1)/2;j<=dimY/2;++j){
      // 	for(int i=0;i<inputData.nCols();++i)
      // 	  inBuffer[indexJ][i]=tmpDSM[abs(j)][i];
      // 	--indexJ;
      // 	// ++indexJ;
      // }
      // for(int y=tmpDSM.nRows()-1;y>=0;--y){
      // 	if(y){//inBuffer already initialized for y=0
      // 	  //erase first line from inBuffer
      // 	  inBuffer.erase(inBuffer.end()-1);
      // 	  // inBuffer.erase(inBuffer.begin());
      // 	  //read extra line and push back to inBuffer if not out of bounds
      // 	  if(y+dimY/2<tmpDSM.nRows()){
      // 	    //allocate buffer
      // 	    // inBuffer.push_back(inBuffer.back());
      // 	    inBuffer.insert(inBuffer.begin(),*(inBuffer.begin()));
      // 	    for(int i=0;i<tmpDSM.nCols();++i) 
      // 	      inBuffer[0][i]=tmpDSM[y-dimY/2][i];
      // 	      // inBuffer[inBuffer.size()-1][i]=tmpDSM[y+dimY/2][i];
      // 	  }
      // 	  else{
      // 	    int over=y+dimY/2-tmpDSM.nRows();
      // 	    int index=(inBuffer.size()-1)-over;
      // 	    assert(index>=0);
      // 	    assert(index<inBuffer.size());
      // 	    inBuffer.push_back(inBuffer[index]);
      // 	  }
      // 	}
      // 	for(int x=tmpDSM.nCols()-1;x>=0;--x){
      // 	  double centerValue=inBuffer[(dimY-1)/2][x];
      // 	  //test
      // 	  cout << "pixel ("  << x << "," << y << "): " << centerValue << endl;
      // 	  short nmasked=0;
      // 	  std::vector<double> neighbors;
      // 	  for(int j=-(dimY-1)/2;j<=dimY/2;++j){
      // 	    for(int i=-(dimX-1)/2;i<=dimX/2;++i){
      // 	      indexI=x+i;
      // 	      //check if out of bounds
      // 	      if(indexI<0)
      // 		indexI=-indexI;
      // 	      else if(indexI>=tmpDSM.nCols())
      // 		indexI=tmpDSM.nCols()-i;
      // 	      if(y+j<0)
      // 		indexJ=-j;
      // 	      else if(y+j>=tmpDSM.nRows())
      // 		indexJ=(dimY>2) ? (dimY-1)/2-j : 0;
      // 	      else
      // 		indexJ=(dimY-1)/2+j;
      // 	      double difference=(centerValue-inBuffer[indexJ][indexI]);
      // 	      //test
      // 	      cout << "centerValue-inBuffer[" << indexJ << "][" << indexI << "]=" << centerValue << " - " << inBuffer[indexJ][indexI] << " = " << difference << endl;
      // 	      if(i||j)//skip centerValue
      // 		neighbors.push_back(inBuffer[indexJ][indexI]);
      // 	      if(difference>hThreshold_opt[iheight])
      // 		++nmasked;
      // 	    }
      // 	  }
      // 	  //test
      // 	  cout << "pixel " << x << ", " << y << ": nmasked is " << nmasked << endl;
      // 	  if(nmasked<=nlimit[iheight]){
      // 	    ++nchange;
      // 	    //reset pixel in outputMask
      // 	    tmpData[y][x]=0;
      // 	    //test
      // 	    cout << "pixel " << x << ", " << y << " is ground" << endl;
      // 	  }
      // 	  else{
      // 	    //reset pixel height in tmpDSM
      // 	    sort(neighbors.begin(),neighbors.end());
      // 	    assert(neighbors.size()>1);
      // 	    inBuffer[(dimY-1)/2][x]=neighbors[1];
      // 	    //test
      // 	    cout << "pixel " << x << ", " << y << " is surface, reset DSM to " << neighbors[1] << endl;
      // 	    /* inBuffer[(dimY-1)/2][x]=stat.mymin(neighbors); */
      // 	  }
      // 	}
      // }
      //to here

      theFilter.dsm2dtm_nwse(inputData,tmpData,hThreshold_opt[iheight],nlimit[iheight],dim_opt[0]);
      if(verbose_opt[0])
      	cout << "filtering NESW" << endl;
      theFilter.dsm2dtm_nesw(inputData,tmpData,hThreshold_opt[iheight],nlimit[iheight],dim_opt[0]);
      if(verbose_opt[0])
      	cout << "filtering SENW" << endl;
      theFilter.dsm2dtm_senw(inputData,tmpData,hThreshold_opt[iheight],nlimit[iheight],dim_opt[0]);
      if(verbose_opt[0])
      	cout << "filtering SWNE" << endl;
      theFilter.dsm2dtm_swne(inputData,tmpData,hThreshold_opt[iheight],nlimit[iheight],dim_opt[0]);
    }
    outputData=tmpData;
    //todo: interpolate
    //outputData.setMask(tmpData,1,0);
  }    
  else if(postFilter_opt[0]=="etew_min"){
    //Elevation Threshold with Expand Window (ETEW) Filter (p.73 from Airborne LIDAR Data Processing and Analysis Tools ALDPAT 1.0)
    //first iteration is performed assuming only minima are selected using options -fir all -comp min
    //increase cells and thresholds until no points from the previous iteration are discarded.
    int dim=dim_opt[0];
    // theFilter.setNoValue(0);
    int iteration=1;
    while(nchange>minChange_opt[0]&&dim<=dim_opt[1]){
      double hThreshold=maxSlope_opt[0]*dim;
      nchange=theFilter.morphology(inputData,outputData,"erode",dim,dim,disc_opt[0],hThreshold);
      inputData=outputData;
      dim+=2;//change from theory: originally double cellCize
      std::cout << "iteration " << iteration << ": " << nchange << " pixels changed" << std::endl;
      ++iteration;
    }
  }    
  else if(postFilter_opt[0]=="promorph"){
    //Progressive morphological filter tgrs2003_zhang vol41 pp 872-882
    //first iteration is performed assuming only minima are selected using options -fir all -comp min
    //increase cells and thresholds until no points from the previous iteration are discarded.
    int dim=dim_opt[0];
    double hThreshold=hThreshold_opt[0];
    int iteration=1;
    while(nchange>minChange_opt[0]&&dim<=dim_opt[1]){
      std::cout << "iteration " << iteration << " with window size " << dim << " and dh_max: " << hThreshold << std::endl;
      try{
        nchange=theFilter.morphology(inputData,outputData,"erode",dim,dim,disc_opt[0],hThreshold);
        theFilter.morphology(outputData,inputData,"dilate",dim,dim,disc_opt[0],hThreshold);
	theFilter.doit(inputData,outputData,"median",dim,dim,1,disc_opt[0]);
	inputData=outputData;
      }
      catch(std::string errorString){
        cout << errorString << endl;
        exit(1);
      }
      int newdim=(dim==1)? 3: 2*(dim-1)+1;
      hThreshold=hThreshold_opt[0]+maxSlope_opt[0]*(newdim-dim)*input.getDeltaX();
      dim=newdim;
      if(hThreshold_opt.size()>1){
	if(hThreshold>hThreshold_opt[1]){
	  hThreshold=hThreshold_opt[1];
	}
      }
      std::cout << "iteration " << iteration << ": " << nchange << " pixels changed" << std::endl;
      ++iteration;
    }
  }    
  else if(postFilter_opt[0]=="open"){
    try{
      theFilter.morphology(inputData,tmpData,"erode",dim_opt[0],dim_opt[0],disc_opt[0],hThreshold_opt[0]);
      theFilter.morphology(tmpData,outputData,"dilate",dim_opt[0],dim_opt[0],disc_opt[0],hThreshold_opt[0]);
      outputData=inputData;
    }
    catch(std::string errorString){
      cout << errorString << endl;
      exit(1);
    }
  }
  else if(postFilter_opt[0]=="close"){
    try{
      theFilter.morphology(inputData,tmpData,"dilate",dim_opt[0],dim_opt[0],disc_opt[0],hThreshold_opt[0]);
      theFilter.morphology(tmpData,outputData,"erode",dim_opt[0],dim_opt[0],disc_opt[0],hThreshold_opt[0]);
    }
    catch(std::string errorString){
      cout << errorString << endl;
      exit(1);
    }
  }
  //write outputData to outputWriter
  outputWriter.writeDataBlock(outputData,GDT_Float64,0,outputData.nCols()-1,0,outputData.nRows()-1);

  // progress=1;
  // pfnProgress(progress,pszMessage,pProgressArg);
  input.close();
  outputWriter.close();
  return 0;
}

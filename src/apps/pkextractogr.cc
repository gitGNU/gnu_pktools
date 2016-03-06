/**********************************************************************
pkextractogr.cc: extract pixel values from raster image from a (vector or raster) sample
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
#include <math.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <algorithm>
#include <ctime>
#include <vector>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

/******************************************************************************/
/*! \page pkextractogr pkextractogr
 extract pixel values from raster image using a vector dataset sample
## SYNOPSIS

<code>
  Usage: pkextractogr -i input [-s sample | -rand number | -grid size] -o output
</code>

<code>

  Options: [-ln layer]* [-c class]* [-t threshold]* [-f format] [-ft fieldType] [-b band]* [-r rule]*

  Advanced options:
  [-sband band -eband band]* [-bndnodata band -srcnodata value]* [-tp threshold] [-buf value [-circ]]
</code>

\section pkextractogr_description Description

The utility pkextractogr extracts pixel values from an input raster dataset, based on the locations you provide via a sample file. Alternatively, a random sample or systematic grid of points can also be extracted. The sample can be a vector file with points or polygons. In the case of polygons, you can either extract the values for all raster pixels that are covered by the polygons, or extract a single value for each polygon such as the centroid, mean, median, etc. As output, a new copy of the vector file is created with an extra attribute for the extracted pixel value. For each raster band in the input image, a separate attribute is created. For instance, if the raster dataset contains three bands, three attributes are created (b0, b1 and b2). 

A typical usage of pkextract is to prepare a training sample for one of the classifiers implemented in pktools.

\anchor pkextractogr_rules 

Overview of the possible extraction rules:

extraction rule | output features
--------------- | ---------------
point | extract a single pixel within the polygon or on each point feature
allpoints | Extract all pixel values covered by the polygon
centroid | Extract pixel value at the centroid of the polygon
mean | Extract average of all pixel values within the polygon
stdev | Extract standard deviation of all pixel values within the polygon
median | Extract median of all pixel values within the polygon
min | Extract minimum value of all pixels within the polygon
max | Extract maximum value of all pixels within the polygon
sum | Extract sum of the values of all pixels within the polygon
mode | Extract the mode of classes within the polygon (classes must be set with the option class)
proportion | Extract proportion of class(es) within the polygon (classes must be set with the option class)
count | Extract count of class(es) within the polygon (classes must be set with the option class).
percentile | Extract percentile as defined by option perc (e.g, 95th percentile of values covered by polygon)

\section pkextract_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Raster input dataset containing band information | 
 | s      | sample               | std::string |       |OGR vector dataset with features to be extracted from input data. Output will contain features with input band information included | 
 | ln     | ln                   | std::string |       |Layer name(s) in sample (leave empty to select all) | 
 | rand   | random               | unsigned int |       |Create simple random sample of points. Provide number of points to generate | 
 | grid   | grid                 | double |       |Create systematic grid of points. Provide cell grid size (in projected units, e.g,. m) | 
 | o      | output               | std::string |       |Output sample dataset | 
 | c      | class                | int  |       |Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset. Make sure to set classes if rule is set to mode, proportion or count | 
 | t      | threshold            | float | 100   |Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold per vector sample layer | 
 | perc   | perc                 | double | 95    |Percentile value used for rule percentile | 
 | f      | f                    | std::string | SQLite |Output sample dataset format | 
 | ft     | ftype                | std::string | Real  |Field type (only Real or Integer) | 
 | b      | band                 | int  |       |Band index(es) to extract (0 based). Leave empty to use all bands | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 | r      | rule                 | std::string | centroid |Rule how to report image information per feature (only for vector sample). point (single point or at centroid if polygon), allpoints (within polygon), centroid, mean, stdev, median, proportion, count, min, max, mode, sum, percentile. | 
 | bndnodata | bndnodata            | int  | 0     |Band(s) in input image to check if pixel is valid (used for srcnodata) | 
 | srcnodata | srcnodata            | double |       |Invalid value(s) for input image | 
 | tp     | thresholdPolygon     | float |       |(absolute) threshold for selecting samples in each polygon | 
 | buf    | buffer               | short | 0     |Buffer for calculating statistics for point features (in number of pixels)  | 
 | circ   | circular             | bool | false |Use a circular disc kernel buffer (for vector point sample datasets only, use in combination with buffer option) | 

Usage: pkextract -i input [-s sample | -rand number | -grid size] -o output -r rule


Examples
========
Some examples how to use pkextract can be found \ref examples_pkextract "here"
**/

namespace rule{
  enum RULE_TYPE {point=0, mean=1, proportion=2, custom=3, min=4, max=5, mode=6, centroid=7, sum=8, median=9, stdev=10, percentile=11, count=12, allpoints=13};
}

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> image_opt("i", "input", "Raster input dataset containing band information");
  Optionpk<string> sample_opt("s", "sample", "OGR vector dataset with features to be extracted from input data. Output will contain features with input band information included. Sample image can also be GDAL raster dataset.");
  Optionpk<string> layer_opt("ln", "ln", "Layer name(s) in sample (leave empty to select all)");
  Optionpk<unsigned int> random_opt("rand", "random", "Create simple random sample of points. Provide number of points to generate");
  Optionpk<double> grid_opt("grid", "grid", "Create systematic grid of points. Provide cell grid size (in projected units, e.g,. m)");
  Optionpk<string> output_opt("o", "output", "Output sample dataset");
  Optionpk<int> class_opt("c", "class", "Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset. Make sure to set classes if rule is set to mode, proportion or count");
  Optionpk<float> threshold_opt("t", "threshold", "Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold per vector sample layer. If using raster land cover maps as a sample dataset, you can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es)", 100);
  Optionpk<double> percentile_opt("perc","perc","Percentile value used for rule percentile",95);
  Optionpk<string> ogrformat_opt("f", "f", "Output sample dataset format","SQLite");
  Optionpk<string> ftype_opt("ft", "ftype", "Field type (only Real or Integer)", "Real");
  Optionpk<int> band_opt("b", "band", "Band index(es) to extract (0 based). Leave empty to use all bands");
  Optionpk<unsigned short> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned short> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<string> rule_opt("r", "rule", "Rule how to report image information per feature (only for vector sample). point (single point within polygon), allpoints (all points within polygon), centroid, mean, stdev, median, proportion, count, min, max, mode, sum, percentile.","centroid");
  Optionpk<double> srcnodata_opt("srcnodata", "srcnodata", "Invalid value(s) for input image");
  Optionpk<int> bndnodata_opt("bndnodata", "bndnodata", "Band(s) in input image to check if pixel is valid (used for srcnodata)", 0);
  Optionpk<float> polythreshold_opt("tp", "thresholdPolygon", "(absolute) threshold for selecting samples in each polygon");
  Optionpk<short> buffer_opt("buf", "buffer", "Buffer for calculating statistics for point features (in number of pixels) ",0);
  Optionpk<bool> disc_opt("circ", "circular", "Use a circular disc kernel buffer (for vector point sample datasets only, use in combination with buffer option)", false);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose mode if > 0", 0,2);

  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  bndnodata_opt.setHide(1);
  srcnodata_opt.setHide(1);
  polythreshold_opt.setHide(1);
  percentile_opt.setHide(1);
  buffer_opt.setHide(1);
  disc_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=image_opt.retrieveOption(argc,argv);
    sample_opt.retrieveOption(argc,argv);
    layer_opt.retrieveOption(argc,argv);
    random_opt.retrieveOption(argc,argv);
    grid_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    percentile_opt.retrieveOption(argc,argv);
    ogrformat_opt.retrieveOption(argc,argv);
    ftype_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    bstart_opt.retrieveOption(argc,argv);
    bend_opt.retrieveOption(argc,argv);
    rule_opt.retrieveOption(argc,argv);
    bndnodata_opt.retrieveOption(argc,argv);
    srcnodata_opt.retrieveOption(argc,argv);
    polythreshold_opt.retrieveOption(argc,argv);
    buffer_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkextract -i input [-s sample | -rand number | -grid size] -o output" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  std::map<std::string, rule::RULE_TYPE> ruleMap;
  //initialize ruleMap
  ruleMap["point"]=rule::point;
  ruleMap["centroid"]=rule::centroid;
  ruleMap["mean"]=rule::mean;
  ruleMap["stdev"]=rule::stdev;
  ruleMap["median"]=rule::median;
  ruleMap["proportion"]=rule::proportion;
  ruleMap["count"]=rule::count;
  ruleMap["min"]=rule::min;
  ruleMap["max"]=rule::max;
  ruleMap["custom"]=rule::custom;
  ruleMap["mode"]=rule::mode;
  ruleMap["sum"]=rule::sum;
  ruleMap["percentile"]=rule::percentile;
  ruleMap["allpoints"]=rule::allpoints;

  if(srcnodata_opt.size()){
    while(srcnodata_opt.size()<bndnodata_opt.size())
      srcnodata_opt.push_back(srcnodata_opt[0]);
  }

  if(verbose_opt[0])
    std::cout << class_opt << std::endl;
  statfactory::StatFactory stat;
  stat.setNoDataValues(srcnodata_opt);
  Vector2d<unsigned int> posdata;
  unsigned long int nsample=0;
  unsigned long int ntotalvalid=0;
  unsigned long int ntotalinvalid=0;

  ImgReaderGdal imgReader;
  if(image_opt.empty()){
    std::cerr << "No image dataset provided (use option -i). Use --help for help information";
      exit(1);
  }
  if(output_opt.empty()){
    std::cerr << "No output dataset provided (use option -o). Use --help for help information";
      exit(1);
  }
  try{
    imgReader.open(image_opt[0]);
  }
  catch(std::string errorstring){
    std::cout << errorstring << std::endl;
    exit(1);
  }

  //check if rule contains allpoints
  if(find(rule_opt.begin(),rule_opt.end(),"allpoints")!=rule_opt.end()){
    //allpoints should be the only rule
    rule_opt.clear();
    rule_opt.push_back("allpoints");
  }

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
  
  int nband=(band_opt.size()) ? band_opt.size() : imgReader.nrOfBand();
  if(class_opt.size()){
    if(nband>1){
      cerr << "Warning: using only first band or multiband image" << endl;
      nband=1;
      band_opt.clear();
      band_opt.push_back(0);
    }
  }

  if(verbose_opt[0]>1)
    std::cout << "Number of bands in input image: " << imgReader.nrOfBand() << std::endl;

  OGRFieldType fieldType;
  int ogr_typecount=11;//hard coded for now!
  if(verbose_opt[0]>1)
    std::cout << "field and label types can be: ";
  for(int iType = 0; iType < ogr_typecount; ++iType){
    if(verbose_opt[0]>1)
      std::cout << " " << OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType);
    if( OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType) != NULL
        && EQUAL(OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType),
                 ftype_opt[0].c_str()))
      fieldType=(OGRFieldType) iType;
  }
  switch( fieldType ){
  case OFTInteger:
  case OFTReal:
  case OFTRealList:
  case OFTString:
    if(verbose_opt[0]>1)
      std::cout << std::endl << "field type is: " << OGRFieldDefn::GetFieldTypeName(fieldType) << std::endl;
    break;
  default:
    cerr << "field type " << OGRFieldDefn::GetFieldTypeName(fieldType) << " not supported" << std::endl;
    exit(1);
    break;
  }
  
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  srand(time(NULL));

  bool sampleIsRaster=false;
  bool sampleIsVirtual=false;

  ImgReaderOgr sampleReaderOgr;
  ImgWriterOgr sampleWriterOgr;

  if(sample_opt.size()){
    try{
      sampleReaderOgr.open(sample_opt[0]);
    }
    catch(string errorString){
      sampleIsRaster=true;
    }
  }
  else{
    try{
      sampleWriterOgr.open("/vsimem/virtual",ogrformat_opt[0]);
      char     **papszOptions=NULL;
      sampleIsVirtual=true;

      if(random_opt.size()){
        //create simple random sampling within boundary
        double ulx,uly,lrx,lry;
        imgReader.getBoundingBox(ulx,uly,lrx,lry);
        if(random_opt[0]>0)
          sampleWriterOgr.createLayer("points", imgReader.getProjection(), wkbPoint, papszOptions);
        OGRPoint pt;
        unsigned int ipoint;
        for(ipoint=0;ipoint<random_opt[0];++ipoint){
          OGRFeature *pointFeature;
          pointFeature=sampleWriterOgr.createFeature();
          double theX=ulx+static_cast<double>(rand())/(RAND_MAX)*(lrx-ulx);
          double theY=uly-static_cast<double>(rand())/(RAND_MAX)*(uly-lry);
          pt.setX(theX);
          pt.setY(theY);
          pointFeature->SetGeometry( &pt ); 
          if(sampleWriterOgr.createFeature(pointFeature) != OGRERR_NONE ){
            string errorString="Failed to create feature in vector dataset";
            throw(errorString);
          }
          // OGRFeature::DestroyFeature(pointFeature);
        }
        if(!ipoint){
          string errorString="Error: no random point created";
          throw(errorString);
        }
      }
      else if(grid_opt.size()){
        //create systematic grid of points 
        double ulx,uly,lrx,lry;
        imgReader.getBoundingBox(ulx,uly,lrx,lry);
        if(uly-grid_opt[0]/2<lry&&ulx+grid_opt[0]/2>lrx){
          string errorString="Error: grid distance too large";
          throw(errorString);
        }
        else if(grid_opt[0]>0)
          sampleWriterOgr.createLayer("points", imgReader.getProjection(), wkbPoint, papszOptions);
        else{
          string errorString="Error: grid distance must be strictly positive number";
          throw(errorString);
        }
        OGRPoint pt;
        unsigned int ipoint=0;
        for(double theY=uly-grid_opt[0]/2;theY>lry;theY-=grid_opt[0]){
          for(double theX=ulx+grid_opt[0]/2;theX<lrx;theX+=grid_opt[0]){
            if(verbose_opt[0]>1)
              cout << "position: " << theX << " " << theY << endl;
            OGRFeature *pointFeature;
            pointFeature=sampleWriterOgr.createFeature();
            pt.setX(theX);
            pt.setY(theY);
            pointFeature->SetGeometry( &pt ); 
            if(sampleWriterOgr.createFeature(pointFeature) != OGRERR_NONE ){
              string errorString="Failed to create feature in vector dataset";
              throw(errorString);
            }
            ++ipoint;
            // OGRFeature::DestroyFeature(pointFeature);
          }
        }
        if(!ipoint){
          string errorString="Error: no points created in grid";
          throw(errorString);
        }
      }
      else{
        std::cerr << "Error: no sample dataset provided (use option -s). Use --help for help information";
        exit(1);
      }
      sampleWriterOgr.close();
      sampleReaderOgr.open("/vsimem/virtual");
    }
    catch(string errorString){
      cerr << errorString << endl;
      exit(1);
    }
  }

  if(sampleIsRaster){
    cerr << "Error: sample must be vector dataset in OGR format";
    exit(1);
  }
  else{//vector dataset
    if(verbose_opt[0]>1)
      std::cout << "creating image sample writer " << output_opt[0] << std::endl;

    ImgWriterOgr ogrWriter;
    double vectords_ulx;
    double vectords_uly;
    double vectords_lrx;
    double vectords_lry;
    bool calculateSpatialStatistics=false;
    try{
      sampleReaderOgr.getExtent(vectords_ulx,vectords_uly,vectords_lrx,vectords_lry);
      bool hasCoverage=((vectords_ulx < imgReader.getLrx())&&(vectords_lrx > imgReader.getUlx())&&(vectords_lry < imgReader.getUly())&&(vectords_uly > imgReader.getLry()));
      if(!hasCoverage){
	ostringstream ess;
	ess << "No coverage in " << image_opt[0] << " for any layer in " << sample_opt[0] << endl;
	throw(ess.str());
      }
      ogrWriter.open(output_opt[0],ogrformat_opt[0]);
      //if class_opt not set, get number of classes from input image for these rules
      for(int irule=0;irule<rule_opt.size();++irule){
        switch(ruleMap[rule_opt[irule]]){
        case(rule::point):
        case(rule::centroid):
        case(rule::allpoints):
          break;
        case(rule::proportion):
        case(rule::count):
        case(rule::custom):
        case(rule::mode):{
          if(class_opt.empty()){
            int theBand=0;
            double minValue=0;
            double maxValue=0;
            if(band_opt.size())
              theBand=band_opt[0];
            imgReader.getMinMax(minValue,maxValue,theBand);
            int nclass=maxValue-minValue+1;
            if(nclass<0&&nclass<256){
              string errorString="Could not automatically define classes, please set class option";
              throw(errorString);
            }
            for(int iclass=minValue;iclass<=maxValue;++iclass)
              class_opt.push_back(iclass);
          }
        }//deliberate fall through: calculate spatial statistics for all non-point like rules
        default:
          calculateSpatialStatistics=true;
          break;
        }
      }
    }
    catch(string errorString){
      cerr << errorString << endl;
      exit(1);
    }
    
    //support multiple layers
    int nlayerRead=sampleReaderOgr.getDataSource()->GetLayerCount();
    int ilayerWrite=0;
    unsigned long int ntotalvalid=0;

    if(verbose_opt[0])
      std::cout << "number of layers: " << nlayerRead << endl;

    for(int ilayer=0;ilayer<nlayerRead;++ilayer){
      OGRLayer *readLayer=sampleReaderOgr.getLayer(ilayer);
      string currentLayername=readLayer->GetName();
      int layerIndex=ilayer;
      if(layer_opt.size()){
	vector<string>::const_iterator it=find(layer_opt.begin(),layer_opt.end(),currentLayername);
	if(it==layer_opt.end())
	  continue;
	else
	  layerIndex=it-layer_opt.begin();
      }
      double layer_ulx;
      double layer_uly;
      double layer_lrx;
      double layer_lry;
      sampleReaderOgr.getExtent(layer_ulx,layer_uly,layer_lrx,layer_lry,ilayer);
      bool hasCoverage=((layer_ulx < imgReader.getLrx())&&(layer_lrx > imgReader.getUlx())&&(layer_lry < imgReader.getUly())&&(layer_uly > imgReader.getLry()));
      if(!hasCoverage)
	continue;


      //read entire block for coverage in memory
      //todo: use different data types
      vector< Vector2d<float> > readValuesReal(nband);
      vector< Vector2d<int> > readValuesInt(nband);

      double layer_uli;
      double layer_ulj;
      double layer_lri;
      double layer_lrj;
      imgReader.geo2image(layer_ulx,layer_uly,layer_uli,layer_ulj);
      imgReader.geo2image(layer_lrx,layer_lry,layer_lri,layer_lrj);

      OGRwkbGeometryType layerGeometry=readLayer->GetLayerDefn()->GetGeomType();

      if(layerGeometry==wkbPoint){
        if(calculateSpatialStatistics){
          if(buffer_opt[0]<1)
            buffer_opt[0]=1;
        }
      }
      
      //extend bounding box with buffer
      if(buffer_opt.size()){
        layer_uli-=buffer_opt[0];
        layer_ulj-=buffer_opt[0];
        layer_lri+=buffer_opt[0];
        layer_lrj+=buffer_opt[0];
      }

      //we already checked there is coverage
      layer_uli=(layer_uli<0)? 0 : static_cast<int>(layer_uli);
      layer_ulj=(layer_ulj<0)? 0 : static_cast<int>(layer_ulj);
      layer_lri=(layer_lri>=imgReader.nrOfCol())? imgReader.nrOfCol()-1 : static_cast<int>(layer_lri);
      layer_lrj=(layer_lrj>=imgReader.nrOfRow())? imgReader.nrOfRow()-1 : static_cast<int>(layer_lrj);

      try{
        for(int iband=0;iband<nband;++iband){
          int theBand=(band_opt.size()) ? band_opt[iband] : iband;
          if(theBand<0){
            string errorString="Error: illegal band (must be positive and starting from 0)";
            throw(errorString);
          }
          if(theBand>=imgReader.nrOfBand()){
            string errorString="Error: illegal band (must be lower than number of bands in input raster dataset)";
            throw(errorString);
          }
          if(verbose_opt[0])
            cout << "reading image band " << theBand << endl;
          switch( fieldType ){
          case OFTInteger:
            imgReader.readDataBlock(readValuesInt[iband],GDT_Int32,layer_uli,layer_lri,layer_ulj,layer_lrj,theBand);
            break;
          case OFTReal:
          default:
            imgReader.readDataBlock(readValuesReal[iband],GDT_Float32,layer_uli,layer_lri,layer_ulj,layer_lrj,theBand);
            break;
          }
        }
      }
      catch(std::string e){
        std::cout << e << std::endl;
        exit(1);
      }

      float theThreshold=(threshold_opt.size()==layer_opt.size())? threshold_opt[layerIndex]: threshold_opt[0];
      cout << "processing layer " << currentLayername << endl;
      
      bool createPolygon=true;
      if(find(rule_opt.begin(),rule_opt.end(),"allpoints")!=rule_opt.end())
        createPolygon=false;

      OGRLayer *writeLayer;
      if(createPolygon){
        //create polygon
	if(verbose_opt[0])
	  std::cout << "create polygons" << std::endl;
	char **papszOptions=NULL;
	writeLayer=ogrWriter.createLayer(readLayer->GetName(), imgReader.getProjection(), wkbPolygon, papszOptions);
      }
      else{
	if(verbose_opt[0])
	  std::cout << "create points in layer " << readLayer->GetName() << std::endl;
	char **papszOptions=NULL;

	writeLayer=ogrWriter.createLayer(readLayer->GetName(), imgReader.getProjection(), wkbPoint, papszOptions);
      }
      if(verbose_opt[0])
	std::cout << "copy fields from layer " << ilayer << std::flush << std::endl;
      ogrWriter.copyFields(sampleReaderOgr,ilayer,ilayerWrite);

      for(int irule=0;irule<rule_opt.size();++irule){
        for(int iband=0;iband<nband;++iband){
          int theBand=(band_opt.size()) ? band_opt[iband] : iband;
          switch(ruleMap[rule_opt[irule]]){
          case(rule::proportion):
          case(rule::count):{//count for each class
            for(int iclass=0;iclass<class_opt.size();++iclass){
              ostringstream fs;
              fs << class_opt[iclass];
              if(nband>1)
                fs << "b" << theBand;
              string fieldname=fs.str();
              ogrWriter.createField(fieldname,fieldType,ilayerWrite);
            }
            break;
          }
          default:{
            ostringstream fs;
            if(rule_opt.size()>1||nband==1)
              fs << rule_opt[irule];
            if(nband>1)
              fs << "b" << theBand;
            string fieldname=fs.str();
            ogrWriter.createField(fieldname,fieldType,ilayerWrite);
          }
          }
        }
      }
      OGRFeature *readFeature;
      unsigned long int ifeature=0;
      unsigned long int nfeatureLayer=sampleReaderOgr.getFeatureCount(ilayer);
      unsigned long int ntotalvalidLayer=0;

      if(nfeatureLayer<=0)
        continue;
      progress=0;
      pfnProgress(progress,pszMessage,pProgressArg);
      readLayer->ResetReading();
      while( (readFeature = readLayer->GetNextFeature()) != NULL ){
	bool validFeature=false;
	if(verbose_opt[0]>2)
	  std::cout << "reading feature " << readFeature->GetFID() << std::endl;
	if(theThreshold>0){//percentual value
	  double p=static_cast<double>(rand())/(RAND_MAX);
	  p*=100.0;
	  if(p>theThreshold){
            continue;//do not select for now, go to next feature
	  }
	}
	else{//absolute value
	  if(threshold_opt.size()==layer_opt.size()){
	    if(ntotalvalidLayer>=-theThreshold){
              continue;//do not select any more pixels, go to next column feature
	    }
	  }
	  else{
	    if(ntotalvalid>=-theThreshold){
              continue;//do not select any more pixels, go to next column feature
	    }
	  }
	}
	if(verbose_opt[0]>2)
	  std::cout << "processing feature " << readFeature->GetFID() << std::endl;
	//get x and y from readFeature
	// double x,y;
	OGRGeometry *poGeometry;
	poGeometry = readFeature->GetGeometryRef();
	assert(poGeometry!=NULL);
	try{
	  if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ){
	    OGRPoint readPoint = *((OGRPoint *) poGeometry);
            
	    double i_centre,j_centre;
            imgReader.geo2image(readPoint.getX(),readPoint.getY(),i_centre,j_centre);
	    //nearest neighbour
	    j_centre=static_cast<int>(j_centre);
	    i_centre=static_cast<int>(i_centre);

	    double uli=i_centre-buffer_opt[0];
	    double ulj=j_centre-buffer_opt[0];
	    double lri=i_centre+buffer_opt[0];
	    double lrj=j_centre+buffer_opt[0];

	    //nearest neighbour
	    ulj=static_cast<int>(ulj);
	    uli=static_cast<int>(uli);
	    lrj=static_cast<int>(lrj);
	    lri=static_cast<int>(lri);

	    //check if j is out of bounds
	    if(static_cast<int>(ulj)<0||static_cast<int>(ulj)>=imgReader.nrOfRow())
	      continue;
	    //check if j is out of bounds
	    if(static_cast<int>(uli)<0||static_cast<int>(lri)>=imgReader.nrOfCol())
	      continue;
            
	    OGRPoint ulPoint,urPoint,llPoint,lrPoint;
	    double ulx;
            double uly;
	    double lrx;
            double lry;

	    OGRPolygon writePolygon;
            OGRPoint writePoint;
	    OGRLinearRing writeRing;
	    OGRFeature *writePolygonFeature;

	    int nPointPolygon=0;
	    if(createPolygon){
	      if(disc_opt[0]){
		double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
		double radius=buffer_opt[0]*sqrt(imgReader.getDeltaX()*imgReader.getDeltaY());
		unsigned short nstep = 25;
		for(int i=0;i<nstep;++i){
		  OGRPoint aPoint;
		  aPoint.setX(readPoint.getX()+imgReader.getDeltaX()/2.0+radius*cos(2*PI*i/nstep));
		  aPoint.setY(readPoint.getY()-imgReader.getDeltaY()/2.0+radius*sin(2*PI*i/nstep));
		  writeRing.addPoint(&aPoint);
		}
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
	      }
	      else{
		double ulx,uly,lrx,lry;
		imgReader.image2geo(uli,ulj,ulx,uly);
		imgReader.image2geo(lri,lrj,lrx,lry);
		ulPoint.setX(ulx-imgReader.getDeltaX()/2.0);
		ulPoint.setY(uly+imgReader.getDeltaY()/2.0);
		lrPoint.setX(lrx+imgReader.getDeltaX()/2.0);
		lrPoint.setY(lry-imgReader.getDeltaY()/2.0);
		urPoint.setX(lrx+imgReader.getDeltaX()/2.0);
		urPoint.setY(uly+imgReader.getDeltaY()/2.0);
		llPoint.setX(ulx-imgReader.getDeltaX()/2.0);
		llPoint.setY(lry-imgReader.getDeltaY()/2.0);

		writeRing.addPoint(&ulPoint);
		writeRing.addPoint(&urPoint);
		writeRing.addPoint(&lrPoint);
		writeRing.addPoint(&llPoint);
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
	      }
              writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
              if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
                cerr << "writing feature failed" << std::endl;
              writePolygonFeature->SetGeometry(&writePolygon);
              if(verbose_opt[0]>1)
                std::cout << "copying new fields write polygon " << std::endl;
              if(verbose_opt[0]>1)
                std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;

              OGRPoint readPoint;
              if(find(rule_opt.begin(),rule_opt.end(),"centroid")!=rule_opt.end()){
                if(verbose_opt[0]>1)
                  std::cout << "get centroid" << std::endl;
                writePolygon.Centroid(&readPoint);
                double i,j;
                imgReader.geo2image(readPoint.getX(),readPoint.getY(),i,j);
                int indexJ=static_cast<int>(j-layer_ulj);
                int indexI=static_cast<int>(i-layer_uli);
                bool valid=true;
                valid=valid&&(indexJ>=0);
                valid=valid&&(indexJ<imgReader.nrOfRow());
                valid=valid&&(indexI>=0);
                valid=valid&&(indexI<imgReader.nrOfCol());
                if(valid&&srcnodata_opt.size()){
                  for(int vband=0;vband<bndnodata_opt.size();++vband){
                    switch( fieldType ){
                    case OFTInteger:{
                      int value;
                      value=((readValuesInt[vband])[indexJ])[indexI];
                      if(value==srcnodata_opt[vband])
                        valid=false;
                      break;
                    }
                    case OFTReal:{
                      double value;
                      value=((readValuesReal[vband])[indexJ])[indexI];
                      if(value==srcnodata_opt[vband])
                        valid=false;
                      break;
                    }
                    }
                    if(!valid)
                      continue;
                    else
                      validFeature=true;
                  }
                }
                if(valid){
                  assert(readValuesReal.size()==nband);
                  for(int iband=0;iband<nband;++iband){
                    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                    //write fields for point on surface and centroid
                    string fieldname;
                    ostringstream fs;
                    if(rule_opt.size()>1||nband==1)
                      fs << "centroid";
                    if(nband>1)
                      fs << "b" << theBand;
                    fieldname=fs.str();
                    switch( fieldType ){
                    case OFTInteger:
                      writePolygonFeature->SetField(fieldname.c_str(),static_cast<int>(((readValuesInt[iband])[indexJ])[indexI]));
                      break;
                    case OFTReal:
                      writePolygonFeature->SetField(fieldname.c_str(),((readValuesReal[iband])[indexJ])[indexI]);
                      break;
                    default://not supported
                      std::string errorString="field type not supported";
                      throw(errorString);
                      break;
                    }
                  }
                }
              }//if centroid
              if(find(rule_opt.begin(),rule_opt.end(),"point")!=rule_opt.end()){
                if(verbose_opt[0]>1)
                  std::cout << "get point on surface" << std::endl;
                if(writePolygon.PointOnSurface(&readPoint)!=OGRERR_NONE)
                  writePolygon.Centroid(&readPoint);
                double i,j;
                imgReader.geo2image(readPoint.getX(),readPoint.getY(),i,j);
                int indexJ=static_cast<int>(j-layer_ulj);
                int indexI=static_cast<int>(i-layer_uli);
                bool valid=true;
                valid=valid&&(indexJ>=0);
                valid=valid&&(indexJ<imgReader.nrOfRow());
                valid=valid&&(indexI>=0);
                valid=valid&&(indexI<imgReader.nrOfCol());
                if(valid&&srcnodata_opt.size()){
                  for(int vband=0;vband<bndnodata_opt.size();++vband){
                    switch( fieldType ){
                    case OFTInteger:{
                      int value;
                      value=((readValuesInt[vband])[indexJ])[indexI];
                      if(value==srcnodata_opt[vband])
                        valid=false;
                      break;
                    }
                    case OFTReal:{
                      double value;
                      value=((readValuesReal[vband])[indexJ])[indexI];
                      if(value==srcnodata_opt[vband])
                        valid=false;
                      break;
                    }
                    }
                    if(!valid)
                      continue;
                    else
                      validFeature=true;
                  }
                }
                if(valid){
                  for(int iband=0;iband<nband;++iband){
                    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                    //write fields for point on surface and centroid
                    string fieldname;
                    ostringstream fs;
                    if(rule_opt.size()>1||nband==1)
                      fs << "point";
                    if(nband>1)
                      fs << "b" << theBand;
                    fieldname=fs.str();
                    switch( fieldType ){
                    case OFTInteger:
                      writePolygonFeature->SetField(fieldname.c_str(),static_cast<int>(((readValuesInt[iband])[indexJ])[indexI]));
                      break;
                    case OFTReal:
                      writePolygonFeature->SetField(fieldname.c_str(),((readValuesReal[iband])[indexJ])[indexI]);
                      break;
                    default://not supported
                      std::string errorString="field type not supported";
                      throw(errorString);
                      break;
                    }
                  }
                }
              }//if point
            }//if createPolygon

            if(calculateSpatialStatistics||!createPolygon){
              Vector2d<double> polyValues;
              vector<double> polyClassValues;
	    
              if(class_opt.size()){
                polyClassValues.resize(class_opt.size());
                //initialize
                for(int iclass=0;iclass<class_opt.size();++iclass)
                  polyClassValues[iclass]=0;
              }
              else
                polyValues.resize(nband);

              OGRPoint thePoint;
              for(int j=ulj;j<=lrj;++j){
                for(int i=uli;i<=lri;++i){
                  //check if within raster image
                  if(i<0||i>=imgReader.nrOfCol())
                    continue;
                  if(j<0||j>=imgReader.nrOfRow())
                    continue;
                  int indexJ=j-layer_ulj;
                  int indexI=i-layer_uli;
                  if(indexJ<0)
                    indexJ=0;
                  if(indexI<0)
                    indexI=0;
                  if(indexJ>=imgReader.nrOfRow())
                    indexJ=imgReader.nrOfRow()-1;
                  if(indexI>=imgReader.nrOfCol())
                    indexI=imgReader.nrOfCol()-1;

                  double theX=0;
                  double theY=0;
                  imgReader.image2geo(i,j,theX,theY);
                  thePoint.setX(theX);
                  thePoint.setY(theY);
                  if(disc_opt[0]&&buffer_opt[0]>0){
                    double radius=buffer_opt[0]*sqrt(imgReader.getDeltaX()*imgReader.getDeltaY());
                    if((theX-readPoint.getX())*(theX-readPoint.getX())+(theY-readPoint.getY())*(theY-readPoint.getY())>radius*radius)
                      continue;
                  }
                  bool valid=true;

                  if(srcnodata_opt.size()){
                    for(int vband=0;vband<bndnodata_opt.size();++vband){
                      switch( fieldType ){
                      case OFTInteger:{
                        int value=((readValuesInt[vband])[indexJ])[indexI];
                        if(value==srcnodata_opt[vband]){
                          valid=false;
                        }
                        break;
                      }
                      default:{
                        float value=((readValuesReal[vband])[indexJ])[indexI];
                        if(value==srcnodata_opt[vband]){
                          valid=false;
                        }
                        break;
                      }
                      }
                    }
                  }
                  if(!valid)
                    continue;
                  else
                    validFeature=true;

                  ++nPointPolygon;
                  OGRFeature *writePointFeature;
                  if(valid&&!createPolygon){//write all points
                    if(polythreshold_opt.size())
                      if(nPointPolygon>polythreshold_opt[0])
                        break;
                    //create feature
                    writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
                    if(verbose_opt[0]>1)
                      std::cout << "copying fields from point feature " << std::endl;
                    if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
                      cerr << "writing feature failed" << std::endl;
                    if(verbose_opt[0]>1)
                      std::cout << "set geometry as point " << std::endl;
                    writePointFeature->SetGeometry(&thePoint);
                    assert(wkbFlatten(writePointFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
                    if(verbose_opt[0]>1){
                      std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields:" << std::endl;
                      for(int iField=0;iField<writePointFeature->GetFieldCount();++iField){
                        std::string fieldname=writeLayer->GetLayerDefn()->GetFieldDefn(iField)->GetNameRef();
                        cout << fieldname << endl;
                      }
                    }
                  }
                  if(valid&&class_opt.size()){
                    short value=0;
                    switch( fieldType ){
                    case OFTInteger:
                      value=((readValuesInt[0])[indexJ])[indexI];
                      break;
                    case OFTReal:
                      value=((readValuesReal[0])[indexJ])[indexI];
                      break;
                    }
                    for(int iclass=0;iclass<class_opt.size();++iclass){
                      if(value==class_opt[iclass])
                        polyClassValues[iclass]+=1;
                    }
                  }
                  else if(valid){
                    for(int iband=0;iband<nband;++iband){
                      int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                      double value=0;
                      switch( fieldType ){
                      case OFTInteger:
                        value=((readValuesInt[iband])[indexJ])[indexI];
                        break;
                      case OFTReal:
                        value=((readValuesReal[iband])[indexJ])[indexI];
                        break;
                      }

                      if(verbose_opt[0]>1)
                        std::cout << ": " << value << std::endl;
                      if(!createPolygon){//write all points within polygon
                        string fieldname;
                        ostringstream fs;
                        if(rule_opt.size()>1||nband==1)
                          fs << "allpoints";
                        if(nband>1)
                          fs << "b" << theBand;
                        fieldname=fs.str();
                        int fieldIndex=writePointFeature->GetFieldIndex(fieldname.c_str());
                        if(fieldIndex<0){
                          cerr << "field " << fieldname << " was not found" << endl;
                          exit(1);
                        }
                        if(verbose_opt[0]>1)
                          std::cout << "set field " << fieldname << " to " << value << std::endl;
                        switch( fieldType ){
                        case OFTInteger:
                        case OFTReal:
                          writePointFeature->SetField(fieldname.c_str(),value);
                          break;
                        default://not supported
                          assert(0);
                          break;
                        }
                      }
                      else{
                        polyValues[iband].push_back(value);
                      }
                    }//iband
                  }//else (not class_opt.size())
                  if(valid&&!createPolygon){
                    //write feature
                    if(verbose_opt[0]>1)
                      std::cout << "creating point feature" << std::endl;
                    if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
                      std::string errorString="Failed to create feature in ogr vector dataset";
                      throw(errorString);
                    }
                    //destroy feature
                    // OGRFeature::DestroyFeature( writePointFeature );
                    ++ntotalvalid;
                    ++ntotalvalidLayer;
                  }
                }//for in i
              }//for int j

              if(createPolygon){
                //do not create if no points found within polygon
                if(!nPointPolygon){
                  if(verbose_opt[0])
                    cout << "no points found in polygon, continuing" << endl;
                  continue;
                }
                //write field attributes to polygon feature
                for(int irule=0;irule<rule_opt.size();++irule){
                  //skip centroid and point
                  if(ruleMap[rule_opt[irule]]==rule::centroid||ruleMap[rule_opt[irule]]==rule::point)
                    continue;
                  for(int iband=0;iband<nband;++iband){
                    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                    double theValue=0;
                    string fieldname;
                    ostringstream fs;
                    if(rule_opt.size()>1||nband==1)
                      fs << rule_opt[irule];
                    if(nband>1)
                      fs << "b" << theBand;
                    fieldname=fs.str();
                    switch(ruleMap[rule_opt[irule]]){
                    case(rule::proportion):
                      stat.normalize_pct(polyClassValues);
                    case(rule::count):{//count for each class
                      for(int index=0;index<polyClassValues.size();++index){
                        theValue=polyClassValues[index];
                        ostringstream fs;
                        fs << class_opt[index];
                        if(nband>1)
                          fs << "b" << theBand;
                        fieldname=fs.str();
                      }
                      break;
                    }
                    case(rule::mode):{
                      //maximum votes in polygon
                      if(verbose_opt[0])
                        std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
                      //search for class with maximum votes
                      int maxClass=stat.mymin(class_opt);
                      vector<double>::iterator maxit;
                      maxit=stat.mymax(polyClassValues,polyClassValues.begin(),polyClassValues.end());
                      int maxIndex=distance(polyClassValues.begin(),maxit);
                      maxClass=class_opt[maxIndex];
                      if(verbose_opt[0]>0)
                        std::cout << "maxClass: " << maxClass << std::endl;
                      theValue=maxClass;
                    }
                    case(rule::mean):
                      theValue=stat.mean(polyValues[iband]);
                      break;
                    case(rule::median):
                      theValue=stat.median(polyValues[iband]);
                      break;
                    case(rule::stdev):
                      theValue=sqrt(stat.var(polyValues[iband]));
                      break;
                    case(rule::percentile):
                      theValue=stat.percentile(polyValues[iband],polyValues[iband].begin(),polyValues[iband].end(),percentile_opt[0]);
                      break;
                    case(rule::sum):
                      theValue=stat.sum(polyValues[iband]);
                      break;
                    case(rule::max):
                      theValue=stat.mymax(polyValues[iband]);
                      break;
                    case(rule::min):
                      theValue=stat.mymin(polyValues[iband]);
                      break;
                    case(rule::centroid):
                      theValue=polyValues[iband].back();
                      break;
                    default://not supported
                      break;
                    }
                  
                    switch( fieldType ){
                    case OFTInteger:
                      writePolygonFeature->SetField(fieldname.c_str(),static_cast<int>(theValue));
                      break;
                    case OFTReal:
                      writePolygonFeature->SetField(fieldname.c_str(),theValue);
                      break;
                    case OFTString:
                      writePolygonFeature->SetField(fieldname.c_str(),type2string<double>(theValue).c_str());
                      break;
                    default://not supported
                      std::string errorString="field type not supported";
                      throw(errorString);
                      break;
                    }
                  }
                }
              }
            }
            if(createPolygon){
              //write polygon feature
              if(verbose_opt[0]>1)
                std::cout << "creating polygon feature" << std::endl;
              if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
                std::string errorString="Failed to create polygon feature in ogr vector dataset";
                throw(errorString);
              }
              //test: no need to destroy anymore?
              // OGRFeature::DestroyFeature( writePolygonFeature );
              ++ntotalvalid;
              ++ntotalvalidLayer;
            }
          }
	  else{
	    OGRPolygon readPolygon;
	    OGRMultiPolygon readMultiPolygon;

            //get envelope
            OGREnvelope* psEnvelope=new OGREnvelope();

	    if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
	      readPolygon = *((OGRPolygon *) poGeometry);
	      readPolygon.closeRings();
	      readPolygon.getEnvelope(psEnvelope);
	    }
	    else if(wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon){
	      readMultiPolygon = *((OGRMultiPolygon *) poGeometry);
	      readMultiPolygon.closeRings();
	      readMultiPolygon.getEnvelope(psEnvelope);
	    }
	    else{
	      std::string test;
	      test=poGeometry->getGeometryName();
	      ostringstream oss;
	      oss << "geometry " << test << " not supported";
	      throw(oss.str());
	    }

	    double ulx,uly,lrx,lry;
	    double uli,ulj,lri,lrj;
            ulx=psEnvelope->MinX;
            uly=psEnvelope->MaxY;
            lrx=psEnvelope->MaxX;
            lry=psEnvelope->MinY;
            delete psEnvelope;

            //check if feature is covered by input raster dataset
            if(!imgReader.covers(ulx,uly,lrx,lry))
              continue;

	    OGRFeature *writePolygonFeature;
	    int nPointPolygon=0;
	    if(createPolygon){
              writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
              writePolygonFeature->SetGeometry(poGeometry);
              //writePolygonFeature and readFeature are both of type wkbPolygon
              if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
                cerr << "writing feature failed" << std::endl;
              if(verbose_opt[0]>1)
                std::cout << "copying new fields write polygon " << std::endl;
              if(verbose_opt[0]>1)
                std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
	    }

            OGRPoint readPoint;
	    if(find(rule_opt.begin(),rule_opt.end(),"centroid")!=rule_opt.end()){
              if(verbose_opt[0]>1)
                std::cout << "get centroid" << std::endl;
	      if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon)
		readPolygon.Centroid(&readPoint);
	      else if(wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon)
		readMultiPolygon.Centroid(&readPoint);

              double i,j;
              imgReader.geo2image(readPoint.getX(),readPoint.getY(),i,j);
              int indexJ=static_cast<int>(j-layer_ulj);
              int indexI=static_cast<int>(i-layer_uli);
              bool valid=true;
              valid=valid&&(indexJ>=0);
              valid=valid&&(indexJ<imgReader.nrOfRow());
              valid=valid&&(indexI>=0);
              valid=valid&&(indexI<imgReader.nrOfCol());
              if(valid&&srcnodata_opt.size()){
                for(int vband=0;vband<bndnodata_opt.size();++vband){
                  switch( fieldType ){
                  case OFTInteger:{
                    int value;
                    value=((readValuesInt[vband])[indexJ])[indexI];
                    if(value==srcnodata_opt[vband])
                      valid=false;
                    break;
                  }
                  case OFTReal:{
                    double value;
                    value=((readValuesReal[vband])[indexJ])[indexI];
                    if(value==srcnodata_opt[vband])
                      valid=false;
                    break;
                  }
                  }
                  if(!valid)
                    continue;
                  else
                    validFeature=true;
                }
              }
              if(valid){
                assert(readValuesReal.size()==nband);
                for(int iband=0;iband<nband;++iband){
                  int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                  //write fields for point on surface and centroid
                  string fieldname;
                  ostringstream fs;
                  if(rule_opt.size()>1||nband==1)
                    fs << "centroid";
                  if(nband>1)
                    fs << "b" << theBand;
                  fieldname=fs.str();
                  switch( fieldType ){
                  case OFTInteger:
                    writePolygonFeature->SetField(fieldname.c_str(),static_cast<int>(((readValuesInt[iband])[indexJ])[indexI]));
                    break;
                  case OFTReal:
                    writePolygonFeature->SetField(fieldname.c_str(),((readValuesReal[iband])[indexJ])[indexI]);
                    break;
                  default://not supported
                    std::string errorString="field type not supported";
                    throw(errorString);
                    break;
                  }
                }
              }
            }
	    if(find(rule_opt.begin(),rule_opt.end(),"point")!=rule_opt.end()){
              if(verbose_opt[0]>1)
                std::cout << "get point on surface" << std::endl;
	      if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
		if(readPolygon.PointOnSurface(&readPoint)!=OGRERR_NONE)
		  readPolygon.Centroid(&readPoint);
	      }
	      else if(wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon){
		// if(readMultiPolygon.PointOnSurface(&readPoint)!=OGRERR_NONE)
		  readMultiPolygon.Centroid(&readPoint);
	      }
              double i,j;
              imgReader.geo2image(readPoint.getX(),readPoint.getY(),i,j);
              int indexJ=static_cast<int>(j-layer_ulj);
              int indexI=static_cast<int>(i-layer_uli);
              bool valid=true;
              valid=valid&&(indexJ>=0);
              valid=valid&&(indexJ<imgReader.nrOfRow());
              valid=valid&&(indexI>=0);
              valid=valid&&(indexI<imgReader.nrOfCol());
              if(valid&&srcnodata_opt.size()){
                for(int vband=0;vband<bndnodata_opt.size();++vband){
                  switch( fieldType ){
                  case OFTInteger:{
                    int value;
                    value=((readValuesInt[vband])[indexJ])[indexI];
                    if(value==srcnodata_opt[vband])
                      valid=false;
                    break;
                  }
                  case OFTReal:{
                    double value;
                    value=((readValuesReal[vband])[indexJ])[indexI];
                    if(value==srcnodata_opt[vband])
                      valid=false;
                    break;
                  }
                  }
                  if(!valid)
                    continue;
                  else
                    validFeature=true;
                }
              }
              if(valid){
                for(int iband=0;iband<nband;++iband){
                  int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                  //write fields for point on surface and centroid
                  string fieldname;
                  ostringstream fs;
                  if(rule_opt.size()>1||nband==1)
                    fs << "point";
                  if(nband>1)
                    fs << "b" << theBand;
                  fieldname=fs.str();
                  switch( fieldType ){
                  case OFTInteger:
                    writePolygonFeature->SetField(fieldname.c_str(),static_cast<int>(((readValuesInt[iband])[indexJ])[indexI]));
                    break;
                  case OFTReal:
                    writePolygonFeature->SetField(fieldname.c_str(),((readValuesReal[iband])[indexJ])[indexI]);
                    break;
                  default://not supported
                    std::string errorString="field type not supported";
                    throw(errorString);
                    break;
                  }
                }
              }
            }
            if(calculateSpatialStatistics||ruleMap[rule_opt[0]]==rule::allpoints){
              imgReader.geo2image(ulx,uly,uli,ulj);
              imgReader.geo2image(lrx,lry,lri,lrj);
              //nearest neighbour
              ulj=static_cast<int>(ulj);
              uli=static_cast<int>(uli);
              lrj=static_cast<int>(lrj);
              lri=static_cast<int>(lri);
              //iterate through all pixels
              if(verbose_opt[0]>1)
                std::cout << "bounding box for polygon feature " << ifeature << ": " << uli << " " << ulj << " " << lri << " " << lrj << std::endl;

              if(uli<0)
                uli=0;
              if(lri<0)
                lri=0;
              if(uli>=imgReader.nrOfCol())
                uli=imgReader.nrOfCol()-1;
              if(lri>=imgReader.nrOfCol())
                lri=imgReader.nrOfCol()-1;
              if(ulj<0)
                ulj=0;
              if(lrj<0)
                lrj=0;
              if(ulj>=imgReader.nrOfRow())
                ulj=imgReader.nrOfRow()-1;
              if(lrj>=imgReader.nrOfRow())
                lrj=imgReader.nrOfRow()-1;

              Vector2d<double> polyValues;
              vector<double> polyClassValues;
	    
              if(class_opt.size()){
                polyClassValues.resize(class_opt.size());
                //initialize
                for(int iclass=0;iclass<class_opt.size();++iclass)
                  polyClassValues[iclass]=0;
              }
              else
                polyValues.resize(nband);

              OGRPoint thePoint;
              for(int j=ulj;j<=lrj;++j){
                for(int i=uli;i<=lri;++i){
                  //check if within raster image
                  if(i<0||i>=imgReader.nrOfCol())
                    continue;
                  if(j<0||j>=imgReader.nrOfRow())
                    continue;
                  int indexJ=j-layer_ulj;
                  int indexI=i-layer_uli;

                  double theX=0;
                  double theY=0;
                  imgReader.image2geo(i,j,theX,theY);
                  thePoint.setX(theX);
                  thePoint.setY(theY);
                  //check if point is on surface
		  if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
		    if(!readPolygon.Contains(&thePoint))
		      continue;
		  }
		  else if(wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon){
		    if(!readMultiPolygon.Contains(&thePoint))
		      continue;
		  }
		  
                  bool valid=true;
                  if(srcnodata_opt.size()){
                    for(int vband=0;vband<bndnodata_opt.size();++vband){
                      switch( fieldType ){
                      case OFTInteger:{
                        int value=((readValuesInt[vband])[indexJ])[indexI];
                        if(value==srcnodata_opt[vband]){
                          valid=false;
                        }
                        break;
                      }
                      default:{
                        float value=((readValuesReal[vband])[indexJ])[indexI];
                        if(value==srcnodata_opt[vband]){
                          valid=false;
                        }
                        break;
                      }
                      }
                    }
                  }
                  if(!valid)
                    continue;
                  else
                    validFeature=true;

                  if(verbose_opt[0]>1)
                    std::cout << "point is on surface:" << thePoint.getX() << "," << thePoint.getY() << std::endl;
                  ++nPointPolygon;

                  OGRFeature *writePointFeature;
                  if(!createPolygon){//write all points within polygon
                    if(polythreshold_opt.size())
                      if(nPointPolygon>polythreshold_opt[0])
                        break;
                    //create feature
                    writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
                    if(verbose_opt[0]>1)
                      std::cout << "copying fields from polygons " << std::endl;
                    if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
                      cerr << "writing feature failed" << std::endl;
                    if(verbose_opt[0]>1)
                      std::cout << "set geometry as point " << std::endl;
                    writePointFeature->SetGeometry(&thePoint);
                    assert(wkbFlatten(writePointFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
                    if(verbose_opt[0]>1){
                      std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields:" << std::endl;
                      for(int iField=0;iField<writePointFeature->GetFieldCount();++iField){
                        std::string fieldname=writeLayer->GetLayerDefn()->GetFieldDefn(iField)->GetNameRef();
                        cout << fieldname << endl;
                      }
                    }
                  }
                  if(class_opt.size()){
                    short value=0;
                    switch( fieldType ){
                    case OFTInteger:
                      value=((readValuesInt[0])[indexJ])[indexI];
                      break;
                    case OFTReal:
                      value=((readValuesReal[0])[indexJ])[indexI];
                      break;
                    }
                    for(int iclass=0;iclass<class_opt.size();++iclass){
                      if(value==class_opt[iclass])
                        polyClassValues[iclass]+=1;
                    }
                  }
                  else{
                    for(int iband=0;iband<nband;++iband){
                      int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                      double value=0;
                      switch( fieldType ){
                      case OFTInteger:
                        value=((readValuesInt[iband])[indexJ])[indexI];
                        break;
                      case OFTReal:
                        value=((readValuesReal[iband])[indexJ])[indexI];
                        break;
                      }

                      if(verbose_opt[0]>1)
                        std::cout << ": " << value << std::endl;
                      if(!createPolygon){//write all points within polygon
                        string fieldname;
                        ostringstream fs;
                        if(rule_opt.size()>1||nband==1)
                          fs << "allpoints";
                        if(nband>1)
                          fs << "b" << theBand;
                        fieldname=fs.str();
                        int fieldIndex=writePointFeature->GetFieldIndex(fieldname.c_str());
                        if(fieldIndex<0){
                          cerr << "field " << fieldname << " was not found" << endl;
                          exit(1);
                        }
                        if(verbose_opt[0]>1)
                          std::cout << "set field " << fieldname << " to " << value << std::endl;
                        switch( fieldType ){
                        case OFTInteger:
                        case OFTReal:
                          writePointFeature->SetField(fieldname.c_str(),value);
                          break;
                        default://not supported
                          assert(0);
                          break;
                        }
                      }
                      else{
                        polyValues[iband].push_back(value);
                      }
                    }//iband
                  }//else (not class_opt.size())
                  if(!createPolygon){
                    //write feature
                    if(verbose_opt[0]>1)
                      std::cout << "creating point feature" << std::endl;
                    if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
                      std::string errorString="Failed to create feature in ogr vector dataset";
                      throw(errorString);
                    }
                    //destroy feature
                    // OGRFeature::DestroyFeature( writePointFeature );
                    ++ntotalvalid;
                    ++ntotalvalidLayer;
                  }
                }//for in i
              }//for int j
              if(createPolygon){
                //do not create if no points found within polygon
                if(!nPointPolygon){
                  if(verbose_opt[0])
                    cout << "no points found in polygon, continuing" << endl;
                  continue;
                }
                //write field attributes to polygon feature
                for(int irule=0;irule<rule_opt.size();++irule){
                  //skip centroid and point
                  if(ruleMap[rule_opt[irule]]==rule::centroid||ruleMap[rule_opt[irule]]==rule::point)
                    continue;
                  for(int iband=0;iband<nband;++iband){
                    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                    double theValue=0;
                    string fieldname;
                    ostringstream fs;
                    if(rule_opt.size()>1||nband==1)
                      fs << rule_opt[irule];
                    if(nband>1)
                      fs << "b" << theBand;
                    fieldname=fs.str();
                    switch(ruleMap[rule_opt[irule]]){
                    case(rule::proportion):
                      stat.normalize_pct(polyClassValues);
                    case(rule::count):{//count for each class
                      for(int index=0;index<polyClassValues.size();++index){
                        theValue=polyClassValues[index];
                        ostringstream fs;
                        fs << class_opt[index];
                        if(nband>1)
                          fs << "b" << theBand;
                        fieldname=fs.str();
                      }
                      break;
                    }
                    case(rule::mode):{
                      //maximum votes in polygon
                      if(verbose_opt[0])
                        std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
                      //search for class with maximum votes
                      int maxClass=stat.mymin(class_opt);
                      vector<double>::iterator maxit;
                      maxit=stat.mymax(polyClassValues,polyClassValues.begin(),polyClassValues.end());
                      int maxIndex=distance(polyClassValues.begin(),maxit);
                      maxClass=class_opt[maxIndex];
                      if(verbose_opt[0]>0)
                        std::cout << "maxClass: " << maxClass << std::endl;
                      theValue=maxClass;
                    }
                    case(rule::mean):
                      theValue=stat.mean(polyValues[iband]);
                      break;
                    case(rule::median):
                      theValue=stat.median(polyValues[iband]);
                      break;
                    case(rule::stdev):
                      theValue=sqrt(stat.var(polyValues[iband]));
                      break;
                    case(rule::percentile):
                      theValue=stat.percentile(polyValues[iband],polyValues[iband].begin(),polyValues[iband].end(),percentile_opt[0]);
                      break;
                    case(rule::sum):
                      theValue=stat.sum(polyValues[iband]);
                      break;
                    case(rule::max):
                      theValue=stat.mymax(polyValues[iband]);
                      break;
                    case(rule::min):
                      theValue=stat.mymin(polyValues[iband]);
                      break;
                    case(rule::centroid):
                      theValue=polyValues[iband].back();
                      break;
                    default://not supported
                      break;
                    }
                  
                    switch( fieldType ){
                    case OFTInteger:
                      writePolygonFeature->SetField(fieldname.c_str(),static_cast<int>(theValue));
                      break;
                    case OFTReal:
                      writePolygonFeature->SetField(fieldname.c_str(),theValue);
                      break;
                    case OFTString:
                      writePolygonFeature->SetField(fieldname.c_str(),type2string<double>(theValue).c_str());
                      break;
                    default://not supported
                      std::string errorString="field type not supported";
                      throw(errorString);
                      break;
                    }
                  }
                }
              }
            }
            if(createPolygon){
              //write polygon feature
              if(verbose_opt[0]>1)
                std::cout << "creating polygon feature" << std::endl;
              if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
                std::string errorString="Failed to create polygon feature in ogr vector dataset";
                throw(errorString);
              }
              //test: no need to destroy anymore?
              // OGRFeature::DestroyFeature( writePolygonFeature );
              ++ntotalvalid;
              ++ntotalvalidLayer;
            }
	  }
	  ++ifeature;
	  if(theThreshold>0){
	    if(threshold_opt.size()==layer_opt.size())
	      progress=(100.0/theThreshold)*static_cast<float>(ntotalvalidLayer)/nfeatureLayer;
	    else
	      progress=static_cast<float>(ntotalvalidLayer)/nfeatureLayer;
	  }
 	  else
	    progress=static_cast<float>(ifeature+1)/(-theThreshold);
	  pfnProgress(progress,pszMessage,pProgressArg);
        }
	catch(std::string e){
	  std::cout << e << std::endl;
	  continue;
        }
	catch(int npoint){
	  if(verbose_opt[0])
	    std::cout << "number of points read in polygon: " << npoint << std::endl;
	  continue;
        }
      }
      // if(rbox_opt[0]>0||cbox_opt[0]>0)
      //   boxWriter.close();
      progress=1.0;
      pfnProgress(progress,pszMessage,pProgressArg);
      ++ilayerWrite;
    }//for ilayer
    sampleReaderOgr.close();
    ogrWriter.close();
  }
  progress=1.0;
  pfnProgress(progress,pszMessage,pProgressArg);
  imgReader.close();
}

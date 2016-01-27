/**********************************************************************
pkextract.cc: extract pixel values from raster image from a (vector or raster) sample
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
/*! \page pkextract pkextract
 extract pixel values from raster image from a (vector or raster) sample
## SYNOPSIS

<code>
  Usage: pkextract -i input [-s sample | -rand number | -grid size] -o output
</code>

<code>

  Options: [-ln layer]* [-c class]* [-t threshold]* [-f format] [-ft fieldType] [-lt labelType] [-polygon] [-b band]* [-r rule]

  Advanced options:
  [-sband band -eband band]* [-bndnodata band -srcnodata value]* [-tp threshold] [-test testSample] [-bn attribute] [-cn attribute] [-geo value] [-down value] [-buf value [-circ]]
</code>

\section pkextract_description Description

The utility pkextract extracts pixel values from an input raster dataset, based on the locations you provide via a sample file. Alternatively, a random sample or systematic grid of points can also be extracted. The sample can be a vector file with points or polygons. In the case of polygons, you can either extract the values for all raster pixels that are covered by the polygons, or extract a single value for each polygon such as the centroid, mean, median, etc. As output, a new copy of the vector file is created with an extra attribute for the extracted pixel value. For each raster band in the input image, a separate attribute is created. For instance, if the raster dataset contains three bands, three attributes are created (b0, b1 and b2). 

Instead of a vector dataset, the sample can also be a raster dataset with categorical values. The typical use case is a land cover map that overlaps the input raster dataset. The utility then extracts pixels from the input raster for the respective land cover classes. To select a random subset of the sample raster dataset you can set the threshold option -t with a percentage value. 

A typical usage of pkextract is to prepare a training sample for one of the classifiers implemented in pktools.

\anchor pkextract_rules 

Overview of the possible extraction rules:

extraction rule | output features
--------------- | ---------------
point | Extract all pixel values covered by the polygon (option -polygon not set) or extract a pixel on the surface option (-polygon set).
centroid | Extract pixel value at the centroid of the polygon.
mean | Extract average of all pixel values within the polygon.
stdev | Extract standard deviation of all pixel values within the polygon.
median | Extract median of all pixel values within the polygon.
min | Extract minimum value of all pixels within the polygon.
max | Extract maximum value of all pixels within the polygon.
sum | Extract sum of the values of all pixels within the polygon.
mode | Extract the mode of classes within the polygon (classes must be set with the option class).
proportion | Extract proportion of class(es) within the polygon (classes must be set with the option class).
count | Extract count of class(es) within the polygon (classes must be set with the option class).
percentile | Extract percentile as defined by option perc (e.g, 95th percentile of values covered by polygon).

\section pkextract_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Raster input dataset containing band information | 
 | s      | sample               | std::string |       |OGR vector dataset with features to be extracted from input data. Output will contain features with input band information included. Sample image can also be GDAL raster dataset. | 
 | ln     | ln                   | std::string |       |Layer name(s) in sample (leave empty to select all) | 
 | rand   | random               | unsigned int |       |Create simple random sample of points. Provide number of points to generate | 
 | grid   | grid                 | double |       |Create systematic grid of points. Provide cell grid size (in projected units, e.g,. m) | 
 | o      | output               | std::string |       |Output sample dataset | 
 | c      | class                | int  |       |Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset. Make sure to set classes if rule is set to mode, proportion or count | 
 | t      | threshold            | float | 100   |Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold per vector sample layer. If using raster land cover maps as a sample dataset, you can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es) | 
 | perc   | perc                 | double | 95    |Percentile value used for rule percentile | 
 | f      | f                    | std::string | SQLite |Output sample dataset format | 
 | ft     | ftype                | std::string | Real  |Field type (only Real or Integer) | 
 | lt     | ltype                | std::string | Integer |Label type: In16 or String | 
 | polygon | polygon              | bool | false |Create OGRPolygon as geometry instead of OGRPoint. | 
 | b      | band                 | int  |       |Band index(es) to extract (0 based). Leave empty to use all bands | 
 | sband  | startband            | unsigned short |      |Start band sequence number | 
 | eband  | endband              | unsigned short |      |End band sequence number   | 
 | r      | rule                 | std::string | centroid |Rule how to report image information per feature (only for vector sample). point (value at each point or at centroid if polygon), centroid, mean, stdev, median, proportion, count, min, max, mode, sum, percentile. | 
 | bndnodata | bndnodata            | int  | 0     |Band(s) in input image to check if pixel is valid (used for srcnodata) | 
 | srcnodata | srcnodata            | double |       |Invalid value(s) for input image | 
 | tp     | thresholdPolygon     | float |       |(absolute) threshold for selecting samples in each polygon | 
 | test   | test                 | std::string |       |Test sample dataset (use this option in combination with threshold<100 to create a training (output) and test set | 
 | bn     | bname                | std::string | b     |For single band input data, this extra attribute name will correspond to the raster values. For multi-band input data, multiple attributes with this prefix will be added (e.g. b0, b1, b2, etc.) | 
 | cn     | cname                | std::string | label |Name of the class label in the output vector dataset | 
 | geo    | geo                  | short | 1     |Use geo coordinates (set to 0 to use image coordinates) | 
 | down   | down                 | short | 1     |Down sampling factor (for raster sample datasets only). Can be used to create grid points | 
 | buf    | buffer               | short |       |Buffer for calculating statistics for point features  | 
 | circ   | circular             | bool | false |Use a circular disc kernel buffer (for vector point sample datasets only, use in combination with buffer option) | 

Usage: pkextract -i input [-s sample | -rand number | -grid size] -o output


Examples
========
Some examples how to use pkextract can be found \ref examples_pkextract "here"
**/

namespace rule{
  enum RULE_TYPE {point=0, mean=1, proportion=2, custom=3, min=4, max=5, mode=6, centroid=7, sum=8, median=9, stdev=10, percentile=11, count=12};
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
  Optionpk<string> ltype_opt("lt", "ltype", "Label type: In16 or String", "Integer");
  Optionpk<bool> polygon_opt("polygon", "polygon", "Create OGRPolygon as geometry instead of OGRPoint.", false);
  Optionpk<int> band_opt("b", "band", "Band index(es) to extract (0 based). Leave empty to use all bands");
  Optionpk<unsigned short> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned short> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<string> rule_opt("r", "rule", "Rule how to report image information per feature (only for vector sample). point (value at each point or at centroid if polygon), centroid, mean, stdev, median, proportion, count, min, max, mode, sum, percentile.", "centroid");
  Optionpk<double> srcnodata_opt("srcnodata", "srcnodata", "Invalid value(s) for input image");
  Optionpk<int> bndnodata_opt("bndnodata", "bndnodata", "Band(s) in input image to check if pixel is valid (used for srcnodata)", 0);
  Optionpk<float> polythreshold_opt("tp", "thresholdPolygon", "(absolute) threshold for selecting samples in each polygon");
  Optionpk<string> test_opt("test", "test", "Test sample dataset (use this option in combination with threshold<100 to create a training (output) and test set");
  Optionpk<string> fieldname_opt("bn", "bname", "For single band input data, this extra attribute name will correspond to the raster values. For multi-band input data, multiple attributes with this prefix will be added (e.g. b0, b1, b2, etc.)", "b");
  Optionpk<string> label_opt("cn", "cname", "Name of the class label in the output vector dataset", "label");
  Optionpk<short> geo_opt("geo", "geo", "Use geo coordinates (set to 0 to use image coordinates)", 1);
  Optionpk<short> down_opt("down", "down", "Down sampling factor (for raster sample datasets only). Can be used to create grid points", 1);
  Optionpk<short> buffer_opt("buf", "buffer", "Buffer for calculating statistics for point features ");
  Optionpk<bool> disc_opt("circ", "circular", "Use a circular disc kernel buffer (for vector point sample datasets only, use in combination with buffer option)", false);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose mode if > 0", 0,2);

  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  bndnodata_opt.setHide(1);
  srcnodata_opt.setHide(1);
  polythreshold_opt.setHide(1);
  percentile_opt.setHide(1);
  test_opt.setHide(1);
  fieldname_opt.setHide(1);
  label_opt.setHide(1);
  geo_opt.setHide(1);
  down_opt.setHide(1);
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
    ltype_opt.retrieveOption(argc,argv);
    polygon_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    bstart_opt.retrieveOption(argc,argv);
    bend_opt.retrieveOption(argc,argv);
    rule_opt.retrieveOption(argc,argv);
    bndnodata_opt.retrieveOption(argc,argv);
    srcnodata_opt.retrieveOption(argc,argv);
    polythreshold_opt.retrieveOption(argc,argv);
    test_opt.retrieveOption(argc,argv);
    fieldname_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    geo_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    buffer_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    // rbox_opt.retrieveOption(argc,argv);
    // cbox_opt.retrieveOption(argc,argv);
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

  if(srcnodata_opt.size()){
    while(srcnodata_opt.size()<bndnodata_opt.size())
      srcnodata_opt.push_back(srcnodata_opt[0]);
    while(bndnodata_opt.size()<srcnodata_opt.size())
      bndnodata_opt.push_back(bndnodata_opt[0]);
  }

  if(verbose_opt[0])
    std::cout << class_opt << std::endl;
  statfactory::StatFactory stat;
  stat.setNoDataValues(srcnodata_opt);
  Vector2d<unsigned int> posdata;
  unsigned long int nsample=0;
  unsigned long int ntotalvalid=0;
  unsigned long int ntotalinvalid=0;
  vector<unsigned long int> nvalid(class_opt.size());
  vector<unsigned long int> ninvalid(class_opt.size());
  for(int it=0;it<nvalid.size();++it){
    nvalid[it]=0;
    ninvalid[it]=0;
  }

  ImgReaderGdal imgReader;
  if(image_opt.empty()){
    std::cerr << "No image dataset provided (use option -i). Use --help for help information";
      exit(0);
  }
  if(output_opt.empty()){
    std::cerr << "No output dataset provided (use option -o). Use --help for help information";
      exit(0);
  }
  try{
    imgReader.open(image_opt[0]);
  }
  catch(std::string errorstring){
    std::cout << errorstring << std::endl;
    exit(0);
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

  if(fieldname_opt.size()<nband){
    std::string bandString=fieldname_opt[0];
    fieldname_opt.clear();
    fieldname_opt.resize(nband);
    for(int iband=0;iband<nband;++iband){
      int theBand=(band_opt.size()) ? band_opt[iband] : iband;
      ostringstream fs;
      fs << bandString << theBand;
      fieldname_opt[iband]=fs.str();
    }
  }

  if(verbose_opt[0])
    std::cout << fieldname_opt << std::endl;
  
  if(verbose_opt[0]>1)
    std::cout << "Number of bands in input image: " << imgReader.nrOfBand() << std::endl;

  OGRFieldType fieldType;
  OGRFieldType labelType;
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
    if( OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType) != NULL
        && EQUAL(OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType),
                 ltype_opt[0].c_str()))
      labelType=(OGRFieldType) iType;
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
    exit(0);
    break;
  }
  switch( labelType ){
  case OFTInteger:
  case OFTReal:
  case OFTRealList:
  case OFTString:
    if(verbose_opt[0]>1)
      std::cout << std::endl << "label type is: " << OGRFieldDefn::GetFieldTypeName(labelType) << std::endl;
    break;
  default:
    cerr << "label type " << OGRFieldDefn::GetFieldTypeName(labelType) << " not supported" << std::endl;
    exit(0);
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
    }
    catch(string errorString){
      cerr << errorString << endl;
    }
    char     **papszOptions=NULL;
    sampleWriterOgr.createLayer("points", imgReader.getProjection(), wkbPoint, papszOptions);
    sampleIsVirtual=true;

    // string fieldName="label";
    // string fieldValue="class";
    // sampleWriterOgr.createField(fieldName,OFTString);
    if(random_opt.size()){
      //create simple random sampling within boundary
      OGRPoint pt;
      double ulx,uly,lrx,lry;
      imgReader.getBoundingBox(ulx,uly,lrx,lry);
      for(unsigned int ipoint=1;ipoint<=random_opt[0];++ipoint){
	OGRFeature *pointFeature;
	pointFeature=sampleWriterOgr.createFeature();
	// pointFeature->SetField(fieldName.c_str(),fieldValue.c_str());
	double theX=ulx+static_cast<double>(rand())/(RAND_MAX)*(lrx-ulx);
	double theY=uly-static_cast<double>(rand())/(RAND_MAX)*(uly-lry);
	pt.setX(theX);
	pt.setY(theY);
	pointFeature->SetGeometry( &pt ); 
	if(sampleWriterOgr.createFeature(pointFeature) != OGRERR_NONE ){
	  cerr << "Failed to create feature in shapefile" << endl;
	  exit( 1 );
	}
	OGRFeature::DestroyFeature(pointFeature);
      }
    }
    else if(grid_opt.size()){
      //create systematic grid of points 
      OGRPoint pt;
      double ulx,uly,lrx,lry;
      imgReader.getBoundingBox(ulx,uly,lrx,lry);
      unsigned int ipoint=0;
      for(double theY=uly-grid_opt[0]/2;theY>lry;theY-=grid_opt[0]){
	for(double theX=ulx+grid_opt[0]/2;theX<lrx;theX+=grid_opt[0]){
	  if(verbose_opt[0]>1)
	    cout << "position: " << theX << " " << theY << endl;
	  OGRFeature *pointFeature;
	  pointFeature=sampleWriterOgr.createFeature();
	  // pointFeature->SetField(fieldName.c_str(),fieldValue.c_str());
	  pt.setX(theX);
	  pt.setY(theY);
	  pointFeature->SetGeometry( &pt ); 
	  if(sampleWriterOgr.createFeature(pointFeature) != OGRERR_NONE ){
	    cerr << "Failed to create feature in shapefile" << endl;
	    exit( 1 );
	  }
	  OGRFeature::DestroyFeature(pointFeature);
	}
      }
    }
    else{
      std::cerr << "No sample dataset provided (use option -s). Use --help for help information";
      exit(0);
    }
    try{
      sampleWriterOgr.close();
      sampleReaderOgr.open("/vsimem/virtual");
    }
    catch(string errorString){
      cerr << errorString << endl;
    }
  }

  if(sampleIsRaster){
    if(class_opt.empty()){
      // std::cout << "Warning: no classes selected, if a classes must be extracted, set to -1 for all classes using option -c -1" << std::endl;
      ImgReaderGdal classReader;
      ImgWriterOgr ogrWriter;
      assert(sample_opt.size());
      classReader.open(sample_opt[0]);
      // vector<int> classBuffer(classReader.nrOfCol());
      vector<double> classBuffer(classReader.nrOfCol());
      Vector2d<double> imgBuffer(nband);//[band][col]
      vector<double> sample(2+nband);//x,y,band values
      Vector2d<double> writeBuffer;
      // vector<int> writeBufferClass;
      vector<double> writeBufferClass;
      vector<int> selectedClass;
      Vector2d<double> selectedBuffer;
      double oldimgrow=-1;
      int irow=0;
      int icol=0;
      if(verbose_opt[0]>1)
        std::cout << "extracting sample from image..." << std::endl;
      progress=0;
      pfnProgress(progress,pszMessage,pProgressArg);
      for(irow=0;irow<classReader.nrOfRow();++irow){
        if(irow%down_opt[0])
          continue;
        // classReader.readData(classBuffer,GDT_Int32,irow);
        classReader.readData(classBuffer,GDT_Float64,irow);
        double x,y;//geo coordinates
        double iimg,jimg;//image coordinates in img image
        for(icol=0;icol<classReader.nrOfCol();++icol){
          if(icol%down_opt[0])
            continue;
          // int theClass=0;
          double theClass=classBuffer[icol];
          // int processClass=-1;
          int processClass=0;
          // if(class_opt[0]<0){//process every class except 0
          //   if(classBuffer[icol]){
          //     processClass=0;
          //     theClass=classBuffer[icol];
          //   }
          // }
          // else{
          //   for(int iclass=0;iclass<class_opt.size();++iclass){
          //     if(classBuffer[icol]==class_opt[iclass]){
          //       processClass=iclass;
          //       theClass=class_opt[iclass];
          //     }
          //   }
          // }
          // if(processClass>=0){
          bool valid=true;
          if(valid){
            if(geo_opt[0]){
              classReader.image2geo(icol,irow,x,y);
              sample[0]=x;
              sample[1]=y;
              if(verbose_opt[0]>1){
                std::cout.precision(12);
                std::cout << theClass << " " << x << " " << y << std::endl;
              }
              //find col in img
              imgReader.geo2image(x,y,iimg,jimg);
              //nearest neighbour
              jimg=static_cast<int>(jimg);
              iimg=static_cast<int>(iimg);
              if(static_cast<int>(iimg)<0||static_cast<int>(iimg)>=imgReader.nrOfCol())
                continue;
            }
            else{
              iimg=icol;
              jimg=irow;
              sample[0]=iimg;
              sample[1]=jimg;
            }
            if(static_cast<int>(jimg)<0||static_cast<int>(jimg)>=imgReader.nrOfRow())
              continue;

            bool valid=true;

            if(static_cast<int>(jimg)!=static_cast<int>(oldimgrow)){
              assert(imgBuffer.size()==nband);
              for(int iband=0;iband<nband;++iband){
		int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                imgReader.readData(imgBuffer[iband],GDT_Float64,static_cast<int>(jimg),theBand);
                assert(imgBuffer[iband].size()==imgReader.nrOfCol());
		if(srcnodata_opt.size()){
		  vector<int>::const_iterator bndit=bndnodata_opt.begin();
		  vector<double>::const_iterator srcit=srcnodata_opt.begin();
		  while(bndit!=bndnodata_opt.end()&&srcit!=srcnodata_opt.end()){
		    if((*bndit==theBand)&&(*srcit==imgBuffer[iband][static_cast<int>(iimg)])){
		      valid=false;
		      break;
		    }
		    else{
		      ++bndit;
		      ++srcit;
		    }
		  }
		}
	      }
              oldimgrow=jimg;
	    }

            if(valid){
              for(int iband=0;iband<imgBuffer.size();++iband){
                if(imgBuffer[iband].size()!=imgReader.nrOfCol()){
                  std::cout << "Error in band " << iband << ": " << imgBuffer[iband].size() << "!=" << imgReader.nrOfCol() << std::endl;
                  assert(imgBuffer[iband].size()==imgReader.nrOfCol());
                }
                sample[iband+2]=imgBuffer[iband][static_cast<int>(iimg)];
              }
              float theThreshold=(threshold_opt.size()>1)?threshold_opt[processClass]:threshold_opt[0];
              if(theThreshold>0){//percentual value
                double p=static_cast<double>(rand())/(RAND_MAX);
                p*=100.0;
                if(p>theThreshold)
		  continue;//do not select for now, go to next column
              }
              else if(nvalid.size()>processClass){//absolute value
                if(nvalid[processClass]>=-theThreshold)
                  continue;//do not select any more pixels for this class, go to next column to search for other classes
              }
	      writeBuffer.push_back(sample);
	      writeBufferClass.push_back(theClass);
	      ++ntotalvalid;
              if(nvalid.size()>processClass)
                ++(nvalid[processClass]);
	    }
            else{
              ++ntotalinvalid;
              if(ninvalid.size()>processClass)
                ++(ninvalid[processClass]);
            }
          }//processClass
        }//icol
        progress=static_cast<float>(irow+1.0)/classReader.nrOfRow();
        pfnProgress(progress,pszMessage,pProgressArg);
      }//irow
      progress=100;
      pfnProgress(progress,pszMessage,pProgressArg);
      if(writeBuffer.size()>0){
        assert(ntotalvalid==writeBuffer.size());
        if(verbose_opt[0]>0)
          std::cout << "creating image sample writer " << output_opt[0] << " with " << writeBuffer.size() << " samples (" << ntotalinvalid << " invalid)" << std::endl;
        ogrWriter.open(output_opt[0],ogrformat_opt[0]);
        char     **papszOptions=NULL;
        ostringstream slayer;
        slayer << "training data";
        std::string layername=slayer.str();
        ogrWriter.createLayer(layername, imgReader.getProjection(), wkbPoint, papszOptions);
        std::string fieldname="fid";//number of the point
        ogrWriter.createField(fieldname,OFTInteger);
        map<std::string,double> pointAttributes;
        ogrWriter.createField(label_opt[0],labelType);
        for(int iband=0;iband<nband;++iband){
	  int theBand=(band_opt.size()) ? band_opt[iband] : iband;
          ogrWriter.createField(fieldname_opt[iband],fieldType);
        }
        std::cout << "writing sample to " << output_opt[0] << "..." << std::endl;
        progress=0;
        pfnProgress(progress,pszMessage,pProgressArg);
        for(int isample=0;isample<writeBuffer.size();++isample){
          if(verbose_opt[0]>1)
            std::cout << "writing sample " << isample << std::endl;
          pointAttributes[label_opt[0]]=writeBufferClass[isample];
          for(int iband=0;iband<writeBuffer[0].size()-2;++iband){
	    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
            // ostringstream fs;
            // if(nband==1)
            //   fs << fieldname_opt[0];
            // else
            //   fs << fieldname_opt[0] << theBand;
            // pointAttributes[fs.str()]=writeBuffer[isample][iband+2];
            pointAttributes[fieldname_opt[iband]]=writeBuffer[isample][iband+2];
          }
          if(verbose_opt[0]>1)
            std::cout << "all bands written" << std::endl;
          ogrWriter.addPoint(writeBuffer[isample][0],writeBuffer[isample][1],pointAttributes,fieldname,isample);
          progress=static_cast<float>(isample+1.0)/writeBuffer.size();
          pfnProgress(progress,pszMessage,pProgressArg);
        }
        ogrWriter.close();
      }
      else{
        std::cout << "No data found for any class " << std::endl;
      }
      classReader.close();
      nsample=writeBuffer.size();
      if(verbose_opt[0])
        std::cout << "total number of samples written: " << nsample << std::endl;
    }
    else{//class_opt.size()!=0
      assert(class_opt[0]);
      //   if(class_opt[0]){
      assert(threshold_opt.size()==1||threshold_opt.size()==class_opt.size());
      ImgReaderGdal classReader;
      ImgWriterOgr ogrWriter;
      if(verbose_opt[0]>1){
        std::cout << "reading position from sample dataset " << std::endl;
        std::cout << "class thresholds: " << std::endl;
        for(int iclass=0;iclass<class_opt.size();++iclass){
          if(threshold_opt.size()>1)
            std::cout << class_opt[iclass] << ": " << threshold_opt[iclass] << std::endl;
          else
            std::cout << class_opt[iclass] << ": " << threshold_opt[0] << std::endl;
        }
      }
      classReader.open(sample_opt[0]);
      vector<int> classBuffer(classReader.nrOfCol());
      // vector<double> classBuffer(classReader.nrOfCol());
      Vector2d<double> imgBuffer(nband);//[band][col]
      vector<double> sample(2+nband);//x,y,band values
      Vector2d<double> writeBuffer;
      vector<int> writeBufferClass;
      // vector<double> writeBufferClass;
      vector<int> selectedClass;
      Vector2d<double> selectedBuffer;
      double oldimgrow=-1;
      int irow=0;
      int icol=0;
      if(verbose_opt[0]>1)
        std::cout << "extracting sample from image..." << std::endl;
      progress=0;
      pfnProgress(progress,pszMessage,pProgressArg);
      for(irow=0;irow<classReader.nrOfRow();++irow){
        if(irow%down_opt[0])
          continue;
        classReader.readData(classBuffer,GDT_Int32,irow);
        // classReader.readData(classBuffer,GDT_Float64,irow);
        double x,y;//geo coordinates
        double iimg,jimg;//image coordinates in img image
        for(icol=0;icol<classReader.nrOfCol();++icol){
          if(icol%down_opt[0])
            continue;
          int theClass=0;
          // double theClass=0;
          int processClass=-1;
          if(class_opt.empty()){//process every class
            if(classBuffer[icol]){
              processClass=0;
              theClass=classBuffer[icol];
            }
          }
          else{
            for(int iclass=0;iclass<class_opt.size();++iclass){
              if(classBuffer[icol]==class_opt[iclass]){
                processClass=iclass;
                theClass=class_opt[iclass];
              }
            }
          }
          if(processClass>=0){
            //         if(classBuffer[icol]==class_opt[0]){
            if(geo_opt[0]){
              classReader.image2geo(icol,irow,x,y);
              sample[0]=x;
              sample[1]=y;
              if(verbose_opt[0]>1){
                std::cout.precision(12);
                std::cout << theClass << " " << x << " " << y << std::endl;
              }
              //find col in img
              imgReader.geo2image(x,y,iimg,jimg);
              //nearest neighbour
              jimg=static_cast<int>(jimg);
              iimg=static_cast<int>(iimg);
              if(static_cast<int>(iimg)<0||static_cast<int>(iimg)>=imgReader.nrOfCol())
                continue;
            }
            else{
              iimg=icol;
              jimg=irow;
              sample[0]=iimg;
              sample[1]=jimg;
            }
            if(static_cast<int>(jimg)<0||static_cast<int>(jimg)>=imgReader.nrOfRow())
              continue;

            bool valid=true;

            if(static_cast<int>(jimg)!=static_cast<int>(oldimgrow)){
              assert(imgBuffer.size()==nband);
              for(int iband=0;iband<nband;++iband){
		int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                imgReader.readData(imgBuffer[iband],GDT_Float64,static_cast<int>(jimg),theBand);
                assert(imgBuffer[iband].size()==imgReader.nrOfCol());

		if(srcnodata_opt.size()){
		  vector<int>::const_iterator bndit=bndnodata_opt.begin();
		  vector<double>::const_iterator srcit=srcnodata_opt.begin();
		  while(bndit!=bndnodata_opt.end()&&srcit!=srcnodata_opt.end()){
		    if((*bndit==theBand)&&(*srcit==imgBuffer[iband][static_cast<int>(iimg)])){
		      valid=false;
		      break;
		    }
		    else{
		      ++bndit;
		      ++srcit;
		    }
		  }
		}
              }
              oldimgrow=jimg;
            }
            if(valid){
              for(int iband=0;iband<imgBuffer.size();++iband){
                if(imgBuffer[iband].size()!=imgReader.nrOfCol()){
                  std::cout << "Error in band " << iband << ": " << imgBuffer[iband].size() << "!=" << imgReader.nrOfCol() << std::endl;
                  assert(imgBuffer[iband].size()==imgReader.nrOfCol());
                }
                sample[iband+2]=imgBuffer[iband][static_cast<int>(iimg)];
              }
              float theThreshold=(threshold_opt.size()>1)?threshold_opt[processClass]:threshold_opt[0];
              if(theThreshold>0){//percentual value
                double p=static_cast<double>(rand())/(RAND_MAX);
                p*=100.0;
                if(p>theThreshold)
                  continue;//do not select for now, go to next column
              }
              else if(nvalid.size()>processClass){//absolute value
                if(nvalid[processClass]>=-theThreshold)
                  continue;//do not select any more pixels for this class, go to next column to search for other classes
              }
              writeBuffer.push_back(sample);
              //             writeBufferClass.push_back(class_opt[processClass]);
              writeBufferClass.push_back(theClass);
              ++ntotalvalid;
              if(nvalid.size()>processClass)
                ++(nvalid[processClass]);
            }
            else{
              ++ntotalinvalid;
              if(ninvalid.size()>processClass)
                ++(ninvalid[processClass]);
            }
          }//processClass
        }//icol
        progress=static_cast<float>(irow+1.0)/classReader.nrOfRow();
        pfnProgress(progress,pszMessage,pProgressArg);
      }//irow
      if(writeBuffer.size()>0){
        assert(ntotalvalid==writeBuffer.size());
        if(verbose_opt[0]>0)
          std::cout << "creating image sample writer " << output_opt[0] << " with " << writeBuffer.size() << " samples (" << ntotalinvalid << " invalid)" << std::endl;
        ogrWriter.open(output_opt[0],ogrformat_opt[0]);
        char     **papszOptions=NULL;
        ostringstream slayer;
        slayer << "training data";
        std::string layername=slayer.str();
        ogrWriter.createLayer(layername, imgReader.getProjection(), wkbPoint, papszOptions);
        std::string fieldname="fid";//number of the point
        ogrWriter.createField(fieldname,OFTInteger);
        map<std::string,double> pointAttributes;
        //         ogrWriter.createField(label_opt[0],OFTInteger);
        ogrWriter.createField(label_opt[0],labelType);
        for(int iband=0;iband<nband;++iband){
	  int theBand=(band_opt.size()) ? band_opt[iband] : iband;
          // ostringstream fs;
          // if(nband==1)
          //   fs << fieldname_opt[0];
          // else
          //   fs << fieldname_opt[0] << theBand;
          // ogrWriter.createField(fs.str(),fieldType);
          ogrWriter.createField(fieldname_opt[iband],fieldType);
        }
        pfnProgress(progress,pszMessage,pProgressArg);
        std::cout << "writing sample to " << output_opt[0] << "..." << std::endl;
        progress=0;
        pfnProgress(progress,pszMessage,pProgressArg);
        for(int isample=0;isample<writeBuffer.size();++isample){
          pointAttributes[label_opt[0]]=writeBufferClass[isample];
          for(int iband=0;iband<writeBuffer[0].size()-2;++iband){
	    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
            // ostringstream fs;
            // if(nband==1)
            //   fs << fieldname_opt[0];
            // else
            //   fs << fieldname_opt[0] << theBand;
            // pointAttributes[fs.str()]=writeBuffer[isample][iband+2];
            pointAttributes[fieldname_opt[iband]]=writeBuffer[isample][iband+2];
          }
          ogrWriter.addPoint(writeBuffer[isample][0],writeBuffer[isample][1],pointAttributes,fieldname,isample);
          progress=static_cast<float>(isample+1.0)/writeBuffer.size();
          pfnProgress(progress,pszMessage,pProgressArg);
        }
        ogrWriter.close();
      }
      else{
        std::cout << "No data found for any class " << std::endl;
      }
      classReader.close();
      nsample=writeBuffer.size();
      if(verbose_opt[0]){
        std::cout << "total number of samples written: " << nsample << std::endl;
        if(nvalid.size()==class_opt.size()){
          for(int iclass=0;iclass<class_opt.size();++iclass)
            std::cout << "class " << class_opt[iclass] << " has " << nvalid[iclass] << " samples" << std::endl;
        }
      }
    }
  }
  else{//vector dataset
    if(verbose_opt[0]>1)
      std::cout << "creating image sample writer " << output_opt[0] << std::endl;
    ImgWriterOgr ogrWriter;
    ImgWriterOgr ogrTestWriter;
    double e_ulx;
    double e_uly;
    double e_lrx;
    double e_lry;
    try{
      sampleReaderOgr.getExtent(e_ulx,e_uly,e_lrx,e_lry);
      bool hasCoverage=((e_ulx < imgReader.getLrx())&&(e_lrx > imgReader.getUlx())&&(e_lry < imgReader.getUly())&&(e_uly > imgReader.getLry()));
      if(!hasCoverage){
	ostringstream ess;
	ess << "No coverage in " << image_opt[0] << " for any layer in " << sample_opt[0] << endl;
	throw(ess.str());
      }
      ogrWriter.open(output_opt[0],ogrformat_opt[0]);
      if(test_opt.size()){
	if(verbose_opt[0]>1)
	  std::cout << "creating image test writer " << test_opt[0] << std::endl;
	ogrTestWriter.open(test_opt[0],ogrformat_opt[0]);
      }

      //if class_opt not set, get number of classes from input image for these rules
      switch(ruleMap[rule_opt[0]]){
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
      sampleReaderOgr.getExtent(e_ulx,e_uly,e_lrx,e_lry,ilayer);
      bool hasCoverage=((e_ulx < imgReader.getLrx())&&(e_lrx > imgReader.getUlx())&&(e_lry < imgReader.getUly())&&(e_uly > imgReader.getLry()));
      if(!hasCoverage)
	continue;
      float theThreshold=(threshold_opt.size()==layer_opt.size())? threshold_opt[layerIndex]: threshold_opt[0];
      cout << "processing layer " << currentLayername << endl;
      
      readLayer->ResetReading();
      OGRLayer *writeLayer;
      OGRLayer *writeTestLayer;

      if(polygon_opt[0]){
	if(verbose_opt[0])
	  std::cout << "create polygons" << std::endl;
	char **papszOptions=NULL;
	writeLayer=ogrWriter.createLayer(readLayer->GetName(), imgReader.getProjection(), wkbPolygon, papszOptions);
	if(test_opt.size())
	  writeTestLayer=ogrTestWriter.createLayer(readLayer->GetName(), imgReader.getProjection(), wkbPolygon, papszOptions);
      }
      else{
	if(verbose_opt[0])
	  std::cout << "create points in layer " << readLayer->GetName() << std::endl;
	char **papszOptions=NULL;

	writeLayer=ogrWriter.createLayer(readLayer->GetName(), imgReader.getProjection(), wkbPoint, papszOptions);
	if(test_opt.size()){
	  char **papszOptions=NULL;
	  writeTestLayer=ogrTestWriter.createLayer(readLayer->GetName(), imgReader.getProjection(), wkbPoint, papszOptions);
	}
      }
      if(verbose_opt[0])
	std::cout << "copy fields from layer " << ilayer << std::flush << std::endl;
      ogrWriter.copyFields(sampleReaderOgr,ilayer,ilayerWrite);

      if(test_opt.size()){
	if(verbose_opt[0])
	  std::cout << "copy fields test writer" << std::flush << std::endl;
	ogrTestWriter.copyFields(sampleReaderOgr,ilayer,ilayerWrite);
      }
      // vector<std::string> fieldnames;
      // if(verbose_opt[0])
      // 	std::cout << "get fields" << std::flush << std::endl;
      // sampleReaderOgr.getFields(fieldnames);
      // assert(fieldnames.size()==ogrWriter.getFieldCount(ilayerWrite));
      // map<std::string,double> pointAttributes;

      if(class_opt.size()){
	switch(ruleMap[rule_opt[0]]){
	case(rule::proportion)://proportion for each class
	case(rule::count):{//count for each class
	  for(int iclass=0;iclass<class_opt.size();++iclass){
	    ostringstream cs;
	    cs << class_opt[iclass];
	    ogrWriter.createField(cs.str(),fieldType,ilayerWrite);
	  }
	  break;
	}
	case(rule::custom):
	case(rule::mode):
	  ogrWriter.createField(label_opt[0],fieldType,ilayerWrite);
	if(test_opt.size())
	  ogrTestWriter.createField(label_opt[0],fieldType,ilayerWrite);
	break;
	}
      }
      else{
	for(int iband=0;iband<nband;++iband){
	  int theBand=(band_opt.size()) ? band_opt[iband] : iband;
	  ostringstream fs;
	  fs << fieldname_opt[iband];
	  if(verbose_opt[0]>1)
	    std::cout << "creating field " << fs.str() << std::endl;

	  ogrWriter.createField(fs.str(),fieldType,ilayerWrite);
	  if(test_opt.size())
	    ogrTestWriter.createField(fs.str(),fieldType,ilayerWrite);
	}
      }
      OGRFeature *readFeature;
      unsigned long int ifeature=0;
      unsigned long int nfeatureLayer=sampleReaderOgr.getFeatureCount(ilayer);
      unsigned long int ntotalvalidLayer=0;
      progress=0;
      pfnProgress(progress,pszMessage,pProgressArg);
      while( (readFeature = readLayer->GetNextFeature()) != NULL ){
	bool validFeature=false;
	bool writeTest=false;//write this feature to test_opt[0] instead of output_opt
	if(verbose_opt[0]>0)
	  std::cout << "reading feature " << readFeature->GetFID() << std::endl;
	if(theThreshold>0){//percentual value
	  // if(!test_opt.size()&&ntotalvalid>threshold_opt[0]/100.0*nfeature)
	  //   break;
	  double p=static_cast<double>(rand())/(RAND_MAX);
	  p*=100.0;
	  if(p>theThreshold){
	    if(test_opt.size())
	      writeTest=true;
	    else
	      continue;//do not select for now, go to next feature
	  }
	}
	else{//absolute value
	  if(threshold_opt.size()==layer_opt.size()){
	    if(ntotalvalidLayer>=-theThreshold){
	      if(test_opt.size())
		writeTest=true;
	      else
		continue;//do not select any more pixels, go to next column feature
	    }
	  }
	  else{
	    if(ntotalvalid>=-theThreshold){
	      if(test_opt.size())
		writeTest=true;
	      else
		continue;//do not select any more pixels, go to next column feature
	    }
	  }
	}
	if(verbose_opt[0]>0)
	  std::cout << "processing feature " << readFeature->GetFID() << std::endl;
	//get x and y from readFeature
	double x,y;
	OGRGeometry *poGeometry;
	poGeometry = readFeature->GetGeometryRef();
	assert(poGeometry!=NULL);
	try{
	  if(wkbFlatten(poGeometry->getGeometryType()) == wkbPoint ){

	    if(!buffer_opt.size()){
	      switch(ruleMap[rule_opt[0]]){
	      case(rule::point):
	      case(rule::centroid):
		buffer_opt.push_back(1);//default
	      break;
	      default:
		buffer_opt.push_back(3);//default
	      }
	    }

	    if(verbose_opt[0]>1)
	      std::cout << "boundary: " << buffer_opt[0] << std::endl;

	    OGRPolygon writePolygon;
	    OGRLinearRing writeRing;
	    OGRPoint writeCentroidPoint;
	    OGRFeature *writePolygonFeature;
	    OGRFeature *writeCentroidFeature;

	    OGRPoint *poPoint = (OGRPoint *) poGeometry;
	    writeCentroidPoint=*poPoint;

	    x=poPoint->getX();
	    y=poPoint->getY();

	    double i_centre,j_centre;
	    if(geo_opt[0])
	      imgReader.geo2image(x,y,i_centre,j_centre);
	    else{
	      i_centre=x;
	      j_centre=y;
	    }
	    //nearest neighbour
	    j_centre=static_cast<int>(j_centre);
	    i_centre=static_cast<int>(i_centre);

	    double uli=i_centre-buffer_opt[0]/2;
	    double ulj=j_centre-buffer_opt[0]/2;
	    double lri=i_centre+buffer_opt[0]/2;
	    double lrj=j_centre+buffer_opt[0]/2;

	    //nearest neighbour
	    ulj=static_cast<int>(ulj);
	    uli=static_cast<int>(uli);
	    lrj=static_cast<int>(lrj);
	    lri=static_cast<int>(lri);

	    // if((polygon_opt[0]&&ruleMap[rule_opt[0]]==rule::point)||(ruleMap[rule_opt[0]]==rule::centroid)){
	    //   uli=i_centre;
	    //   ulj=j_centre;
	    //   lri=i_centre;
	    //   lrj=j_centre;
	    // }

	    //check if j is out of bounds
	    if(static_cast<int>(ulj)<0||static_cast<int>(ulj)>=imgReader.nrOfRow())
	      continue;
	    //check if j is out of bounds
	    if(static_cast<int>(uli)<0||static_cast<int>(lri)>=imgReader.nrOfCol())
	      continue;

	    OGRPoint ulPoint,urPoint,llPoint,lrPoint;
	    double ulx,uly;
	    double urx,ury;

	    if(polygon_opt[0]){
	      if(disc_opt[0]){
		double gt[6];// { 444720, 30, 0, 3751320, 0, -30 };
		double radius=buffer_opt[0]/2.0*sqrt(imgReader.getDeltaX()*imgReader.getDeltaY());
		unsigned short nstep = 25;
		for(int i=0;i<nstep;++i){
		  OGRPoint aPoint;
		  aPoint.setX(x+imgReader.getDeltaX()/2.0+radius*cos(2*PI*i/nstep));
		  aPoint.setY(y-imgReader.getDeltaY()/2.0+radius*sin(2*PI*i/nstep));
		  writeRing.addPoint(&aPoint);
		}
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
	      }
	      else{
		double llx,lly;
		double lrx,lry;
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
	    }

	    if((polygon_opt[0]&&ruleMap[rule_opt[0]]==rule::point)||(ruleMap[rule_opt[0]]==rule::centroid)){
	      uli=i_centre;
	      ulj=j_centre;
	      lri=i_centre;
	      lrj=j_centre;
	    }

	    int nPointWindow=0;//similar to nPointPolygon for polygons
	    if(polygon_opt[0]){
	      if(writeTest)
		writePolygonFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    else if(ruleMap[rule_opt[0]]!=rule::point){
	      if(writeTest)
		writeCentroidFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writeCentroidFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    Vector2d<double> windowValues;
	    vector<double> windowClassValues;

	    if(class_opt.size()){
	      windowClassValues.resize(class_opt.size());
	      //initialize
	      for(int iclass=0;iclass<class_opt.size();++iclass)
		windowClassValues[iclass]=0;
	    }
	    else
	      windowValues.resize(nband);
	    vector< Vector2d<double> > readValues(nband);
	    for(int iband=0;iband<nband;++iband){
	      int theBand=(band_opt.size()) ? band_opt[iband] : iband;
	      imgReader.readDataBlock(readValues[iband],GDT_Float64,uli,lri,ulj,lrj,theBand);
	    }

	    OGRPoint thePoint;
	    for(int j=ulj;j<=lrj;++j){
	      for(int i=uli;i<=lri;++i){
		//check if within raster image
		if(i<0||i>=imgReader.nrOfCol())
		  continue;
		if(j<0||j>=imgReader.nrOfRow())
		  continue;
		//no need to check if point is on surface
		double theX=0;
		double theY=0;
		imgReader.image2geo(i,j,theX,theY);
		thePoint.setX(theX);
		thePoint.setY(theY);

		if(disc_opt[0]&&buffer_opt[0]>1){
		  double radius=buffer_opt[0]/2.0*sqrt(imgReader.getDeltaX()*imgReader.getDeltaY());
		  if((theX-x)*(theX-x)+(theY-y)*(theY-y)>radius*radius)
		    continue;
		}
		bool valid=true;

		if(srcnodata_opt.size()){
		  for(int vband=0;vband<bndnodata_opt.size();++vband){
		    double value=((readValues[bndnodata_opt[vband]])[j-ulj])[i-uli];
		    if(value==srcnodata_opt[vband]){
		      valid=false;
		      break;
		    }
		  }
		}

		if(!valid)
		  continue;
		else
		  validFeature=true;

		// writeRing.addPoint(&thePoint);//already done

		++nPointWindow;
		OGRFeature *writePointFeature;
		if(!polygon_opt[0]){
		  //create feature
		  if(ruleMap[rule_opt[0]]==rule::point){//do not create in case of mean, stdev, median, sum or centroid (only create point at centroid)
		    if(writeTest)
		      writePointFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
		    else
		      writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
		    if(verbose_opt[0]>1)
		      std::cout << "copying fields from polygons " << std::endl;
		    //Geometry of readFeature and writePointFeature are both wkbPoint
		    //attributes AND geometry are copied with SetFrom
		    //test
		    // writePointFeature->SetGeometry(&thePoint);
		    if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
		      cerr << "writing feature failed" << std::endl;

		    assert(wkbFlatten(writePointFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
		    // OGRGeometry *updateGeometry;
		    // updateGeometry = writePointFeature->GetGeometryRef();
		    // OGRPoint *poPoint = (OGRPoint *) updateGeometry;
		    if(verbose_opt[0]>1)
		      std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields" << std::endl;
		  }
		}
		if(class_opt.size()){
		  short value=((readValues[0])[j-ulj])[i-uli];
		  for(int iclass=0;iclass<class_opt.size();++iclass){
		    if(value==class_opt[iclass])
		      windowClassValues[iclass]+=1;
		  }
		}
		else{
		  for(int iband=0;iband<nband;++iband){
		    int theBand=(band_opt.size()) ? band_opt[iband] : iband;
		    assert(j-ulj>=0);
		    assert(j-ulj<readValues[iband].size());
		    assert(i-uli>=0);
		    assert(i-uli<((readValues[iband])[j-ulj]).size());
		    double value=((readValues[iband])[j-ulj])[i-uli];
		    // imgReader.readData(value,GDT_Float64,i,j,theBand);
		    if(verbose_opt[0]>1)
		      std::cout << ": " << value << std::endl;
		    if(polygon_opt[0]||ruleMap[rule_opt[0]]!=rule::point){
		      windowValues[iband].push_back(value);
		    }
		    else{
		      try{
			if(verbose_opt[0]>1)
			  std::cout << "set field " << fieldname_opt[iband] << " to " << value << std::endl;
			switch( fieldType ){
			case OFTInteger:
			case OFTReal:
			  writePointFeature->SetField(fieldname_opt[iband].c_str(),value);
			  break;
			case OFTString:
			  writePointFeature->SetField(fieldname_opt[iband].c_str(),type2string<double>(value).c_str());
			  break;
			default://not supported
			  std::string errorString="field type not supported";
			  throw(errorString);
			  break;
			}
		      }
		      catch(std::string e){
			std::cout << e << std::endl;
			exit(1);
		      }
		    }//else
		  }//iband
		}//else (class_opt.size())
		if(!polygon_opt[0]){
		  if(ruleMap[rule_opt[0]]==rule::point){//do not create in case of mean or median value (only at centroid)
		    //write feature
		    if(verbose_opt[0]>1)
		      std::cout << "creating point feature" << std::endl;
		    if(writeTest){
		      if(writeTestLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			std::string errorString="Failed to create feature in test ogr vector dataset";
			throw(errorString);
		      }
		    }
		    else{
		      if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			std::string errorString="Failed to create feature in ogr vector dataset";
			throw(errorString);
		      }
		    }
 		    //destroy feature
		    OGRFeature::DestroyFeature( writePointFeature );
		    ++ntotalvalid;
		    ++ntotalvalidLayer;
		  }
		}
	      }
	    }
	    if(polygon_opt[0]||ruleMap[rule_opt[0]]!=rule::point){
	      //do not create if no points found within polygon
	      if(!nPointWindow){
		if(verbose_opt[0])
		  cout << "no points found in window, continuing" << endl;
		continue;
	      }
	      //add ring to polygon
	      if(polygon_opt[0]){
		// writePolygon.addRing(&writeRing);//already done
		// writePolygon.closeRings();//already done
		//write geometry of writePolygon
		//test
		// writePolygonFeature->SetGeometry(&writePolygon);
		if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		//test
		writePolygonFeature->SetGeometry(&writePolygon);
		assert(wkbFlatten(writePolygonFeature->GetGeometryRef()->getGeometryType()) == wkbPolygon);

		if(verbose_opt[0]>1)
		  std::cout << "copying new fields write polygon " << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
		//write polygon feature
	      }
	      else{//write value of polygon to centroid point
		//create feature
		if(verbose_opt[0]>1)
		  std::cout << "copying fields from polygons " << std::endl;
		//test
		//Geometry of readFeature and writeCentroidFeature are both wkbPoint
		//attributes AND geometry are copied with SetFrom
		// writeCentroidFeature->SetGeometry(&writeCentroidPoint);
		if(writeCentroidFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		assert(wkbFlatten(writeCentroidFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
		//test
		// OGRGeometry *updateGeometry;
		// updateGeometry = writeCentroidFeature->GetGeometryRef();
		// assert(wkbFlatten(updateGeometry->getGeometryType()) == wkbPoint );
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writeCentroidFeature->GetFieldCount() << " fields" << std::endl;
	      }
	      if(class_opt.empty()){
		if(ruleMap[rule_opt[0]]==rule::point){//value at centroid of polygon
		  if(verbose_opt[0])
		    std::cout << "number of points in window: " << nPointWindow << std::endl;
		  for(int index=0;index<windowValues.size();++index){
		    //test
		    if(windowValues[index].size()!=1){
		      cerr << "Error: windowValues[index].size()=" << windowValues[index].size() << endl;
		      assert(windowValues[index].size()==1);
		    }
		    double theValue=windowValues[index].back();

		    if(verbose_opt[0])
		      std::cout << "number of points in window: " << nPointWindow << std::endl;
		    int theBand=(band_opt.size()) ? band_opt[index] : index;

		    try{
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
			break;
		      case OFTString:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			break;
		      default://not supported
			std::string errorString="field type not supported";
			throw(errorString);
			break;
		      }
		    }
		    catch(std::string e){
		      std::cout << e << std::endl;
		      exit(1);
		    }
		  }
		}
		else{//ruleMap[rule_opt[0]] is not rule::point
		  double theValue=0;
		  for(int index=0;index<windowValues.size();++index){
		    try{
		      if(ruleMap[rule_opt[0]]==rule::mean)
			theValue=stat.mean(windowValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::stdev)
			theValue=sqrt(stat.var(windowValues[index]));
		      else if(ruleMap[rule_opt[0]]==rule::median)
			theValue=stat.median(windowValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::percentile)
			theValue=stat.percentile(windowValues[index],windowValues[index].begin(),windowValues[index].end(),percentile_opt[0]);
		      else if(ruleMap[rule_opt[0]]==rule::sum)
			theValue=stat.sum(windowValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::max)
			theValue=stat.mymax(windowValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::min)
			theValue=stat.mymin(windowValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::centroid){
			if(verbose_opt[0])
			  std::cout << "number of points in polygon: " << nPointWindow << std::endl;
			assert(nPointWindow<=1);
			assert(nPointWindow==windowValues[index].size());
			theValue=windowValues[index].back();
		      }
		      else{
			std::string errorString="rule not supported";
			throw(errorString);
		      }
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
			break;
		      case OFTString:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			break;
		      default://not supported
			std::string errorString="field type not supported";
			throw(errorString);
			break;
		      }
		    }
		    catch(std::string e){
		      std::cout << e << std::endl;
		      exit(1);
		    }
		  }
		}
	      }
	      else{//class_opt is set
		if(ruleMap[rule_opt[0]]==rule::proportion||ruleMap[rule_opt[0]]==rule::count){
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointWindow << std::endl;
		  if(ruleMap[rule_opt[0]]==rule::proportion)
		    stat.normalize_pct(windowClassValues);
		  for(int index=0;index<windowClassValues.size();++index){
		    double theValue=windowClassValues[index];
		    ostringstream fs;
		    fs << class_opt[index];
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		  }
		}
		else if(ruleMap[rule_opt[0]]==rule::custom){
		  assert(polygon_opt[0]);//not implemented for points
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointWindow << std::endl;
		  stat.normalize_pct(windowClassValues);
		  assert(windowClassValues.size()==2);//11:broadleaved, 12:coniferous
		  if(windowClassValues[0]>=75)//broadleaved
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(11));
		  else if(windowClassValues[1]>=75)//coniferous
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(12));
		  else if(windowClassValues[0]>25&&windowClassValues[1]>25)//mixed
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(13));
		  else{
		    if(verbose_opt[0]){
		      std::cout << "No valid value in windowClassValues..." << std::endl;
		      for(int index=0;index<windowClassValues.size();++index){
			double theValue=windowClassValues[index];
			std::cout << theValue << " ";
		      }
		      std::cout << std::endl;
		    }
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(20));
		  }
		}
		else if(ruleMap[rule_opt[0]]==rule::mode){
		  //maximum votes in polygon
		  if(verbose_opt[0])
		    std::cout << "number of points in window: " << nPointWindow << std::endl;
		  //search for class with maximum votes
		  int maxClass=stat.mymin(class_opt);
		  vector<double>::iterator maxit;
		  maxit=stat.mymax(windowClassValues,windowClassValues.begin(),windowClassValues.end());
		  int maxIndex=distance(windowClassValues.begin(),maxit);
		  maxClass=class_opt[maxIndex];
		  if(verbose_opt[0]>0)
		    std::cout << "maxClass: " << maxClass << std::endl;
		  if(polygon_opt[0])
		    writePolygonFeature->SetField(label_opt[0].c_str(),maxClass);
		  else
		    writeCentroidFeature->SetField(label_opt[0].c_str(),maxClass);
		}
	      }
	      if(polygon_opt[0]){
		if(verbose_opt[0]>1)
		  std::cout << "creating polygon feature" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		else{
		  if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writePolygonFeature );
		++ntotalvalid;
		++ntotalvalidLayer;
	      }
	      else{
		if(verbose_opt[0]>1)
		  std::cout << "creating point feature in centroid" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		else{
		  //test
		  assert(validFeature);
		  if(writeLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writeCentroidFeature );
		++ntotalvalid;
		++ntotalvalidLayer;
	      }
	    }
	  }//if wkbPoint
	  else if(wkbFlatten(poGeometry->getGeometryType()) == wkbPolygon){
            
	    OGRPolygon readPolygon = *((OGRPolygon *) poGeometry);
	    OGRPolygon writePolygon;
	    OGRLinearRing writeRing;
	    OGRPoint writeCentroidPoint;
	    OGRFeature *writePolygonFeature;
	    OGRFeature *writeCentroidFeature;

	    readPolygon.closeRings();

	    if(verbose_opt[0]>1)
	      std::cout << "get point on polygon" << std::endl;
	    if(ruleMap[rule_opt[0]]==rule::centroid)
	      readPolygon.Centroid(&writeCentroidPoint);
	    else if(readPolygon.PointOnSurface(&writeCentroidPoint)!=OGRERR_NONE){
	      // cerr << "function PointOnSurface failed, trying centroid instead" << endl;
	      readPolygon.Centroid(&writeCentroidPoint);
	    }
	    double ulx,uly,lrx,lry;
	    double uli,ulj,lri,lrj;
	    if((polygon_opt[0]&&ruleMap[rule_opt[0]]==rule::point)||(ruleMap[rule_opt[0]]==rule::centroid)){
	      ulx=writeCentroidPoint.getX();
	      uly=writeCentroidPoint.getY();
	      lrx=ulx;
	      lry=uly;
	    }
	    else{
	      //get envelope
	      if(verbose_opt[0])
		std::cout << "reading envelope for polygon " << ifeature << std::endl;
	      OGREnvelope* psEnvelope=new OGREnvelope();
	      readPolygon.getEnvelope(psEnvelope);
	      ulx=psEnvelope->MinX;
	      uly=psEnvelope->MaxY;
	      lrx=psEnvelope->MaxX;
	      lry=psEnvelope->MinY;
	      delete psEnvelope;
	    }
	    if(geo_opt[0]){
	      imgReader.geo2image(ulx,uly,uli,ulj);
	      imgReader.geo2image(lrx,lry,lri,lrj);
	    }
	    else{
	      uli=ulx;
	      ulj=uly;
	      lri=lrx;
	      lrj=lry;
	    }
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
	    // if(uli<0||lri>=imgReader.nrOfCol()||ulj<0||lrj>=imgReader.nrOfRow())
	    //   continue;

	    int nPointPolygon=0;

	    if(polygon_opt[0]){
	      if(writeTest)
		writePolygonFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    else if(ruleMap[rule_opt[0]]!=rule::point){
	      if(writeTest)
		writeCentroidFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writeCentroidFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    // vector<double> polyValues;
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
	    vector< Vector2d<double> > readValues(nband);
	    for(int iband=0;iband<nband;++iband){
	      int theBand=(band_opt.size()) ? band_opt[iband] : iband;
	      imgReader.readDataBlock(readValues[iband],GDT_Float64,uli,lri,ulj,lrj,theBand);
	    }

	    OGRPoint thePoint;
	    for(int j=ulj;j<=lrj;++j){
	      for(int i=uli;i<=lri;++i){
		//check if within raster image
		if(i<0||i>=imgReader.nrOfCol())
		  continue;
		if(j<0||j>=imgReader.nrOfRow())
		  continue;
		//check if point is on surface
		double theX=0;
		double theY=0;
		imgReader.image2geo(i,j,theX,theY);
		thePoint.setX(theX);
		thePoint.setY(theY);
		if(ruleMap[rule_opt[0]]!=rule::centroid&&!readPolygon.Contains(&thePoint))
		  continue;

		bool valid=true;

		if(srcnodata_opt.size()){
		  for(int vband=0;vband<bndnodata_opt.size();++vband){
		    double value=((readValues[bndnodata_opt[vband]])[j-ulj])[i-uli];
		    if(value==srcnodata_opt[vband]){
		      valid=false;
		      break;
		    }
		  }
		}

		if(!valid)
		  continue;
		else
		  validFeature=true;

		writeRing.addPoint(&thePoint);//todo: check if I need to add all interior points to ring or do I need to check if point is on ring first?
		// if(writeRing.isPointOnRingBoundary(&thePoint))
		//    writeRing.addPoint(&thePoint);
		if(verbose_opt[0]>1)
		  std::cout << "point is on surface:" << thePoint.getX() << "," << thePoint.getY() << std::endl;
		++nPointPolygon;

		if(polythreshold_opt.size())
		  if(nPointPolygon>polythreshold_opt[0])
		    continue;
		// throw(nPointPolygon);
		OGRFeature *writePointFeature;
		if(!polygon_opt[0]){
		  //create feature
		  if(ruleMap[rule_opt[0]]==rule::point){//do not create in case of mean, stdev, median, sum or centroid (only create point at centroid)
		    if(writeTest)
		      writePointFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
		    else
		      writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
		    if(verbose_opt[0]>1)
		      std::cout << "copying fields from polygons " << std::endl;
		    if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
		      cerr << "writing feature failed" << std::endl;
		    if(verbose_opt[0]>1)
		      std::cout << "set geometry as point " << std::endl;
		    //test
		    writePointFeature->SetGeometry(&thePoint);
		    assert(wkbFlatten(writePointFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
		    // OGRGeometry *updateGeometry;
		    // updateGeometry = writePointFeature->GetGeometryRef();
		    // OGRPoint *poPoint = (OGRPoint *) updateGeometry;
		    if(verbose_opt[0]>1)
		      std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields" << std::endl;
		  }
		}
		if(class_opt.size()){
		  short value=((readValues[0])[j-ulj])[i-uli];
		  for(int iclass=0;iclass<class_opt.size();++iclass){
		    if(value==class_opt[iclass])
		      polyClassValues[iclass]+=1;
		  }
		}
		else{
		  for(int iband=0;iband<nband;++iband){
		    double value=((readValues[iband])[j-ulj])[i-uli];
		    if(verbose_opt[0]>1)
		      std::cout << ": " << value << std::endl;
		    if(polygon_opt[0]||ruleMap[rule_opt[0]]!=rule::point)
		      polyValues[iband].push_back(value);
		    else{
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[iband] << " to " << value << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			writePointFeature->SetField(fieldname_opt[iband].c_str(),value);
			break;
		      case OFTString:
			writePointFeature->SetField(fieldname_opt[iband].c_str(),type2string<double>(value).c_str());
			break;
		      default://not supported
			assert(0);
			break;
		      }
		    }//else
		  }//iband
		}//else (class_opt.size())
		if(!polygon_opt[0]){
		  if(ruleMap[rule_opt[0]]==rule::point){//do not create in case of mean or median value (only at centroid)
		    //write feature
		    if(verbose_opt[0]>1)
		      std::cout << "creating point feature" << std::endl;
		    if(writeTest){
		      if(writeTestLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			std::string errorString="Failed to create feature in test ogr vector dataset";
			throw(errorString);
		      }
		    }
		    else{
		      if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			std::string errorString="Failed to create feature in ogr vector dataset";
			throw(errorString);
		      }
		    }
		    //destroy feature
		    OGRFeature::DestroyFeature( writePointFeature );
		    ++ntotalvalid;
		    ++ntotalvalidLayer;
		  }
		}
	      }
	    }
	    if(polygon_opt[0]||ruleMap[rule_opt[0]]!=rule::point){
	      //do not create if no points found within polygon
	      if(!nPointPolygon){
		if(verbose_opt[0])
		  cout << "no points found in polygon, continuing" << endl;
		continue;
	      }
	      //add ring to polygon
	      if(polygon_opt[0]){
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
		//write geometry of writePolygon
		//test
		//writePolygonFeature and readFeature are both of type wkbPolygon
		// writePolygonFeature->SetGeometry(&writePolygon);
		if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "copying new fields write polygon " << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
		//write polygon feature
	      }
	      else{//write value of polygon to centroid point
		//create feature
		if(verbose_opt[0]>1)
		  std::cout << "copying fields from polygons " << std::endl;
		if(writeCentroidFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		writeCentroidFeature->SetGeometry(&writeCentroidPoint);
		assert(wkbFlatten(writeCentroidFeature->GetGeometryRef()->getGeometryType()) == wkbPoint );
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writeCentroidFeature->GetFieldCount() << " fields" << std::endl;
	      }
	      if(class_opt.empty()){
		if(ruleMap[rule_opt[0]]==rule::point){//value at centroid of polygon
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  for(int index=0;index<polyValues.size();++index){
		    assert(polyValues[index].size()==1);
		    double theValue=polyValues[index].back();

		    if(verbose_opt[0])
		      std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		    int theBand=(band_opt.size()) ? band_opt[index] : index;
		    try{
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
			break;
		      case OFTString:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			break;
		      default://not supported
			std::string errorString="field type not supported";
			throw(errorString);
			break;
		      }
		    }
		    catch(std::string e){
		      std::cout << e << std::endl;
		      exit(1);
		    }
		  }
		}
		else{//ruleMap[rule_opt[0]] is not rule::point
		  double theValue=0;
		  for(int index=0;index<polyValues.size();++index){
		    try{
		      if(ruleMap[rule_opt[0]]==rule::mean)
			theValue=stat.mean(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::stdev)
			theValue=sqrt(stat.var(polyValues[index]));
		      else if(ruleMap[rule_opt[0]]==rule::median)
			theValue=stat.median(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::percentile)
			theValue=stat.percentile(polyValues[index],polyValues[index].begin(),polyValues[index].end(),percentile_opt[0]);
		      else if(ruleMap[rule_opt[0]]==rule::sum)
			theValue=stat.sum(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::max)
			theValue=stat.mymax(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::min)
			theValue=stat.mymin(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::centroid){
			if(verbose_opt[0])
			  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
			assert(nPointPolygon<=1);
			assert(nPointPolygon==polyValues[index].size());
			theValue=polyValues[index].back();
		      }
		      else{
			std::string errorString="rule not supported";
			throw(errorString);
		      }
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
			break;
		      case OFTString:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			break;
		      default:
			std::string errorString="field type not supported";
			throw(errorString);
			break;
		      }
		    }
		    catch(std::string e){
		      std::cout << e << std::endl;
		      exit(1);
		    }
		  }
		}
	      }
	      else{//class_opt is set
		if(ruleMap[rule_opt[0]]==rule::proportion||ruleMap[rule_opt[0]]==rule::count){
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  if(ruleMap[rule_opt[0]]==rule::proportion)
		    stat.normalize_pct(polyClassValues);
		  for(int index=0;index<polyClassValues.size();++index){
		    double theValue=polyClassValues[index];
		    ostringstream fs;
		    fs << class_opt[index];
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		  }
		}
		else if(ruleMap[rule_opt[0]]==rule::custom){
		  assert(polygon_opt[0]);//not implemented for points
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  stat.normalize_pct(polyClassValues);
		  assert(polyClassValues.size()==2);//11:broadleaved, 12:coniferous
		  if(polyClassValues[0]>=75)//broadleaved
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(11));
		  else if(polyClassValues[1]>=75)//coniferous
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(12));
		  else if(polyClassValues[0]>25&&polyClassValues[1]>25)//mixed
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(13));
		  else{
		    if(verbose_opt[0]){
		      std::cout << "No valid value in polyClassValues..." << std::endl;
		      for(int index=0;index<polyClassValues.size();++index){
			double theValue=polyClassValues[index];
			std::cout << theValue << " ";
		      }
		      std::cout << std::endl;
		    }
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(20));
		  }
		}
		else if(ruleMap[rule_opt[0]]==rule::mode){
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
		  if(polygon_opt[0])
		    writePolygonFeature->SetField(label_opt[0].c_str(),maxClass);
		  else
		    writeCentroidFeature->SetField(label_opt[0].c_str(),maxClass);
		}
	      }
	      if(polygon_opt[0]){
		if(verbose_opt[0]>1)
		  std::cout << "creating polygon feature" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		else{
		  if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writePolygonFeature );
		++ntotalvalid;
		++ntotalvalidLayer;
	      }
	      else{
		if(verbose_opt[0]>1)
		  std::cout << "creating point feature in centroid" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		else{
		  //test
		  assert(validFeature);
		  if(writeLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writeCentroidFeature );
		++ntotalvalid;
		++ntotalvalidLayer;
	      }
	    }
	  }
	  else if(wkbFlatten(poGeometry->getGeometryType()) == wkbMultiPolygon){//todo: try to use virtual OGRGeometry instead of OGRMultiPolygon and OGRPolygon
	    OGRMultiPolygon readPolygon = *((OGRMultiPolygon *) poGeometry);
	    OGRPolygon writePolygon;
	    OGRLinearRing writeRing;
	    OGRPoint writeCentroidPoint;
	    OGRFeature *writePolygonFeature;
	    OGRFeature *writeCentroidFeature;

	    readPolygon.closeRings();

	    if(verbose_opt[0]>1)
	      std::cout << "get centroid point from polygon" << std::endl;
	    readPolygon.Centroid(&writeCentroidPoint);

	    double ulx,uly,lrx,lry;
	    double uli,ulj,lri,lrj;
	    if((polygon_opt[0]&&ruleMap[rule_opt[0]]==rule::point)||(ruleMap[rule_opt[0]]==rule::centroid)){
	      ulx=writeCentroidPoint.getX();
	      uly=writeCentroidPoint.getY();
	      lrx=ulx;
	      lry=uly;
	    }
	    else{
	      //get envelope
	      if(verbose_opt[0])
		std::cout << "reading envelope for polygon " << ifeature << std::endl;
	      OGREnvelope* psEnvelope=new OGREnvelope();
	      readPolygon.getEnvelope(psEnvelope);
	      ulx=psEnvelope->MinX;
	      uly=psEnvelope->MaxY;
	      lrx=psEnvelope->MaxX;
	      lry=psEnvelope->MinY;
	      delete psEnvelope;
	    }
	    // if(geo_opt[0]){
	      imgReader.geo2image(ulx,uly,uli,ulj);
	      imgReader.geo2image(lrx,lry,lri,lrj);
	    // }
	    // else{
	    //   uli=ulx;
	    //   ulj=uly;
	    //   lri=lrx;
	    //   lrj=lry;
	    // }
	    //nearest neighbour
	    ulj=static_cast<int>(ulj);
	    uli=static_cast<int>(uli);
	    lrj=static_cast<int>(lrj);
	    lri=static_cast<int>(lri);
	    //iterate through all pixels
	    if(verbose_opt[0]>1)
	      std::cout << "bounding box for multipologon feature " << ifeature << ": " << uli << " " << ulj << " " << lri << " " << lrj << std::endl;

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
	    // if(uli<0||lri>=imgReader.nrOfCol()||ulj<0||lrj>=imgReader.nrOfRow())
	    //   continue;

	    int nPointPolygon=0;
	    if(polygon_opt[0]){
	      if(writeTest)
		writePolygonFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    else if(ruleMap[rule_opt[0]]!=rule::point){
	      if(writeTest)
		writeCentroidFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writeCentroidFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    // vector<double> polyValues;
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
	    vector< Vector2d<double> > readValues(nband);
	    for(int iband=0;iband<nband;++iband){
	      int theBand=(band_opt.size()) ? band_opt[iband] : iband;
	      imgReader.readDataBlock(readValues[iband],GDT_Float64,uli,lri,ulj,lrj,theBand);
	    }

	    OGRPoint thePoint;
	    for(int j=ulj;j<=lrj;++j){
	      for(int i=uli;i<=lri;++i){
		//check if within raster image
		if(i<0||i>=imgReader.nrOfCol())
		  continue;
		if(j<0||j>=imgReader.nrOfRow())
		  continue;
		//check if point is on surface
		double theX=0;
		double theY=0;
		imgReader.image2geo(i,j,theX,theY);
		thePoint.setX(theX);
		thePoint.setY(theY);

		if(ruleMap[rule_opt[0]]!=rule::centroid&&!readPolygon.Contains(&thePoint))
		  continue;

		bool valid=true;

		if(srcnodata_opt.size()){
		  for(int vband=0;vband<bndnodata_opt.size();++vband){
		    double value=((readValues[bndnodata_opt[vband]])[j-ulj])[i-uli];
		    if(value==srcnodata_opt[vband]){
		      valid=false;
		      break;
		    }
		  }
		}

		if(!valid)
		  continue;
		else
		  validFeature=true;
		
		writeRing.addPoint(&thePoint);
		// if(writeRing.isPointOnRingBoundary(&thePoint))
		//    writeRing.addPoint(&thePoint);
		if(verbose_opt[0]>1)
		    std::cout << "point is on surface:" << thePoint.getX() << "," << thePoint.getY() << std::endl;
		  ++nPointPolygon;

		  if(polythreshold_opt.size())
		    if(nPointPolygon>polythreshold_opt[0])
		      continue;
		  // throw(nPointPolygon);
		  OGRFeature *writePointFeature;
		  if(!polygon_opt[0]){
		    //create feature
		    if(ruleMap[rule_opt[0]]==rule::point){//do not create in case of mean, stdev, median or sum (only create point at centroid)
		      if(writeTest)
			writePointFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
		      else
			writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
		      if(verbose_opt[0]>1)
			std::cout << "copying fields from polygons " << std::endl;
		      if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
			cerr << "writing feature failed" << std::endl;
		      if(verbose_opt[0]>1)
			std::cout << "set geometry as point " << std::endl;
		      //test
		      writePointFeature->SetGeometry(&thePoint);
		      assert(wkbFlatten(writePointFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
		      // OGRGeometry *updateGeometry;
		      // updateGeometry = writePointFeature->GetGeometryRef();
		      // OGRPoint *poPoint = (OGRPoint *) updateGeometry;
		      if(verbose_opt[0]>1)
			std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields" << std::endl;
		    }
		  }
		  if(class_opt.size()){
		    short value=((readValues[0])[j-ulj])[i-uli];
		    for(int iclass=0;iclass<class_opt.size();++iclass){
		      if(value==class_opt[iclass])
			polyClassValues[iclass]+=1;
		    }
		  }
		  else{
		    for(int iband=0;iband<nband;++iband){
		      double value=((readValues[iband])[j-ulj])[i-uli];
		      if(verbose_opt[0]>1)
			std::cout << ": " << value << std::endl;
		      if(polygon_opt[0]||ruleMap[rule_opt[0]]!=rule::point)
			polyValues[iband].push_back(value);
		      else{
			if(verbose_opt[0]>1)
			  std::cout << "set field " << fieldname_opt[iband] << " to " << value << std::endl;
			switch( fieldType ){
			case OFTInteger:
			case OFTReal:
			  writePointFeature->SetField(fieldname_opt[iband].c_str(),value);
			  break;
			case OFTString:
			  writePointFeature->SetField(fieldname_opt[iband].c_str(),type2string<double>(value).c_str());
			  break;
			default://not supported
			  assert(0);
			  break;
			}
		      }//else
		    }//iband
		  }//else (class_opt.size())
		  if(!polygon_opt[0]){
		    if(ruleMap[rule_opt[0]]==rule::point){//do not create in case of mean, stdev or median value (only at centroid)
		      //write feature
		      if(verbose_opt[0]>1)
			std::cout << "creating point feature" << std::endl;
		      if(writeTest){
			if(writeTestLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			  std::string errorString="Failed to create feature in ogr vector dataset";
			  throw(errorString);
			}
		      }
		      else{
			if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			  std::string errorString="Failed to create feature in ogr vector dataset";
			  throw(errorString);
			}
		      }
		      //destroy feature
		      OGRFeature::DestroyFeature( writePointFeature );
		    }
		  }
		  // ++isample;
		  ++ntotalvalid;
		  ++ntotalvalidLayer;
	      }
	    }
	    if(!validFeature)
	      continue;
	    if(polygon_opt[0]||ruleMap[rule_opt[0]]!=rule::point){
	      //do not create if no points found within polygon
	      if(!nPointPolygon)
		continue;
	      //add ring to polygon
	      if(polygon_opt[0]){
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
		//write geometry of writePolygon
		//test
		//writePolygon and readFeature are from geometry type wkbMultiPolygon
		// writePolygonFeature->SetGeometry(&writePolygon);
		if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		assert(writePolygonFeature->GetGeometryRef()->getGeometryType()==wkbMultiPolygon);
		if(verbose_opt[0]>1)
		  std::cout << "copying new fields write polygon " << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
		//write polygon feature
	      }
	      else{//write band information of polygon to centroid point
		//create feature
		if(verbose_opt[0]>1)
		  std::cout << "copying fields from polygons " << std::endl;
		if(writeCentroidFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		writeCentroidFeature->SetGeometry(&writeCentroidPoint);
		assert(wkbFlatten(writeCentroidFeature->GetGeometryRef()->getGeometryType()) == wkbPoint);
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writeCentroidFeature->GetFieldCount() << " fields" << std::endl;
	      }
	      if(class_opt.empty()){
		if(ruleMap[rule_opt[0]]==rule::point){//value at centroid of polygon
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  for(int index=0;index<polyValues.size();++index){
		    //test
		    assert(polyValues[index].size()==1);
		    double theValue=polyValues[index].back();
		    if(verbose_opt[0])
		      std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		    int theBand=(band_opt.size()) ? band_opt[index] : index;
		    try{
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
			break;
		      case OFTString:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			break;
		      default://not supported
			std::string errorString="field type not supported";
			throw(errorString);
			break;
		      }
		    }
		    catch(std::string e){
		      std::cout << e << std::endl;
		      exit(1);
		    }
		  }
		}
		else{//ruleMap[rule_opt[0]] is not rule::point
		  double theValue=0;
		  for(int index=0;index<polyValues.size();++index){
		    try{
		      if(ruleMap[rule_opt[0]]==rule::mean)
			theValue=stat.mean(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::stdev)
			theValue=sqrt(stat.var(polyValues[index]));
		      else if(ruleMap[rule_opt[0]]==rule::median)
			theValue=stat.median(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::percentile)
			theValue=stat.percentile(polyValues[index],polyValues[index].begin(),polyValues[index].end(),percentile_opt[0]);
		      else if(ruleMap[rule_opt[0]]==rule::sum)
			theValue=stat.sum(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::max)
			theValue=stat.mymax(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::min)
			theValue=stat.mymin(polyValues[index]);
		      else if(ruleMap[rule_opt[0]]==rule::centroid){
			if(verbose_opt[0])
			  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
			assert(nPointPolygon<=1);
			assert(nPointPolygon==polyValues[index].size());
			theValue=polyValues[index].back();
		      }
		      else{
			std::string errorString="rule not supported";
			throw(errorString);
		      }
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
		      case OFTReal:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
			break;
		      case OFTString:
			if(polygon_opt[0])
			  writePolygonFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			else
			  writeCentroidFeature->SetField(fieldname_opt[index].c_str(),type2string<double>(theValue).c_str());
			break;
		      default://not supported
			std::string errorString="field type not supported";
			throw(errorString);
			break;
		      }
		    }
		    catch(std::string e){
		      std::cout << e << std::endl;
		      exit(1);
		    }
		  }
		}
	      }
	      else{//class_opt is set
		if(ruleMap[rule_opt[0]]==rule::proportion||ruleMap[rule_opt[0]]==rule::count){
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  if(ruleMap[rule_opt[0]]==rule::proportion)
		    stat.normalize_pct(polyClassValues);
		  for(int index=0;index<polyClassValues.size();++index){
		    double theValue=polyClassValues[index];
		    ostringstream fs;
		    fs << class_opt[index];
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		  }
		}
		else if(ruleMap[rule_opt[0]]==rule::custom){
		  assert(polygon_opt[0]);//not implemented for points
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  stat.normalize_pct(polyClassValues);
		  assert(polyClassValues.size()==2);//11:broadleaved, 12:coniferous
		  if(polyClassValues[0]>=75)//broadleaved
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(11));
		  else if(polyClassValues[1]>=75)//coniferous
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(12));
		  else if(polyClassValues[0]>25&&polyClassValues[1]>25)//mixed
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(13));
		  else{
		    if(verbose_opt[0]){
		      std::cout << "No valid value in polyClassValues..." << std::endl;
		      for(int index=0;index<polyClassValues.size();++index){
			double theValue=polyClassValues[index];
			std::cout << theValue << " ";
		      }
		      std::cout << std::endl;
		    }
		    writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(20));
		  }
		}
		else if(ruleMap[rule_opt[0]]==rule::mode){
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
		  if(polygon_opt[0])
		    writePolygonFeature->SetField(label_opt[0].c_str(),maxClass);
		  else
		    writeCentroidFeature->SetField(label_opt[0].c_str(),maxClass);
		}
	      }

	      if(polygon_opt[0]){
		if(verbose_opt[0]>1)
		  std::cout << "creating polygon feature" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		else{
		  if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writePolygonFeature );
		++ntotalvalid;
		++ntotalvalidLayer;
	      }
	      else{
		if(verbose_opt[0]>1)
		  std::cout << "creating point feature in centroid" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		else{
		  //test
		  assert(validFeature);
		  if(writeLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in ogr vector dataset";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writeCentroidFeature );
		++ntotalvalid;
		++ntotalvalidLayer;
	      }
	    }
	  }
	  else{
	    std::string test;
	    test=poGeometry->getGeometryName();
	    ostringstream oss;
	    oss << "geometry " << test << " not supported";
	    throw(oss.str());
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
      }//end of getNextFeature
      // if(rbox_opt[0]>0||cbox_opt[0]>0)
      //   boxWriter.close();
      progress=1.0;
      pfnProgress(progress,pszMessage,pProgressArg);
      ++ilayerWrite;
    }//for ilayer
    sampleReaderOgr.close();
    ogrWriter.close();
    if(test_opt.size())
      ogrTestWriter.close();
  }//else (vector)
  progress=1.0;
  pfnProgress(progress,pszMessage,pProgressArg);
  imgReader.close();
}
  

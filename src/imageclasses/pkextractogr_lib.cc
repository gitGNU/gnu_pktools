/**********************************************************************
pkextractogr_lib.cc: extract pixel values from raster image from a vector sample
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
#include <sstream>
#include <iostream>
#include <algorithm>
#include <ctime>
#include <vector>
#include <memory>
#include "imageclasses/ImgRaster.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

namespace rule{
  enum RULE_TYPE {point=0, mean=1, proportion=2, custom=3, min=4, max=5, mode=6, centroid=7, sum=8, median=9, stdev=10, percentile=11, count=12, allpoints=13};
}



/**
 * @param app application specific option arguments
 * 
 * @return CE_None if success, CE_Failure if failure
 */
CPLErr ImgRaster::extractOgr(const AppFactory& app){
  // Optionpk<string> image_opt("i", "input", "Raster input dataset containing band information");
  Optionpk<string> sample_opt("s", "sample", "OGR vector dataset with features to be extracted from input data. Output will contain features with input band information included. Sample image can also be GDAL raster dataset.");
  Optionpk<string> layer_opt("ln", "ln", "Layer name(s) in sample (leave empty to select all)");
  Optionpk<unsigned int> random_opt("rand", "random", "Create simple random sample of points. Provide number of points to generate");
  Optionpk<double> grid_opt("grid", "grid", "Create systematic grid of points. Provide cell grid size (in projected units, e.g,. m)");
  Optionpk<string> output_opt("o", "output", "Output sample dataset");
  Optionpk<int> class_opt("c", "class", "Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset. Make sure to set classes if rule is set to mode, proportion or count");
  Optionpk<float> threshold_opt("t", "threshold", "Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold per vector sample layer. If using raster land cover maps as a sample dataset, you can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es)", 100);
  Optionpk<double> percentile_opt("perc","perc","Percentile value(s) used for rule percentile",95);
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
  // Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",1000,1);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose mode if > 0", 0,2);

  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  bndnodata_opt.setHide(1);
  srcnodata_opt.setHide(1);
  polythreshold_opt.setHide(1);
  percentile_opt.setHide(1);
  buffer_opt.setHide(1);
  disc_opt.setHide(1);
  // memory_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=sample_opt.retrieveOption(app.getArgc(),app.getArgv());
    layer_opt.retrieveOption(app.getArgc(),app.getArgv());
    random_opt.retrieveOption(app.getArgc(),app.getArgv());
    grid_opt.retrieveOption(app.getArgc(),app.getArgv());
    output_opt.retrieveOption(app.getArgc(),app.getArgv());
    class_opt.retrieveOption(app.getArgc(),app.getArgv());
    threshold_opt.retrieveOption(app.getArgc(),app.getArgv());
    percentile_opt.retrieveOption(app.getArgc(),app.getArgv());
    ogrformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    ftype_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    bstart_opt.retrieveOption(app.getArgc(),app.getArgv());
    bend_opt.retrieveOption(app.getArgc(),app.getArgv());
    rule_opt.retrieveOption(app.getArgc(),app.getArgv());
    bndnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    srcnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    polythreshold_opt.retrieveOption(app.getArgc(),app.getArgv());
    buffer_opt.retrieveOption(app.getArgc(),app.getArgv());
    disc_opt.retrieveOption(app.getArgc(),app.getArgv());
    // memory_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
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

    statfactory::StatFactory stat;
    if(srcnodata_opt.size()){
      while(srcnodata_opt.size()<bndnodata_opt.size())
        srcnodata_opt.push_back(srcnodata_opt[0]);
      stat.setNoDataValues(srcnodata_opt);
    }
    // if(verbose_opt[0])
    //   std::cout << class_opt << std::endl;
    Vector2d<unsigned int> posdata;

    // ImgRaster imgReader;
    // ImgRaster imgReader;
    // if(image_opt.empty()){
    //   std::cerr << "No image dataset provided (use option -i). Use --help for help information";
    //     exit(1);
    // }
    if(output_opt.empty()){
      std::cerr << "No output dataset provided (use option -o). Use --help for help information";
      return(CE_Failure);
    }
    // try{
    //   imgReader.open(image_opt[0],memory_opt[0]);
    // }
    // catch(std::string errorstring){
    //   std::cout << errorstring << std::endl;
    //   exit(1);
    // }

    //check if rule contains allpoints
    if(find(rule_opt.begin(),rule_opt.end(),"allpoints")!=rule_opt.end()){
      //allpoints should be the only rule
      rule_opt.clear();
      rule_opt.push_back("allpoints");
    }

    //convert start and end band options to vector of band indexes
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
    int nband=(band_opt.size()) ? band_opt.size() : this->nrOfBand();
    if(class_opt.size()){
      if(nband>1){
        cerr << "Warning: using only first band or multiband image" << endl;
        nband=1;
        band_opt.clear();
        band_opt.push_back(0);
      }
    }

    if(verbose_opt[0]>1)
      std::cout << "Number of bands in input image: " << this->nrOfBand() << std::endl;

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
      return(CE_Failure);
      break;
    }
  
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    srand(time(NULL));

    bool sampleIsRaster=false;

    ImgReaderOgr sampleReaderOgr;
    ImgWriterOgr sampleWriterOgr;

    if(sample_opt.size()){
      try{
        sampleReaderOgr.open(sample_opt[0]);
      }
      catch(string errorString){
        //todo: sampleIsRaster will not work from GDAL 2.0!!?? (unification of driver for raster and vector datasets)
        sampleIsRaster=true;
      }
    }
    else{
      sampleWriterOgr.open("/vsimem/virtual",ogrformat_opt[0]);
      char     **papszOptions=NULL;

      if(random_opt.size()){
        //create simple random sampling within boundary
        double ulx,uly,lrx,lry;
        this->getBoundingBox(ulx,uly,lrx,lry);
        if(random_opt[0]>0)
          sampleWriterOgr.createLayer("points", this->getProjection(), wkbPoint, papszOptions);
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
        this->getBoundingBox(ulx,uly,lrx,lry);
        if(uly-grid_opt[0]/2<lry&&ulx+grid_opt[0]/2>lrx){
          string errorString="Error: grid distance too large";
          throw(errorString);
        }
        else if(grid_opt[0]>0)
          sampleWriterOgr.createLayer("points", this->getProjection(), wkbPoint, papszOptions);
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
        string errorString="Error: no sample dataset provided (use option -s). Use --help for help information";
        throw(errorString);
      }
      sampleWriterOgr.close();
      sampleReaderOgr.open("/vsimem/virtual");
    }

    if(sampleIsRaster){
      cerr << "Error: sample must be vector dataset in OGR format";
      return(CE_Failure);
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

      sampleReaderOgr.getExtent(vectords_ulx,vectords_uly,vectords_lrx,vectords_lry);
      bool hasCoverage=((vectords_ulx < this->getLrx())&&(vectords_lrx > this->getUlx())&&(vectords_lry < this->getUly())&&(vectords_uly > this->getLry()));
      if(!hasCoverage){
        ostringstream ess;
        ess << "No coverage for any layer in " << sample_opt[0] << endl;
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
            this->getMinMax(minValue,maxValue,theBand);
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
        bool hasCoverage=((layer_ulx < this->getLrx())&&(layer_lrx > this->getUlx())&&(layer_lry < this->getUly())&&(layer_uly > this->getLry()));
        if(!hasCoverage)
          continue;

        //align bounding box to input image
        layer_ulx-=fmod(layer_ulx-this->getUlx(),this->getDeltaX());
        layer_lrx+=fmod(this->getLrx()-layer_lrx,this->getDeltaX());
        layer_uly+=fmod(this->getUly()-layer_uly,this->getDeltaY());
        layer_lry-=fmod(layer_lry-this->getLry(),this->getDeltaY());

        //do not read outside input image
        if(layer_ulx<this->getUlx())
          layer_ulx=this->getUlx();
        if(layer_lrx>this->getLrx())
          layer_lrx=this->getLrx();
        if(layer_uly>this->getUly())
          layer_uly=this->getUly();
        if(layer_lry<this->getLry())
          layer_lry=this->getLry();

        //read entire block for coverage in memory
        //todo: use different data types
        vector< Vector2d<float> > readValuesReal(nband);
        vector< Vector2d<int> > readValuesInt(nband);

        double layer_uli;
        double layer_ulj;
        double layer_lri;
        double layer_lrj;
        this->geo2image(layer_ulx,layer_uly,layer_uli,layer_ulj);
        this->geo2image(layer_lrx,layer_lry,layer_lri,layer_lrj);

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
        layer_lri=(layer_lri>=this->nrOfCol())? this->nrOfCol()-1 : static_cast<int>(layer_lri);
        layer_lrj=(layer_lrj>=this->nrOfRow())? this->nrOfRow()-1 : static_cast<int>(layer_lrj);

        for(int iband=0;iband<nband;++iband){
          int theBand=(band_opt.size()) ? band_opt[iband] : iband;
          if(theBand<0){
            string errorString="Error: illegal band (must be positive and starting from 0)";
            throw(errorString);
          }
          if(theBand>=this->nrOfBand()){
            string errorString="Error: illegal band (must be lower than number of bands in input raster dataset)";
            throw(errorString);
          }
          if(verbose_opt[0])
            cout << "reading image band " << theBand << " block rows " << layer_ulj << "-" << layer_lrj << ", cols " << layer_uli << "-" << layer_lri << endl;
          switch( fieldType ){
          case OFTInteger:
            this->readDataBlock(readValuesInt[iband],layer_uli,layer_lri,layer_ulj,layer_lrj,theBand);
            break;
          case OFTReal:
          default:
            this->readDataBlock(readValuesReal[iband],layer_uli,layer_lri,layer_ulj,layer_lrj,theBand);
            break;
          }
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
          writeLayer=ogrWriter.createLayer(readLayer->GetName(), this->getProjection(), wkbPolygon, papszOptions);
        }
        else{
          if(verbose_opt[0])
            std::cout << "create points in layer " << readLayer->GetName() << std::endl;
          char **papszOptions=NULL;

          writeLayer=ogrWriter.createLayer(readLayer->GetName(), this->getProjection(), wkbPoint, papszOptions);
        }
        if(verbose_opt[0])
          std::cout << "copy fields from layer " << ilayer << std::flush << std::endl;
        ogrWriter.copyFields(sampleReaderOgr,ilayer,ilayerWrite);

        for(int irule=0;irule<rule_opt.size();++irule){
          for(int iband=0;iband<nband;++iband){
            int theBand=(band_opt.size()) ? band_opt[iband] : iband;
            ostringstream fs;
            if(rule_opt.size()>1||nband==1)
              fs << rule_opt[irule];
            if(nband>1)
              fs << "b" << theBand;
            switch(ruleMap[rule_opt[irule]]){
            case(rule::proportion):
            case(rule::count):{//count for each class
              for(int iclass=0;iclass<class_opt.size();++iclass){
                ostringstream fsclass;
                fsclass << fs.str() << "class" << class_opt[iclass];
                ogrWriter.createField(fsclass.str(),fieldType,ilayerWrite);
              }
              break;
            }
            case(rule::percentile):{//for each percentile
              for(int iperc=0;iperc<percentile_opt.size();++iperc){
                ostringstream fsperc;
                fsperc << fs.str() << percentile_opt[iperc];
                ogrWriter.createField(fsperc.str(),fieldType,ilayerWrite);
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
              this->geo2image(readPoint.getX(),readPoint.getY(),i_centre,j_centre);
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
              if(static_cast<int>(ulj)<0||static_cast<int>(ulj)>=this->nrOfRow())
                continue;
              //check if j is out of bounds
              if(static_cast<int>(uli)<0||static_cast<int>(lri)>=this->nrOfCol())
                continue;
            
              OGRPoint ulPoint,urPoint,llPoint,lrPoint;
              OGRPolygon writePolygon;
              OGRPoint writePoint;
              OGRLinearRing writeRing;
              OGRFeature *writePolygonFeature;

              int nPointPolygon=0;
              if(createPolygon){
                if(disc_opt[0]){
                  double radius=buffer_opt[0]*sqrt(this->getDeltaX()*this->getDeltaY());
                  unsigned short nstep = 25;
                  for(int i=0;i<nstep;++i){
                    OGRPoint aPoint;
                    aPoint.setX(readPoint.getX()+this->getDeltaX()/2.0+radius*cos(2*PI*i/nstep));
                    aPoint.setY(readPoint.getY()-this->getDeltaY()/2.0+radius*sin(2*PI*i/nstep));
                    writeRing.addPoint(&aPoint);
                  }
                  writePolygon.addRing(&writeRing);
                  writePolygon.closeRings();
                }
                else{
                  double ulx,uly,lrx,lry;
                  this->image2geo(uli,ulj,ulx,uly);
                  this->image2geo(lri,lrj,lrx,lry);
                  ulPoint.setX(ulx-this->getDeltaX()/2.0);
                  ulPoint.setY(uly+this->getDeltaY()/2.0);
                  lrPoint.setX(lrx+this->getDeltaX()/2.0);
                  lrPoint.setY(lry-this->getDeltaY()/2.0);
                  urPoint.setX(lrx+this->getDeltaX()/2.0);
                  urPoint.setY(uly+this->getDeltaY()/2.0);
                  llPoint.setX(ulx-this->getDeltaX()/2.0);
                  llPoint.setY(lry-this->getDeltaY()/2.0);

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
                  this->geo2image(readPoint.getX(),readPoint.getY(),i,j);
                  int indexJ=static_cast<int>(j-layer_ulj);
                  int indexI=static_cast<int>(i-layer_uli);
                  bool valid=true;
                  valid=valid&&(indexJ>=0);
                  valid=valid&&(indexJ<this->nrOfRow());
                  valid=valid&&(indexI>=0);
                  valid=valid&&(indexI<this->nrOfCol());

                  if(valid){
                    if(srcnodata_opt.empty())
                      validFeature=true;
                    else{
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
                  }
                  if(valid){
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
                  this->geo2image(readPoint.getX(),readPoint.getY(),i,j);
                  int indexJ=static_cast<int>(j-layer_ulj);
                  int indexI=static_cast<int>(i-layer_uli);
                  bool valid=true;
                  valid=valid&&(indexJ>=0);
                  valid=valid&&(indexJ<this->nrOfRow());
                  valid=valid&&(indexI>=0);
                  valid=valid&&(indexI<this->nrOfCol());

                  if(valid){
                    if(srcnodata_opt.empty())
                      validFeature=true;
                    else{
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
                    if(i<0||i>=this->nrOfCol())
                      continue;
                    if(j<0||j>=this->nrOfRow())
                      continue;
                    int indexJ=j-layer_ulj;
                    int indexI=i-layer_uli;
                    if(indexJ<0)
                      indexJ=0;
                    if(indexI<0)
                      indexI=0;
                    if(indexJ>=this->nrOfRow())
                      indexJ=this->nrOfRow()-1;
                    if(indexI>=this->nrOfCol())
                      indexI=this->nrOfCol()-1;

                    double theX=0;
                    double theY=0;
                    this->image2geo(i,j,theX,theY);
                    thePoint.setX(theX);
                    thePoint.setY(theY);
                    if(disc_opt[0]&&buffer_opt[0]>0){
                      double radius=buffer_opt[0]*sqrt(this->getDeltaX()*this->getDeltaY());
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
                            return(CE_Failure);
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
                      vector<double> theValue;
                      vector<string> fieldname;
                      ostringstream fs;
                      if(rule_opt.size()>1||nband==1)
                        fs << rule_opt[irule];
                      if(nband>1)
                        fs << "b" << theBand;
                      switch(ruleMap[rule_opt[irule]]){
                      case(rule::proportion)://deliberate fall through
                        stat.normalize_pct(polyClassValues);
                      case(rule::count):{//count for each class
                        for(int index=0;index<polyClassValues.size();++index){
                          theValue.push_back(polyClassValues[index]);
                          ostringstream fsclass;
                          fsclass << fs.str() << "class" << class_opt[index];
                          fieldname.push_back(fsclass.str());
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
                        theValue.push_back(maxClass);
                        fieldname.push_back(fs.str());
                        break;
                      }
                      case(rule::mean):
                        theValue.push_back(stat.mean(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::median):
                        theValue.push_back(stat.median(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::stdev):
                        theValue.push_back(sqrt(stat.var(polyValues[iband])));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::percentile):{
                        for(int iperc=0;iperc<percentile_opt.size();++iperc){
                          theValue.push_back(stat.percentile(polyValues[iband],polyValues[iband].begin(),polyValues[iband].end(),percentile_opt[iperc]));
                          ostringstream fsperc;
                          fsperc << fs.str() << percentile_opt[iperc];
                          fieldname.push_back(fsperc.str());
                        }
                        break;
                      }
                      case(rule::sum):
                        theValue.push_back(stat.sum(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::max):
                        theValue.push_back(stat.mymax(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::min):
                        theValue.push_back(stat.mymin(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::point):
                      case(rule::centroid):
                        theValue.push_back(polyValues[iband].back());
		      fieldname.push_back(fs.str());
                      break;
                      default://not supported
                        break;
                      }
                      for(int ivalue=0;ivalue<theValue.size();++ivalue){
                        switch( fieldType ){
                        case OFTInteger:
                          writePolygonFeature->SetField(fieldname[ivalue].c_str(),static_cast<int>(theValue[ivalue]));
                          break;
                        case OFTReal:
                          writePolygonFeature->SetField(fieldname[ivalue].c_str(),theValue[ivalue]);
                          break;
                        case OFTString:
                          writePolygonFeature->SetField(fieldname[ivalue].c_str(),type2string<double>(theValue[ivalue]).c_str());
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
              }
              //test
              if(createPolygon&&validFeature){
                // if(createPolygon){
                //write polygon feature
                //todo: create only in case of valid feature
                if(verbose_opt[0]>1)
                  std::cout << "creating polygon feature (1)" << std::endl;
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
              if(!this->covers(ulx,uly,lrx,lry))
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
                  std::cout << "write polygon feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
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
                this->geo2image(readPoint.getX(),readPoint.getY(),i,j);
                int indexJ=static_cast<int>(j-layer_ulj);
                int indexI=static_cast<int>(i-layer_uli);
                bool valid=true;
                valid=valid&&(indexJ>=0);
                valid=valid&&(indexJ<this->nrOfRow());
                valid=valid&&(indexI>=0);
                valid=valid&&(indexI<this->nrOfCol());
                if(valid){
                  if(srcnodata_opt.empty())
                    validFeature=true;
                  else{
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
                }
                if(valid){
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
                this->geo2image(readPoint.getX(),readPoint.getY(),i,j);
                int indexJ=static_cast<int>(j-layer_ulj);
                int indexI=static_cast<int>(i-layer_uli);
                bool valid=true;
                valid=valid&&(indexJ>=0);
                valid=valid&&(indexJ<this->nrOfRow());
                valid=valid&&(indexI>=0);
                valid=valid&&(indexI<this->nrOfCol());
                if(valid){
                  if(srcnodata_opt.empty())
                    validFeature=true;
                  else{
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
                this->geo2image(ulx,uly,uli,ulj);
                this->geo2image(lrx,lry,lri,lrj);
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
                if(uli>=this->nrOfCol())
                  uli=this->nrOfCol()-1;
                if(lri>=this->nrOfCol())
                  lri=this->nrOfCol()-1;
                if(ulj<0)
                  ulj=0;
                if(lrj<0)
                  lrj=0;
                if(ulj>=this->nrOfRow())
                  ulj=this->nrOfRow()-1;
                if(lrj>=this->nrOfRow())
                  lrj=this->nrOfRow()-1;

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
                    if(i<0||i>=this->nrOfCol())
                      continue;
                    if(j<0||j>=this->nrOfRow())
                      continue;
                    int indexJ=j-layer_ulj;
                    int indexI=i-layer_uli;

                    double theX=0;
                    double theY=0;
                    this->image2geo(i,j,theX,theY);
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
                        std::cout << "write point feature has " << writePointFeature->GetFieldCount() << " fields:" << std::endl;
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
                            return(CE_Failure);
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
                      //todo: only if valid feature?
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
                      vector<double> theValue;
                      vector<string> fieldname;
                      ostringstream fs;
                      if(rule_opt.size()>1||nband==1)
                        fs << rule_opt[irule];
                      if(nband>1)
                        fs << "b" << theBand;
                      switch(ruleMap[rule_opt[irule]]){
                      case(rule::proportion):
                        stat.normalize_pct(polyClassValues);
                      case(rule::count):{//count for each class
                        for(int index=0;index<polyClassValues.size();++index){
                          theValue.push_back(polyClassValues[index]);
                          ostringstream fsclass;
                          fsclass << fs.str() << "class" << class_opt[index];
                          fieldname.push_back(fsclass.str());
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
                        theValue.push_back(maxClass);
                        fieldname.push_back(fs.str());
                        break;
                      }
                      case(rule::mean):
                        theValue.push_back(stat.mean(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::median):
                        theValue.push_back(stat.median(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::stdev):
                        theValue.push_back(sqrt(stat.var(polyValues[iband])));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::percentile):{
                        for(int iperc=0;iperc<percentile_opt.size();++iperc){
                          theValue.push_back(stat.percentile(polyValues[iband],polyValues[iband].begin(),polyValues[iband].end(),percentile_opt[iperc]));
                          ostringstream fsperc;
                          fsperc << fs.str() << percentile_opt[iperc];
                          fieldname.push_back(fsperc.str());
                        }
                        break;
                      }
                      case(rule::sum):
                        theValue.push_back(stat.sum(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::max):
                        theValue.push_back(stat.mymax(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::min):
                        theValue.push_back(stat.mymin(polyValues[iband]));
                        fieldname.push_back(fs.str());
                        break;
                      case(rule::centroid):
                      case(rule::point):
                        theValue.push_back(polyValues[iband].back());
		      fieldname.push_back(fs.str());
                      break;
                      default://not supported
                        break;
                      }
                      for(int ivalue=0;ivalue<theValue.size();++ivalue){
                        switch( fieldType ){
                        case OFTInteger:
                          writePolygonFeature->SetField(fieldname[ivalue].c_str(),static_cast<int>(theValue[ivalue]));
                          break;
                        case OFTReal:
                          writePolygonFeature->SetField(fieldname[ivalue].c_str(),theValue[ivalue]);
                          break;
                        case OFTString:
                          writePolygonFeature->SetField(fieldname[ivalue].c_str(),type2string<double>(theValue[ivalue]).c_str());
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
              }
              if(createPolygon&&validFeature){
                //todo: only create if valid feature?
                //write polygon feature
                if(verbose_opt[0]>1)
                  std::cout << "creating polygon feature (2)" << std::endl;
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
    // this->close();
    return(CE_None);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}


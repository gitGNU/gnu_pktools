/**********************************************************************
pkextract.cc: extract pixel values from raster image from a (vector or raster) sample
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
#include <math.h>
#include <stdlib.h>
#ifdef WIN32
#define random rand
#define srandom srand
#endif
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

namespace rule{
  enum RULE_TYPE {point=0, mean=1, proportion=2, custom=3, minimum=4, maximum=5, maxvote=6, centroid=7, sum=8};
}

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> image_opt("i", "image", "Input image file");
  Optionpk<string> sample_opt("s", "sample", "Input sample vector file or class file (e.g. Corine CLC) if class option is set");
  Optionpk<string> mask_opt("m", "mask", "Mask image file");
  Optionpk<int> msknodata_opt("msknodata", "msknodata", "Mask value where image is invalid. If a single mask is used, more nodata values can be set. If more masks are used, use one value for each mask.", 1);
  Optionpk<int> class_opt("c", "class", "Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample file");
  Optionpk<string> output_opt("o", "output", "Output sample file (image file)");
  Optionpk<string> ogrformat_opt("f", "f", "Output sample file format","ESRI Shapefile");
  Optionpk<string> test_opt("test", "test", "Test sample file (use this option in combination with threshold<100 to create a training (output) and test set");
  // Optionpk<string> bufferOutput_opt("bu", "bu", "Buffer output shape file");
  Optionpk<short> geo_opt("g", "geo", "geo coordinates", 1);
  Optionpk<short> down_opt("down", "down", "down sampling factor. Can be used to create grid points", 1);
  Optionpk<float> threshold_opt("t", "threshold", "threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold for vector sample files. If using raster land cover maps as a sample file, you can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es)", 100);
  Optionpk<double> min_opt("min", "min", "minimum number of samples to select (0)", 0);
  Optionpk<short> boundary_opt("bo", "boundary", "boundary for selecting the sample", 1);
  // Optionpk<short> rbox_opt("rb", "rbox", "rectangular boundary box (total width in m) to draw around the selected pixel. Can not combined with class option. Use multiple rbox options for multiple boundary boxes. Use value 0 for no box)", 0);
  // Optionpk<short> cbox_opt("cbox", "cbox", "circular boundary (diameter in m) to draw around the selected pixel. Can not combined with class option. Use multiple cbox options for multiple boundary boxes. Use value 0 for no box)", 0);
  Optionpk<short> disc_opt("circ", "circular", "circular disc kernel boundary", 0);
  Optionpk<string> ftype_opt("ft", "ftype", "Field type (only Real or Integer)", "Real");
  Optionpk<string> ltype_opt("lt", "ltype", "Label type: In16 or String", "Integer");
  Optionpk<string> fieldname_opt("bn", "bname", "Attribute field name of extracted raster band", "B");
  Optionpk<string> label_opt("cn", "cname", "name of the class label in the output vector file", "label");
  Optionpk<bool> polygon_opt("polygon", "polygon", "create OGRPolygon as geometry instead of OGRPoint. Only if sample option is also of polygon type.", false);
  Optionpk<int> band_opt("b", "band", "band index to crop. Use -1 to use all bands)", -1);
  Optionpk<string> rule_opt("r", "rule", "rule how to report image information per feature. point (value at each point or at centroid if polygon), centroid, mean (of polygon), proportion, minimum (of polygon), maximum (of polygon), maxvote, sum.", "point");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=image_opt.retrieveOption(argc,argv);
    sample_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    ogrformat_opt.retrieveOption(argc,argv);
    test_opt.retrieveOption(argc,argv);
    // bufferOutput_opt.retrieveOption(argc,argv);
    geo_opt.retrieveOption(argc,argv);
    down_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    min_opt.retrieveOption(argc,argv);
    boundary_opt.retrieveOption(argc,argv);
    // rbox_opt.retrieveOption(argc,argv);
    // cbox_opt.retrieveOption(argc,argv);
    disc_opt.retrieveOption(argc,argv);
    ftype_opt.retrieveOption(argc,argv);
    ltype_opt.retrieveOption(argc,argv);
    fieldname_opt.retrieveOption(argc,argv);
    label_opt.retrieveOption(argc,argv);
    polygon_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    rule_opt.retrieveOption(argc,argv);
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

  std::map<std::string, rule::RULE_TYPE> ruleMap;
  //initialize ruleMap
  ruleMap["point"]=rule::point;
  ruleMap["centroid"]=rule::centroid;
  ruleMap["mean"]=rule::mean;
  ruleMap["proportion"]=rule::proportion;
  ruleMap["minimum"]=rule::minimum;
  ruleMap["maximum"]=rule::maximum;
  ruleMap["custom"]=rule::custom;
  ruleMap["maxvote"]=rule::maxvote;
  ruleMap["sum"]=rule::sum;

  if(verbose_opt[0])
    std::cout << class_opt << std::endl;
  statfactory::StatFactory stat;
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
  short theDim=boundary_opt[0];
  if(verbose_opt[0]>1)
    std::cout << "boundary: " << boundary_opt[0] << std::endl;
  ImgReaderGdal imgReader;
  if(image_opt.empty()){
    std::cerr << "No image file provided (use option -i). Use --help for help information";
      exit(0);
  }
  if(output_opt.empty()){
    std::cerr << "No output file provided (use option -o). Use --help for help information";
      exit(0);
  }
  if(sample_opt.empty()){
    std::cerr << "No sample file provided (use option -s). Use --help for help information";
    exit(0);
  }
  try{
    imgReader.open(image_opt[0]);
  }
  catch(std::string errorstring){
    std::cout << errorstring << std::endl;
    exit(0);
  }
  int nband=(band_opt[0]<0)?imgReader.nrOfBand():band_opt.size();

  if(fieldname_opt.size()<nband){
    std::string bandString=fieldname_opt[0];
    fieldname_opt.clear();
    fieldname_opt.resize(nband);
    for(int iband=0;iband<nband;++iband){
      int theBand=(band_opt[0]<0)?iband:band_opt[iband];
      ostringstream fs;
      fs << bandString << theBand;
      fieldname_opt[iband]=fs.str();
    }
  }

  if(verbose_opt[0])
    std::cout << fieldname_opt << std::endl;
  vector<ImgReaderGdal> maskReader;
  if(mask_opt.size()){
    maskReader.resize(mask_opt.size());
    for(int imask=0;imask<mask_opt.size();++imask){
      if(verbose_opt[0]>1)
        std::cout << "opening mask image file " << mask_opt[imask] << std::endl;
      maskReader[imask].open(mask_opt[0]);
      if(imgReader.isGeoRef())
        assert(maskReader[imask].isGeoRef());
    }
  }

  Vector2d<int> maskBuffer;
  if(mask_opt.size()){
    maskBuffer.resize(mask_opt.size());
    for(int imask=0;imask<maskReader.size();++imask)
      maskBuffer[imask].resize(maskReader[imask].nrOfCol());
  }
  vector<double> oldmaskrow(mask_opt.size());
  for(int imask=0;imask<mask_opt.size();++imask)
    oldmaskrow[imask]=-1;
  
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
  srandom(time(NULL));

  bool sampleIsRaster=false;
  ImgReaderOgr sampleReaderOgr;
  try{
    sampleReaderOgr.open(sample_opt[0]);
    //test
    // int nlayer=sampleReaderOgr.getDataSource()->GetLayerCount();
    // cout << "nlayer: " << nlayer << endl;
  }
  catch(string errorString){
    sampleIsRaster=true;
  }

  if(sampleIsRaster){
    if(class_opt.empty()){
      // std::cout << "Warning: no classes selected, if a classes must be extracted, set to -1 for all classes using option -c -1" << std::endl;
      ImgReaderGdal classReader;
      ImgWriterOgr ogrWriter;
      // if(verbose_opt[0]>1){
      //   std::cout << "reading position from " << sample_opt[0] << std::endl;
      //   std::cout << "class thresholds: " << std::endl;
      //   for(int iclass=0;iclass<class_opt.size();++iclass){
      //     if(threshold_opt.size()>1)
      //       std::cout << class_opt[iclass] << ": " << threshold_opt[iclass] << std::endl;
      //     else
      //       std::cout << class_opt[iclass] << ": " << threshold_opt[0] << std::endl;
      //   }
      // }
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
            if(static_cast<int>(jimg)!=static_cast<int>(oldimgrow)){
              assert(imgBuffer.size()==nband);
              for(int iband=0;iband<nband;++iband){
                int theBand=(band_opt[0]<0)?iband:band_opt[iband];
                imgReader.readData(imgBuffer[iband],GDT_Float64,static_cast<int>(jimg),theBand);
                assert(imgBuffer[iband].size()==imgReader.nrOfCol());
              }
              oldimgrow=jimg;
            }
            // bool valid=true;
            for(int imask=0;imask<mask_opt.size();++imask){
              double colMask,rowMask;//image coordinates in mask image
              if(mask_opt.size()>1){//multiple masks
                if(geo_opt[0])
                  maskReader[imask].geo2image(x,y,colMask,rowMask);
                else{
                  colMask=icol;
                  rowMask=irow;
                }
                //nearest neighbour
                rowMask=static_cast<int>(rowMask);
                colMask=static_cast<int>(colMask);
                if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[imask].nrOfCol())
                  continue;
                if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[imask])){
                  if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[imask].nrOfRow())
                    continue;
                  else{
                    maskReader[imask].readData(maskBuffer[imask],GDT_Int32,static_cast<int>(rowMask));
                    oldmaskrow[imask]=rowMask;
                  }
                }
                int ivalue=0;
                if(mask_opt.size()==msknodata_opt.size())//one invalid value for each mask
                  ivalue=static_cast<int>(msknodata_opt[imask]);
                else//use same invalid value for each mask
                  ivalue=static_cast<int>(msknodata_opt[0]);
                if(maskBuffer[imask][colMask]==ivalue){
                  valid=false;
                  break;
                }
              }
              else if(maskReader.size()){
                if(geo_opt[0])
                  maskReader[0].geo2image(x,y,colMask,rowMask);
                else{
                  colMask=icol;
                  rowMask=irow;
                }
                //nearest neighbour
                rowMask=static_cast<int>(rowMask);
                colMask=static_cast<int>(colMask);
                if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[0].nrOfCol())
                  continue;
                if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[0])){
                  if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[0].nrOfRow())
                    continue;
                  else{
                    maskReader[0].readData(maskBuffer[0],GDT_Int32,static_cast<int>(rowMask));
                    oldmaskrow[0]=rowMask;
                  }
                }
                for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
                  if(maskBuffer[0][colMask]==static_cast<int>(msknodata_opt[ivalue])){
                    valid=false;
                    break;
                  }
                }
              }
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
                double p=static_cast<double>(random())/(RAND_MAX);
                p*=100.0;
                if(p>theThreshold)
		  continue;//do not select for now, go to next column
              }
              else if(nvalid.size()>processClass){//absolute value
                if(nvalid[processClass]>-theThreshold)
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
          int theBand=(band_opt[0]<0)?iband:band_opt[iband];
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
            int theBand=(band_opt[0]<0)?iband:band_opt[iband];
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
        std::cout << "reading position from " << sample_opt[0] << std::endl;
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
            if(static_cast<int>(jimg)!=static_cast<int>(oldimgrow)){
              assert(imgBuffer.size()==nband);
              for(int iband=0;iband<nband;++iband){
                int theBand=(band_opt[0]<0)?iband:band_opt[iband];
                imgReader.readData(imgBuffer[iband],GDT_Float64,static_cast<int>(jimg),theBand);
                assert(imgBuffer[iband].size()==imgReader.nrOfCol());
              }
              oldimgrow=jimg;
            }
            bool valid=true;
            for(int imask=0;imask<mask_opt.size();++imask){
              double colMask,rowMask;//image coordinates in mask image
              if(mask_opt.size()>1){//multiple masks
                if(geo_opt[0])
                  maskReader[imask].geo2image(x,y,colMask,rowMask);
                else{
                  colMask=icol;
                  rowMask=irow;
                }
                //nearest neighbour
                rowMask=static_cast<int>(rowMask);
                colMask=static_cast<int>(colMask);
                if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[imask].nrOfCol())
                  continue;
                if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[imask])){
                  if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[imask].nrOfRow())
                    continue;
                  else{
                    maskReader[imask].readData(maskBuffer[imask],GDT_Int32,static_cast<int>(rowMask));
                    oldmaskrow[imask]=rowMask;
                  }
                }
                int ivalue=0;
                if(mask_opt.size()==msknodata_opt.size())//one invalid value for each mask
                  ivalue=static_cast<int>(msknodata_opt[imask]);
                else//use same invalid value for each mask
                  ivalue=static_cast<int>(msknodata_opt[0]);
                if(maskBuffer[imask][colMask]==ivalue){
                  valid=false;
                  break;
                }
              }
              else if(maskReader.size()){
                if(geo_opt[0])
                  maskReader[0].geo2image(x,y,colMask,rowMask);
                else{
                  colMask=icol;
                  rowMask=irow;
                }
                //nearest neighbour
                rowMask=static_cast<int>(rowMask);
                colMask=static_cast<int>(colMask);
                if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[0].nrOfCol())
                  continue;
                if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[0])){
                  if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[0].nrOfRow())
                    continue;
                  else{
                    maskReader[0].readData(maskBuffer[0],GDT_Int32,static_cast<int>(rowMask));
                    oldmaskrow[0]=rowMask;
                  }
                }
                for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
                  if(maskBuffer[0][colMask]==static_cast<int>(msknodata_opt[ivalue])){
                    valid=false;
                    break;
                  }
                }
              }
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
                double p=static_cast<double>(random())/(RAND_MAX);
                p*=100.0;
                if(p>theThreshold)
                  continue;//do not select for now, go to next column
              }
              else if(nvalid.size()>processClass){//absolute value
                if(nvalid[processClass]>-theThreshold)
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
          int theBand=(band_opt[0]<0)?iband:band_opt[iband];
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
            int theBand=(band_opt[0]<0)?iband:band_opt[iband];
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
  else{//vector file
    if(verbose_opt[0]>1)
      std::cout << "creating image sample writer " << output_opt[0] << std::endl;
    ImgWriterOgr ogrWriter;
    ImgWriterOgr ogrTestWriter;
    ogrWriter.open(output_opt[0],ogrformat_opt[0]);
    if(test_opt.size()){
      if(verbose_opt[0]>1)
	std::cout << "creating image test writer " << test_opt[0] << std::endl;
      ogrTestWriter.open(test_opt[0],ogrformat_opt[0]);
    }
      
    //support multiple layers
    int nlayer=sampleReaderOgr.getDataSource()->GetLayerCount();
    if(verbose_opt[0])
      std::cout << "number of layers: " << nlayer << endl;
      
    for(int ilayer=0;ilayer<nlayer;++ilayer){
      OGRLayer *readLayer=sampleReaderOgr.getLayer(ilayer);
      cout << "processing layer " << readLayer->GetName() << endl;
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
	std::cout << "copy fields" << std::flush << std::endl;
      ogrWriter.copyFields(sampleReaderOgr,ilayer);

      if(test_opt.size()){
	if(verbose_opt[0])
	  std::cout << "copy fields test writer" << std::flush << std::endl;
	ogrTestWriter.copyFields(sampleReaderOgr,ilayer);
      }
      vector<std::string> fieldnames;
      if(verbose_opt[0])
	std::cout << "get fields" << std::flush << std::endl;
      sampleReaderOgr.getFields(fieldnames);
      assert(fieldnames.size()==ogrWriter.getFieldCount(ilayer));
      map<std::string,double> pointAttributes;

      switch(ruleMap[rule_opt[0]]){
	// switch(rule_opt[0]){
      case(rule::proportion):{//proportion for each class
	assert(class_opt.size());
	for(int iclass=0;iclass<class_opt.size();++iclass){
	  ostringstream cs;
	  cs << class_opt[iclass];
	  ogrWriter.createField(cs.str(),fieldType,ilayer);
	}
	break;
      }
      case(rule::custom):
      case(rule::minimum):
      case(rule::maximum):
      case(rule::maxvote):
	assert(class_opt.size());
      ogrWriter.createField(label_opt[0],fieldType,ilayer);
      if(test_opt.size())
	ogrTestWriter.createField(label_opt[0],fieldType,ilayer);
      break;
      case(rule::point):
      case(rule::mean):
      case(rule::sum):
      case(rule::centroid):
      default:{
	for(int windowJ=-theDim/2;windowJ<(theDim+1)/2;++windowJ){
	  for(int windowI=-theDim/2;windowI<(theDim+1)/2;++windowI){
	    if(disc_opt[0]&&(windowI*windowI+windowJ*windowJ>(theDim/2)*(theDim/2)))
	      continue;
	    for(int iband=0;iband<nband;++iband){
	      int theBand=(band_opt[0]<0)?iband:band_opt[iband];
	      ostringstream fs;
	      if(theDim>1)
		fs << fieldname_opt[iband] << "_" << windowJ << "_" << windowI;
	      else
		fs << fieldname_opt[iband];
	      if(verbose_opt[0]>1)
		std::cout << "creating field " << fs.str() << std::endl;

	      ogrWriter.createField(fs.str(),fieldType,ilayer);
	      if(test_opt.size())
		ogrTestWriter.createField(fs.str(),fieldType,ilayer);
	    }
	  }
	}
	break;
      }
      }
      OGRFeature *readFeature;
      unsigned long int ifeature=0;
      unsigned long int nfeature=sampleReaderOgr.getFeatureCount();
      // ImgWriterOgr boxWriter;
      // if(rbox_opt[0]>0||cbox_opt[0]>0){
      // 	assert(bufferOutput_opt.size());
      // 	assert(test_opt.empty());//not implemented
      //   if(verbose_opt[0]>1)
      //     std::cout << "opening box writer " << bufferOutput_opt[0] << std::endl;
      //   boxWriter.open(bufferOutput_opt[0]);
      //   std::string layername="buffer";
      //   boxWriter.createLayer(layername, imgReader.getProjection(), wkbPolygon);
      //   std::string fieldname="fid";//number of the point
      //   if(verbose_opt[0]>1)
      //     std::cout << "creating field " << fieldname << std::endl;
      //   //       ogrWriter.createField(fieldname,OFTInteger,ilayer);
      //   boxWriter.createField(fieldname,OFTInteger,ilayer);
      // }
      progress=0;
      pfnProgress(progress,pszMessage,pProgressArg);
      while( (readFeature = readLayer->GetNextFeature()) != NULL ){
	bool validFeature=false;
	bool writeTest=false;//write this feature to test_opt[0] instead of output_opt
	if(verbose_opt[0]>0)
	  std::cout << "reading feature " << readFeature->GetFID() << std::endl;
	if(threshold_opt[0]>0){//percentual value
	  double p=static_cast<double>(random())/(RAND_MAX);
	  p*=100.0;
	  if(p>threshold_opt[0]){
	    if(test_opt.size())
	      writeTest=true;
	    else
	      continue;//do not select for now, go to next feature
	  }
	}
	else{//absolute value
	  if(ntotalvalid>-threshold_opt[0]){
	    if(test_opt.size())
	      writeTest=true;
	    else
	      continue;//do not select any more pixels, go to next column feature
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
	    assert(class_opt.size()<=1);//class_opt not implemented for point yet
	    OGRPoint *poPoint = (OGRPoint *) poGeometry;
	    x=poPoint->getX();
	    y=poPoint->getY();
	    bool valid=true;
	    for(int imask=0;imask<mask_opt.size();++imask){
	      double colMask,rowMask;//image coordinates in mask image
	      if(mask_opt.size()>1){//multiple masks
		maskReader[imask].geo2image(x,y,colMask,rowMask);
		//nearest neighbour
		rowMask=static_cast<int>(rowMask);
		colMask=static_cast<int>(colMask);
		if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[imask].nrOfCol())
		  continue;
		if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[imask])){
		  if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[imask].nrOfRow())
		    continue;
		  else{
		    maskReader[imask].readData(maskBuffer[imask],GDT_Int32,static_cast<int>(rowMask));
		    oldmaskrow[imask]=rowMask;
		    assert(maskBuffer.size()==maskReader[imask].nrOfBand());
		  }
		}
		//               char ivalue=0;
		int ivalue=0;
		if(mask_opt.size()==msknodata_opt.size())//one invalid value for each mask
		  ivalue=static_cast<int>(msknodata_opt[imask]);
		else//use same invalid value for each mask
		  ivalue=static_cast<int>(msknodata_opt[0]);
		if(maskBuffer[imask][colMask]==ivalue){
		  valid=false;
		  break;
		}
	      }
	      else if(maskReader.size()){
		maskReader[0].geo2image(x,y,colMask,rowMask);
		//nearest neighbour
		rowMask=static_cast<int>(rowMask);
		colMask=static_cast<int>(colMask);
		if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[0].nrOfCol()){
		  continue;
		  // cerr << colMask << " out of mask col range!" << std::endl;
		  // cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
		  // assert(static_cast<int>(colMask)>=0&&static_cast<int>(colMask)<maskReader[0].nrOfCol());
		}
              
		if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[0])){
		  if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[0].nrOfRow()){
		    continue;
		    // cerr << rowMask << " out of mask row range!" << std::endl;
		    // cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
		    // assert(static_cast<int>(rowMask)>=0&&static_cast<int>(rowMask)<imgReader.nrOfRow());
		  }
		  else{
		    maskReader[0].readData(maskBuffer[0],GDT_Int32,static_cast<int>(rowMask));
		    oldmaskrow[0]=rowMask;
		  }
		}
		for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
		  if(maskBuffer[0][colMask]==static_cast<int>(msknodata_opt[ivalue])){
		    valid=false;
		    break;
		  }
		}
	      }
	    }
	    if(!valid)
	      continue;
	    else
	      validFeature=true;

	    double value;
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
	    //check if j_centre is out of bounds
	    if(static_cast<int>(j_centre)<0||static_cast<int>(j_centre)>=imgReader.nrOfRow())
	      continue;
	    //check if i_centre is out of bounds
	    if(static_cast<int>(i_centre)<0||static_cast<int>(i_centre)>=imgReader.nrOfCol())
	      continue;

	    // if(rbox_opt[0]){
	    //   assert(test_opt.empty());//not implemented
	    //   vector< vector<OGRPoint*> > points;
	    //   points.resize(rbox_opt.size());
	    //   if(verbose_opt[0]>1)
	    //     std::cout << "creating rectangular box for sample " << isample << ": ";
	    //   for(int ibox=0;ibox<rbox_opt.size();++ibox){
	    //     int npoint=4;
	    //     if(verbose_opt[0]>1)
	    //       std::cout << ibox << " ";
	    //     points[ibox].resize(npoint+1);
	    //     vector<OGRPoint> pbPoint(npoint+1);
	    //     pbPoint[0].setX(x-0.5*rbox_opt[ibox]);
	    //     pbPoint[0].setY(y+0.5*rbox_opt[ibox]);
	    //     points[ibox][0]=&(pbPoint[0]);//start point UL
	    //     points[ibox][4]=&(pbPoint[0]);//end point
	    //     pbPoint[1].setX(x+0.5*rbox_opt[ibox]);
	    //     pbPoint[1].setY(y+0.5*rbox_opt[ibox]);
	    //     points[ibox][1]=&(pbPoint[1]);//UR
	    //     pbPoint[2].setX(x+0.5*rbox_opt[ibox]);
	    //     pbPoint[2].setY(y-0.5*rbox_opt[ibox]);
	    //     points[ibox][2]=&(pbPoint[2]);//LR
	    //     pbPoint[3].setX(x-0.5*rbox_opt[ibox]);
	    //     pbPoint[3].setY(y-0.5*rbox_opt[ibox]);
	    //     points[ibox][3]=&(pbPoint[3]);//LL
	    //     std::string fieldname="fid";//number of the point
	    //     boxWriter.addRing(points[ibox],fieldname,isample);
	    //     // boxWriter.addLineString(points[ibox],fieldname,isample);
	    //   }
	    //   if(verbose_opt[0]>1)
	    //     std::cout << std::endl;
	    // }
	    // if(cbox_opt[0]>0){
	    //   vector< vector<OGRPoint*> > points;
	    //   points.resize(cbox_opt.size());
	    //   if(verbose_opt[0]>1)
	    //     std::cout << "creating circular box ";
	    //   for(int ibox=0;ibox<cbox_opt.size();++ibox){
	    //     int npoint=50;
	    //     if(verbose_opt[0]>1)
	    //       std::cout << ibox << " ";
	    //     points[ibox].resize(npoint+1);
	    //     vector<OGRPoint> pbPoint(npoint+1);
	    //     double radius=cbox_opt[ibox]/2.0;
	    //     double alpha=0;
	    //     for(int ipoint=0;ipoint<npoint;++ipoint){
	    //       alpha=ipoint*2.0*PI/static_cast<double>(npoint);
	    //       pbPoint[ipoint].setX(x+radius*cos(alpha));
	    //       pbPoint[ipoint].setY(y+radius*sin(alpha));
	    //       points[ibox][ipoint]=&(pbPoint[ipoint]);
	    //     }
	    //     alpha=0;
	    //     pbPoint[npoint].setX(x+radius*cos(alpha));
	    //     pbPoint[npoint].setY(y+radius*sin(alpha));
	    //     points[ibox][npoint]=&(pbPoint[npoint]);
	    //     std::string fieldname="fid";//number of the point
	    //     boxWriter.addRing(points[ibox],fieldname,isample);
	    //     // boxWriter.addLineString(points[ibox],fieldname,isample);
	    //   }
	    //   if(verbose_opt[0]>1)
	    //     std::cout << std::endl;
	    // }
      
	    OGRFeature *writeFeature;
	    if(verbose_opt[0]>1)
	      std::cout << "create feature " << sample_opt[0] << std::endl;
	    if(writeTest)
	      writeFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	    else
	      writeFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    if(verbose_opt[0]>1)
	      std::cout << "copying fields from points " << sample_opt[0] << std::endl;
	    if(writeFeature->SetFrom(readFeature)!= OGRERR_NONE)
	      cerr << "writing feature failed" << std::endl;

	    if(verbose_opt[0]>1)
	      std::cout << "write feature has " << writeFeature->GetFieldCount() << " fields" << std::endl;

	    vector<double> windowBuffer;
	    for(int windowJ=-theDim/2;windowJ<(theDim+1)/2;++windowJ){
	      for(int windowI=-theDim/2;windowI<(theDim+1)/2;++windowI){
		if(disc_opt[0]&&(windowI*windowI+windowJ*windowJ>(theDim/2)*(theDim/2)))
		  continue;
		int j=j_centre+windowJ;
		//check if j is out of bounds
		if(static_cast<int>(j)<0||static_cast<int>(j)>=imgReader.nrOfRow())
		  continue;
		int i=i_centre+windowI;
		//check if i is out of bounds
		if(static_cast<int>(i)<0||static_cast<int>(i)>=imgReader.nrOfCol())
		  continue;
		if(verbose_opt[0]>1)
		  std::cout << "reading image value at " << i << "," << j;
		for(int iband=0;iband<nband;++iband){
		  int theBand=(band_opt[0]<0)?iband:band_opt[iband];
		  imgReader.readData(value,GDT_Float64,i,j,theBand);
		  if(verbose_opt[0]>1)
		    std::cout << ": " << value << std::endl;
		  ostringstream fs;
		  if(theDim>1)
		    fs << fieldname_opt[iband] << "_" << windowJ << "_" << windowI;
		  else
		    fs << fieldname_opt[iband];
		  if(verbose_opt[0]>1)
		    std::cout << "set field " << fs.str() << " to " << value << std::endl;
		  switch( fieldType ){
		  case OFTInteger:
		    writeFeature->SetField(fs.str().c_str(),static_cast<int>(value));
		    break;
		  case OFTString:
		    {
		      ostringstream os;
		      os << value;
		      writeFeature->SetField(fs.str().c_str(),os.str().c_str());
		      break;
		    }
		  case OFTReal:
		    writeFeature->SetField(fs.str().c_str(),value);
		    break;
		  case OFTRealList:{
		    int fieldIndex=writeFeature->GetFieldIndex(fs.str().c_str());
		    int nCount;
		    const double *theList;
		    theList=writeFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		    vector<double> vectorList(nCount+1);
		    for(int index=0;index<nCount;++index)
		      vectorList[nCount]=value;
		    writeFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    break;
		  }
		  default://not supported
		    assert(0);
		    break;
		  }
		}
	      }
	    }
	    if(verbose_opt[0]>1)
	      std::cout << "creating point feature" << std::endl;
	    if(writeTest){
	      if(writeTestLayer->CreateFeature( writeFeature ) != OGRERR_NONE ){
		std::string errorString="Failed to create feature in shapefile";
		throw(errorString);
	      }
	    }
	    else{
	      if(writeLayer->CreateFeature( writeFeature ) != OGRERR_NONE ){
		std::string errorString="Failed to create feature in shapefile";
		throw(errorString);
	      }
	    }
	    OGRFeature::DestroyFeature( writeFeature );
	    // ++isample;
	    ++ntotalvalid;
	    if(verbose_opt[0])
	      std::cout << "ntotalvalid: " << ntotalvalid << std::endl;
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

	    if(uli<0||lri>=imgReader.nrOfCol()||ulj<0||ulj>=imgReader.nrOfRow())
	      continue;

	    int nPointPolygon=0;
	    if(polygon_opt[0]){
	      if(writeTest)
		writePolygonFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    else if(ruleMap[rule_opt[0]]==rule::mean||ruleMap[rule_opt[0]]==rule::centroid||ruleMap[rule_opt[0]]==rule::sum){
	      if(writeTest)
		writeCentroidFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writeCentroidFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    //previously here
	    vector<double> polyValues;
	    switch(ruleMap[rule_opt[0]]){
	    case(rule::point):
	    case(rule::mean):
	    case(rule::sum):
	    case(rule::centroid):
	    default:
	      polyValues.resize(nband);
	    break;
	    case(rule::proportion):
	    case(rule::custom):
	    case(rule::minimum):
	    case(rule::maximum):
	    case(rule::maxvote):
	      assert(class_opt.size());
	    polyValues.resize(class_opt.size());
	    break;
	    }
	    for(int index=0;index<polyValues.size();++index)
	      polyValues[index]=0;
	    OGRPoint thePoint;
	    for(int j=ulj;j<=lrj;++j){
	      for(int i=uli;i<=lri;++i){
		//check if point is on surface
		double x=0;
		double y=0;
		imgReader.image2geo(i,j,x,y);
		thePoint.setX(x);
		thePoint.setY(y);
		if(readPolygon.Contains(&thePoint)){
		  bool valid=true;
		  for(int imask=0;imask<mask_opt.size();++imask){
		    double colMask,rowMask;//image coordinates in mask image
		    if(mask_opt.size()>1){//multiple masks
		      maskReader[imask].geo2image(x,y,colMask,rowMask);
		      //nearest neighbour
		      rowMask=static_cast<int>(rowMask);
		      colMask=static_cast<int>(colMask);
		      if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[imask].nrOfCol())
			continue;
		      // {
		      //   cerr << colMask << " out of mask col range!" << std::endl;
		      //   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
		      //   assert(static_cast<int>(colMask)>=0&&static_cast<int>(colMask)<maskReader[imask].nrOfCol());
		      // }
              
		      if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[imask])){
			if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[imask].nrOfRow())
			  continue;
			// {
			//   cerr << rowMask << " out of mask row range!" << std::endl;
			//   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
			//   assert(static_cast<int>(rowMask)>=0&&static_cast<int>(rowMask)<imgReader.nrOfRow());
			// }
			else{
			  maskReader[imask].readData(maskBuffer[imask],GDT_Int32,static_cast<int>(rowMask));
			  oldmaskrow[imask]=rowMask;
			  assert(maskBuffer.size()==maskReader[imask].nrOfBand());
			}
		      }
		      //               char ivalue=0;
		      int ivalue=0;
		      if(mask_opt.size()==msknodata_opt.size())//one invalid value for each mask
			ivalue=static_cast<int>(msknodata_opt[imask]);
		      else//use same invalid value for each mask
			ivalue=static_cast<int>(msknodata_opt[0]);
		      if(maskBuffer[imask][colMask]==ivalue){
			valid=false;
			break;
		      }
		    }
		    else if(maskReader.size()){
		      maskReader[0].geo2image(x,y,colMask,rowMask);
		      //nearest neighbour
		      rowMask=static_cast<int>(rowMask);
		      colMask=static_cast<int>(colMask);
		      if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[0].nrOfCol())
			continue;
		      // {
		      //   cerr << colMask << " out of mask col range!" << std::endl;
		      //   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
		      //   assert(static_cast<int>(colMask)>=0&&static_cast<int>(colMask)<maskReader[0].nrOfCol());
		      // }
              
		      if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[0])){
			if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[0].nrOfRow())
			  continue;
			// {
			//   cerr << rowMask << " out of mask row range!" << std::endl;
			//   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
			//   assert(static_cast<int>(rowMask)>=0&&static_cast<int>(rowMask)<imgReader.nrOfRow());
			// }
			else{
			  maskReader[0].readData(maskBuffer[0],GDT_Int32,static_cast<int>(rowMask));
			  oldmaskrow[0]=rowMask;
			}
		      }
		      for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
			if(maskBuffer[0][colMask]==static_cast<int>(msknodata_opt[ivalue])){
			  valid=false;
			  break;
			}
		      }
		    }
		  }
		  if(!valid)
		    continue;
		  else
		    validFeature=true;
		  //check if within raster image
		  if(i<0||i>=imgReader.nrOfCol())
		    continue;
		  if(j<0||j>=imgReader.nrOfRow())
		    continue;
		  writeRing.addPoint(&thePoint);
		  if(verbose_opt[0]>1)
		    std::cout << "point is on surface:" << thePoint.getX() << "," << thePoint.getY() << std::endl;
		  ++nPointPolygon;
		  OGRFeature *writePointFeature;
		  if(!polygon_opt[0]){
		    //create feature
		    if(ruleMap[rule_opt[0]]!=rule::mean&&ruleMap[rule_opt[0]]!=rule::centroid&&ruleMap[rule_opt[0]]!=rule::sum){//do not create in case of mean, sum or centroid (only create point at centroid)
		      if(writeTest)
			writePointFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
		      else
			writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
		      if(verbose_opt[0]>1)
			std::cout << "copying fields from polygons " << sample_opt[0] << std::endl;
		      if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
			cerr << "writing feature failed" << std::endl;
		      writePointFeature->SetGeometry(&thePoint);
		      OGRGeometry *updateGeometry;
		      updateGeometry = writePointFeature->GetGeometryRef();
		      OGRPoint *poPoint = (OGRPoint *) updateGeometry;
		      if(verbose_opt[0]>1)
			std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields" << std::endl;
		    }
		  }
		  if(verbose_opt[0]>1)
		    std::cout << "reading image value within polygon at position " << i << "," << j;
		  for(int iband=0;iband<nband;++iband){
		    int theBand=(band_opt[0]<0)?iband:band_opt[iband];
		    double value=0;
		    imgReader.readData(value,GDT_Float64,i,j,theBand);
		    if(verbose_opt[0]>1)
		      std::cout << ": " << value << std::endl;
		    if(polygon_opt[0]||ruleMap[rule_opt[0]]==rule::mean||ruleMap[rule_opt[0]]==rule::centroid||ruleMap[rule_opt[0]]==rule::sum){
		      int iclass=0;
		      switch(ruleMap[rule_opt[0]]){
		      case(rule::point)://in centroid if polygon_opt==true or all values as points if polygon_opt!=true
		      case(rule::centroid):
		      default:
			polyValues[iband]=value;
		      break;
		      case(rule::sum):
		      case(rule::mean)://mean as polygon if polygon_opt==true or as point in centroid if polygon_opt!=true
			polyValues[iband]+=value;
			break;
		      case(rule::proportion):
		      case(rule::custom):
		      case(rule::minimum):
		      case(rule::maximum):
		      case(rule::maxvote):
			for(iclass=0;iclass<class_opt.size();++iclass){
			  if(value==class_opt[iclass]){
			    assert(polyValues.size()>iclass);
			    polyValues[iclass]+=1;//value
			    break;
			  }
			}
		      break;
		      }
		    }
		    else{
		      // ostringstream fs;
		      // if(imgReader.nrOfBand()==1)
		      //   fs << fieldname_opt[0];
		      // else
		      //   fs << fieldname_opt[0] << theBand;
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[iband] << " to " << value << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
			writePointFeature->SetField(fieldname_opt[iband].c_str(),static_cast<int>(value));
			break;
		      case OFTString:
			{
			  ostringstream os;
			  os << value;
			  writePointFeature->SetField(fieldname_opt[iband].c_str(),os.str().c_str());
			  break;
			}
		      case OFTReal:
			writePointFeature->SetField(fieldname_opt[iband].c_str(),value);
			break;
		      case OFTRealList:{
			int fieldIndex=writePointFeature->GetFieldIndex(fieldname_opt[iband].c_str());
			int nCount;
			const double *theList;
			theList=writePointFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
			vector<double> vectorList(nCount+1);
			for(int index=0;index<nCount;++index)
			  vectorList[nCount]=value;
			writePointFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
			break;
		      }
		      default://not supported
			assert(0);
			break;
		      }
		    }
		  }
		  if(!polygon_opt[0]){
		    if(ruleMap[rule_opt[0]]!=rule::mean&&ruleMap[rule_opt[0]]!=rule::centroid&&ruleMap[rule_opt[0]]!=rule::sum){//do not create in case of mean value (only at centroid)
		      //write feature
		      if(verbose_opt[0]>1)
			std::cout << "creating point feature" << std::endl;
		      if(writeTest){
			if(writeTestLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			  std::string errorString="Failed to create feature in test shapefile";
			  throw(errorString);
			}
		      }
		      else{
			if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			  std::string errorString="Failed to create feature in shapefile";
			  throw(errorString);
			}
		      }
		      //destroy feature
		      OGRFeature::DestroyFeature( writePointFeature );
		      ++ntotalvalid;
		      if(verbose_opt[0])
			std::cout << "ntotalvalid(2): " << ntotalvalid << std::endl;
		    }
		  }
		  // ++isample;
		}
	      }
	    }
	    if(polygon_opt[0]||ruleMap[rule_opt[0]]==rule::mean||ruleMap[rule_opt[0]]==rule::centroid||ruleMap[rule_opt[0]]==rule::sum){
	      //do not create if no points found within polygon
	      if(!nPointPolygon)
		continue;
	      //add ring to polygon
	      if(polygon_opt[0]){
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
		//write geometry of writePolygon
		writePolygonFeature->SetGeometry(&writePolygon);
		if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "copying new fields write polygon " << sample_opt[0] << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
		//write polygon feature
	      }
	      else{//write mean value of polygon to centroid point (ruleMap[rule_opt[0]]==rule::mean)
		//create feature
		if(verbose_opt[0]>1)
		  std::cout << "copying fields from polygons " << sample_opt[0] << std::endl;
		if(writeCentroidFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		writeCentroidFeature->SetGeometry(&writeCentroidPoint);
		OGRGeometry *updateGeometry;
		updateGeometry = writeCentroidFeature->GetGeometryRef();
		assert(wkbFlatten(updateGeometry->getGeometryType()) == wkbPoint );
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writeCentroidFeature->GetFieldCount() << " fields" << std::endl;
	      }
	      switch(ruleMap[rule_opt[0]]){
	      case(rule::point)://value at each point (or at centroid of polygon if line is set)
	      default:{
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		for(int index=0;index<polyValues.size();++index){
		  double theValue=polyValues[index];
		  // ostringstream fs;
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  int theBand=(band_opt[0]<0)?index:band_opt[index];
		  // if(nband==1)
		  //   fs << fieldname_opt[0];
		  // else
		  //   fs << fieldname_opt[0] << theBand;

		  if(verbose_opt[0]>1)
		    std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		  switch( fieldType ){
		  case OFTInteger:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    break;
		  case OFTString:
		    {
		      ostringstream os;
		      os << theValue;
		      if(polygon_opt[0])
			writePolygonFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      else
			writeCentroidFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      break;
		    }
		  case OFTReal:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    break;
		  case OFTRealList:{
		    int fieldIndex;
		    int nCount;
		    const double *theList;
		    vector<double> vectorList;
		    if(polygon_opt[0]){
		      fieldIndex=writePolygonFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writePolygonFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writePolygonFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    else{
		      fieldIndex=writeCentroidFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writeCentroidFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writeCentroidFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    break;
		  }
		  default://not supported
		    std::cout << "field type not supported yet..." << std::endl;
		    break;
		  }
		}
		break;
	      }
	      case(rule::mean):
	      case(rule::sum):
	      case(rule::centroid):{//mean value (written to centroid of polygon if line is not set)
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//test
		if(ruleMap[rule_opt[0]]==rule::centroid)
		  assert(nPointPolygon<=1);
		for(int index=0;index<polyValues.size();++index){
		  double theValue=polyValues[index];
		  // ostringstream fs;
		  if(ruleMap[rule_opt[0]]==rule::mean)
		    theValue/=nPointPolygon;
		  int theBand=(band_opt[0]<0)?index:band_opt[index];
		  // if(nband==1)
		  //   fs << fieldname_opt[0];
		  // else
		  //   fs << fieldname_opt[0] << theBand;
		  if(verbose_opt[0]>1)
		    std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		  switch( fieldType ){
		  case OFTInteger:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    break;
		  case OFTString:
		    {
		      ostringstream os;
		      os << theValue;
		      if(polygon_opt[0])
			writePolygonFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      else
			writeCentroidFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      break;
		    }
		  case OFTReal:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    break;
		  case OFTRealList:{
		    int fieldIndex;
		    int nCount;
		    const double *theList;
		    vector<double> vectorList;
		    if(polygon_opt[0]){
		      fieldIndex=writePolygonFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writePolygonFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writePolygonFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    else{
		      fieldIndex=writeCentroidFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writeCentroidFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writeCentroidFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    break;
		  }
		  default://not supported
		    std::cout << "field type not supported yet..." << std::endl;
		    break;
		  }
		}
		break;
	      }
	      case(rule::proportion):{//proportion classes
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		stat.normalize_pct(polyValues);
		// stat.sum(polyValues);
		for(int index=0;index<polyValues.size();++index){
		  double theValue=polyValues[index];
		  ostringstream fs;
		  fs << class_opt[index];
		  if(polygon_opt[0])
		    writePolygonFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		  else
		    writeCentroidFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		}
		break;
	      }
	      case(rule::custom):{//custom
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		stat.normalize_pct(polyValues);
		assert(polyValues.size()==2);//11:broadleaved, 12:coniferous
		if(polyValues[0]>=75)//broadleaved
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(11));
		else if(polyValues[1]>=75)//coniferous
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(12));
		else if(polyValues[0]>25&&polyValues[1]>25)//mixed
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(13));
		else{
		  if(verbose_opt[0]){
		    std::cout << "No valid value in polyValues..." << std::endl;
		    for(int index=0;index<polyValues.size();++index){
		      double theValue=polyValues[index];
		      std::cout << theValue << " ";
		    }
		    std::cout << std::endl;
		  }
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(20));
		}
		break;
	      }
	      case(rule::minimum):{
		//minimum of polygon
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//search for min class
		int minClass=stat.mymax(class_opt);
		for(int iclass=0;iclass<class_opt.size();++iclass){
		  if(polyValues[iclass]>0){
		    if(verbose_opt[0]>1)
		      std::cout << class_opt[iclass] << ": " << polyValues[iclass] << std::endl;
		    if(class_opt[iclass]<minClass)
		      minClass=class_opt[iclass];
		  }
		}
		if(verbose_opt[0]>0)
		  std::cout << "minClass: " << minClass << std::endl;
		writePolygonFeature->SetField(label_opt[0].c_str(),minClass);
		break;
	      }
	      case(rule::maximum):{
		//maximum of polygon
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//search for max class
		int maxClass=stat.mymin(class_opt);
		for(int iclass=0;iclass<class_opt.size();++iclass){
		  if(polyValues[iclass]>0){
		    if(verbose_opt[0]>1)
		      std::cout << class_opt[iclass] << ": " << polyValues[iclass] << std::endl;
		    if(class_opt[iclass]>maxClass)
		      maxClass=class_opt[iclass];
		  }
		}
		if(verbose_opt[0]>0)
		  std::cout << "maxClass: " << maxClass << std::endl;
		writePolygonFeature->SetField(label_opt[0].c_str(),maxClass);
		break;
	      }
	      case(rule::maxvote):{
		//maximum votes in polygon
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//search for class with maximum votes
		int maxClass=stat.mymin(class_opt);
		vector<double>::iterator maxit;
		maxit=stat.mymax(polyValues,polyValues.begin(),polyValues.end());
		int maxIndex=distance(polyValues.begin(),maxit);
		maxClass=class_opt[maxIndex];
		if(verbose_opt[0]>0)
		  std::cout << "maxClass: " << maxClass << std::endl;
		writePolygonFeature->SetField(label_opt[0].c_str(),maxClass);
		break;
	      }
	      }
	      if(polygon_opt[0]){
		if(verbose_opt[0]>1)
		  std::cout << "creating polygon feature" << std::endl;
		if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		  std::string errorString="Failed to create polygon feature in shapefile";
		  throw(errorString);
		}
		OGRFeature::DestroyFeature( writePolygonFeature );
		++ntotalvalid;
		if(verbose_opt[0])
		  std::cout << "ntotalvalid(1): " << ntotalvalid << std::endl;
	      }
	      else{
		if(verbose_opt[0]>1)
		  std::cout << "creating point feature in centroid" << std::endl;
		if(writeLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		  std::string errorString="Failed to create point feature in shapefile";
		  throw(errorString);
		}
		OGRFeature::DestroyFeature( writeCentroidFeature );
		++ntotalvalid;
		if(verbose_opt[0])
		  std::cout << "ntotalvalid: " << ntotalvalid << std::endl;
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
	      std::cout << "bounding box for feature " << ifeature << ": " << uli << " " << ulj << " " << lri << " " << lrj << std::endl;

	    if(uli<0||lri>=imgReader.nrOfCol()||ulj<0||ulj>=imgReader.nrOfRow())
	      continue;

	    int nPointPolygon=0;
	    if(polygon_opt[0]){
	      if(writeTest)
		writePolygonFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writePolygonFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    else if(ruleMap[rule_opt[0]]==rule::mean||ruleMap[rule_opt[0]]==rule::centroid||ruleMap[rule_opt[0]]==rule::sum){
	      if(writeTest)
		writeCentroidFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
	      else
		writeCentroidFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
	    }
	    //previously here
	    vector<double> polyValues;
	    switch(ruleMap[rule_opt[0]]){
	    case(rule::point):
	    case(rule::mean):
	    case(rule::sum):
	    case(rule::centroid):
	    default:
	      polyValues.resize(nband);
	    break;
	    case(rule::proportion):
	    case(rule::custom):
	    case(rule::minimum):
	    case(rule::maximum):
	    case(rule::maxvote):
	      assert(class_opt.size());
	    polyValues.resize(class_opt.size());
	    break;
	    }
	    for(int index=0;index<polyValues.size();++index)
	      polyValues[index]=0;
	    OGRPoint thePoint;
	    for(int j=ulj;j<=lrj;++j){
	      for(int i=uli;i<=lri;++i){
		//check if point is on surface
		double x=0;
		double y=0;
		imgReader.image2geo(i,j,x,y);
		thePoint.setX(x);
		thePoint.setY(y);
		if(readPolygon.Contains(&thePoint)){
		  bool valid=true;
		  for(int imask=0;imask<mask_opt.size();++imask){
		    double colMask,rowMask;//image coordinates in mask image
		    if(mask_opt.size()>1){//multiple masks
		      maskReader[imask].geo2image(x,y,colMask,rowMask);
		      //nearest neighbour
		      rowMask=static_cast<int>(rowMask);
		      colMask=static_cast<int>(colMask);
		      if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[imask].nrOfCol())
			continue;
		      // {
		      //   cerr << colMask << " out of mask col range!" << std::endl;
		      //   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
		      //   assert(static_cast<int>(colMask)>=0&&static_cast<int>(colMask)<maskReader[imask].nrOfCol());
		      // }
              
		      if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[imask])){
			if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[imask].nrOfRow())
			  continue;
			// {
			//   cerr << rowMask << " out of mask row range!" << std::endl;
			//   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
			//   assert(static_cast<int>(rowMask)>=0&&static_cast<int>(rowMask)<imgReader.nrOfRow());
			// }
			else{
			  maskReader[imask].readData(maskBuffer[imask],GDT_Int32,static_cast<int>(rowMask));
			  oldmaskrow[imask]=rowMask;
			  assert(maskBuffer.size()==maskReader[imask].nrOfBand());
			}
		      }
		      //               char ivalue=0;
		      int ivalue=0;
		      if(mask_opt.size()==msknodata_opt.size())//one invalid value for each mask
			ivalue=static_cast<int>(msknodata_opt[imask]);
		      else//use same invalid value for each mask
			ivalue=static_cast<int>(msknodata_opt[0]);
		      if(maskBuffer[imask][colMask]==ivalue){
			valid=false;
			break;
		      }
		    }
		    else if(maskReader.size()){
		      maskReader[0].geo2image(x,y,colMask,rowMask);
		      //nearest neighbour
		      rowMask=static_cast<int>(rowMask);
		      colMask=static_cast<int>(colMask);
		      if(static_cast<int>(colMask)<0||static_cast<int>(colMask)>=maskReader[0].nrOfCol())
			continue;
		      // {
		      //   cerr << colMask << " out of mask col range!" << std::endl;
		      //   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
		      //   assert(static_cast<int>(colMask)>=0&&static_cast<int>(colMask)<maskReader[0].nrOfCol());
		      // }
              
		      if(static_cast<int>(rowMask)!=static_cast<int>(oldmaskrow[0])){
			if(static_cast<int>(rowMask)<0||static_cast<int>(rowMask)>=maskReader[0].nrOfRow())
			  continue;
			// {
			//   cerr << rowMask << " out of mask row range!" << std::endl;
			//   cerr << x << " " << y << " " << colMask << " " << rowMask << std::endl;
			//   assert(static_cast<int>(rowMask)>=0&&static_cast<int>(rowMask)<imgReader.nrOfRow());
			// }
			else{
			  maskReader[0].readData(maskBuffer[0],GDT_Int32,static_cast<int>(rowMask));
			  oldmaskrow[0]=rowMask;
			}
		      }
		      for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
			if(maskBuffer[0][colMask]==static_cast<int>(msknodata_opt[ivalue])){
			  valid=false;
			  break;
			}
		      }
		    }
		  }
		  if(!valid)
		    continue;
		  else
		    validFeature=true;
		  //check if within raster image
		  if(i<0||i>=imgReader.nrOfCol())
		    continue;
		  if(j<0||j>=imgReader.nrOfRow())
		    continue;
		  writeRing.addPoint(&thePoint);
		  if(verbose_opt[0]>1)
		    std::cout << "point is on surface:" << thePoint.getX() << "," << thePoint.getY() << std::endl;
		  ++nPointPolygon;
		  OGRFeature *writePointFeature;
		  if(!polygon_opt[0]){
		    //create feature
		    if(ruleMap[rule_opt[0]]!=rule::mean&&ruleMap[rule_opt[0]]!=rule::centroid&&ruleMap[rule_opt[0]]!=rule::sum){//do not create in case of mean or sum (only create point at centroid)
		      if(writeTest)
			writePointFeature = OGRFeature::CreateFeature(writeTestLayer->GetLayerDefn());
		      else
			writePointFeature = OGRFeature::CreateFeature(writeLayer->GetLayerDefn());
		      if(verbose_opt[0]>1)
			std::cout << "copying fields from polygons " << sample_opt[0] << std::endl;
		      if(writePointFeature->SetFrom(readFeature)!= OGRERR_NONE)
			cerr << "writing feature failed" << std::endl;
		      writePointFeature->SetGeometry(&thePoint);
		      OGRGeometry *updateGeometry;
		      updateGeometry = writePointFeature->GetGeometryRef();
		      OGRPoint *poPoint = (OGRPoint *) updateGeometry;
		      if(verbose_opt[0]>1)
			std::cout << "write feature has " << writePointFeature->GetFieldCount() << " fields" << std::endl;
		    }
		  }
		  if(verbose_opt[0]>1)
		    std::cout << "reading image value withinin polygon at position " << i << "," << j;
		  for(int iband=0;iband<nband;++iband){
		    int theBand=(band_opt[0]<0)?iband:band_opt[iband];
		    double value=0;
		    imgReader.readData(value,GDT_Float64,i,j,theBand);
		    if(verbose_opt[0]>1)
		      std::cout << ": " << value << std::endl;
		    if(polygon_opt[0]||ruleMap[rule_opt[0]]==rule::mean||ruleMap[rule_opt[0]]==rule::centroid||ruleMap[rule_opt[0]]==rule::sum){
		      int iclass=0;
		      switch(ruleMap[rule_opt[0]]){
		      case(rule::point)://in centroid if polygon_opt==true or all values as points if polygon_opt!=true
		      case(rule::centroid):
		      default:
			polyValues[iband]=value;
		      break;
		      case(rule::sum):
		      case(rule::mean)://mean or sum polygon if polygon_opt==true or as point in centroid if polygon_opt!=true
			polyValues[iband]+=value;
			break;
		      case(rule::proportion):
		      case(rule::custom):
		      case(rule::minimum):
		      case(rule::maximum):
		      case(rule::maxvote):
			for(iclass=0;iclass<class_opt.size();++iclass){
			  if(value==class_opt[iclass]){
			    assert(polyValues.size()>iclass);
			    polyValues[iclass]+=1;//value
			    break;
			  }
			}
		      break;
		      }
		    }
		    else{
		      ostringstream fs;
		      // if(imgReader.nrOfBand()==1)
		      //   fs << fieldname_opt[0];
		      // else
		      //   fs << fieldname_opt[0] << theBand;
		      if(verbose_opt[0]>1)
			std::cout << "set field " << fieldname_opt[iband] << " to " << value << std::endl;
		      switch( fieldType ){
		      case OFTInteger:
			writePointFeature->SetField(fieldname_opt[iband].c_str(),static_cast<int>(value));
			break;
		      case OFTString:
			{
			  ostringstream os;
			  os << value;
			  writePointFeature->SetField(fieldname_opt[iband].c_str(),os.str().c_str());
			  break;
			}
		      case OFTReal:
			writePointFeature->SetField(fieldname_opt[iband].c_str(),value);
			break;
		      case OFTRealList:{
			int fieldIndex=writePointFeature->GetFieldIndex(fieldname_opt[iband].c_str());
			int nCount;
			const double *theList;
			theList=writePointFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
			vector<double> vectorList(nCount+1);
			for(int index=0;index<nCount;++index)
			  vectorList[nCount]=value;
			writePointFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
			break;
		      }
		      default://not supported
			assert(0);
			break;
		      }
		    }
		  }
		  if(!polygon_opt[0]){
		    if(ruleMap[rule_opt[0]]!=rule::mean&&ruleMap[rule_opt[0]]!=rule::centroid&&ruleMap[rule_opt[0]]!=rule::sum){//do not create in case of mean value (only at centroid)
		      //write feature
		      if(verbose_opt[0]>1)
			std::cout << "creating point feature" << std::endl;
		      if(writeTest){
			if(writeTestLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			  std::string errorString="Failed to create feature in shapefile";
			  throw(errorString);
			}
		      }
		      else{
			if(writeLayer->CreateFeature( writePointFeature ) != OGRERR_NONE ){
			  std::string errorString="Failed to create feature in shapefile";
			  throw(errorString);
			}
		      }
		      //destroy feature
		      OGRFeature::DestroyFeature( writePointFeature );
		    }
		  }
		  // ++isample;
		  ++ntotalvalid;
		  if(verbose_opt[0])
		    std::cout << "ntotalvalid: " << ntotalvalid << std::endl;
		}
	      }
	    }
	    //test
	    if(!validFeature)
	      continue;
	    if(polygon_opt[0]||ruleMap[rule_opt[0]]==rule::mean||ruleMap[rule_opt[0]]==rule::centroid||ruleMap[rule_opt[0]]==rule::sum){
	      //add ring to polygon
	      if(polygon_opt[0]){
		writePolygon.addRing(&writeRing);
		writePolygon.closeRings();
		//write geometry of writePolygon
		writePolygonFeature->SetGeometry(&writePolygon);
		if(writePolygonFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "copying new fields write polygon " << sample_opt[0] << std::endl;
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writePolygonFeature->GetFieldCount() << " fields" << std::endl;
		//write polygon feature
	      }
	      else{//write mean value of polygon to centroid point (ruleMap[rule_opt[0]]==rule::mean)
		//create feature
		if(verbose_opt[0]>1)
		  std::cout << "copying fields from polygons " << sample_opt[0] << std::endl;
		if(writeCentroidFeature->SetFrom(readFeature)!= OGRERR_NONE)
		  cerr << "writing feature failed" << std::endl;
		writeCentroidFeature->SetGeometry(&writeCentroidPoint);
		OGRGeometry *updateGeometry;
		updateGeometry = writeCentroidFeature->GetGeometryRef();
		assert(wkbFlatten(updateGeometry->getGeometryType()) == wkbPoint );
		if(verbose_opt[0]>1)
		  std::cout << "write feature has " << writeCentroidFeature->GetFieldCount() << " fields" << std::endl;
	      }
	      switch(ruleMap[rule_opt[0]]){
	      case(rule::point)://value at each point (or at centroid of polygon if line is set)
	      default:{
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		for(int index=0;index<polyValues.size();++index){
		  double theValue=polyValues[index];
		  ostringstream fs;
		  if(verbose_opt[0])
		    std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		  int theBand=(band_opt[0]<0)?index:band_opt[index];
		  // if(nband==1)
		  //   fs << fieldname_opt[0];
		  // else
		  //   fs << fieldname_opt[0] << theBand;
		  if(verbose_opt[0]>1)
		    std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		  switch( fieldType ){
		  case OFTInteger:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    break;
		  case OFTString:
		    {
		      ostringstream os;
		      os << theValue;
		      if(polygon_opt[0])
			writePolygonFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      else
			writeCentroidFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      break;
		    }
		  case OFTReal:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    break;
		  case OFTRealList:{
		    int fieldIndex;
		    int nCount;
		    const double *theList;
		    vector<double> vectorList;
		    if(polygon_opt[0]){
		      fieldIndex=writePolygonFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writePolygonFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writePolygonFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    else{
		      fieldIndex=writeCentroidFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writeCentroidFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writeCentroidFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    break;
		  }//case OFTRealList
		  }//switch(fieldType)
		}//for index
		break;
	      }//case 0 and default
	      case(rule::mean):
	      case(rule::sum):
	      case(rule::centroid):{//mean value (written to centroid of polygon if line is not set)
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//test
		if(ruleMap[rule_opt[0]]==rule::centroid)
		  assert(nPointPolygon<=1);
		for(int index=0;index<polyValues.size();++index){
		  double theValue=polyValues[index];
		  // ostringstream fs;
		  if(ruleMap[rule_opt[0]]==rule::mean)
		    theValue/=nPointPolygon;
		  int theBand=(band_opt[0]<0)?index:band_opt[index];
		  // if(nband==1)
		  //   fs << fieldname_opt[0];
		  // else
		  //   fs << fieldname_opt[0] << theBand;
		  if(verbose_opt[0]>1)
		    std::cout << "set field " << fieldname_opt[index] << " to " << theValue << std::endl;
		  switch( fieldType ){
		  case OFTInteger:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),static_cast<int>(theValue));
		    break;
		  case OFTString:
		    {
		      ostringstream os;
		      os << theValue;
		      if(polygon_opt[0])
			writePolygonFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      else
			writeCentroidFeature->SetField(fieldname_opt[index].c_str(),os.str().c_str());
		      break;
		    }
		  case OFTReal:
		    if(polygon_opt[0])
		      writePolygonFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    else
		      writeCentroidFeature->SetField(fieldname_opt[index].c_str(),theValue);
		    break;
		  case OFTRealList:{
		    int fieldIndex;
		    int nCount;
		    const double *theList;
		    vector<double> vectorList;
		    if(polygon_opt[0]){
		      fieldIndex=writePolygonFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writePolygonFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writePolygonFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    else{
		      fieldIndex=writeCentroidFeature->GetFieldIndex(fieldname_opt[index].c_str());
		      theList=writeCentroidFeature->GetFieldAsDoubleList(fieldIndex,&nCount);
		      vectorList.resize(nCount+1);
		      for(int index=0;index<nCount;++index)
			vectorList[index]=theList[index];
		      vectorList[nCount]=theValue;
		      writeCentroidFeature->SetField(fieldIndex,vectorList.size(),&(vectorList[0]));
		    }
		    break;
		  }
		  }
		}
		break;
	      }
	      case(rule::proportion):{//proportion classes
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		stat.normalize_pct(polyValues);
		// stat.sum(polyValues);
		for(int index=0;index<polyValues.size();++index){
		  double theValue=polyValues[index];
		  ostringstream fs;
		  fs << class_opt[index];
		  if(polygon_opt[0])
		    writePolygonFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		  else
		    writeCentroidFeature->SetField(fs.str().c_str(),static_cast<int>(theValue));
		}
		break;
	      }
	      case(rule::custom):{//custom
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		stat.normalize_pct(polyValues);
		assert(polyValues.size()==2);//11:broadleaved, 12:coniferous
		if(polyValues[0]>=75)//broadleaved
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(11));
		else if(polyValues[1]>=75)//coniferous
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(12));
		else if(polyValues[0]>25&&polyValues[1]>25)//mixed
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(13));
		else{
		  if(verbose_opt[0]){
		    std::cout << "No valid value in polyValues..." << std::endl;
		    for(int index=0;index<polyValues.size();++index){
		      double theValue=polyValues[index];
		      std::cout << theValue << " ";
		    }
		    std::cout << std::endl;
		  }
		  writePolygonFeature->SetField(label_opt[0].c_str(),static_cast<int>(20));
		}
		break;
	      }
	      case(rule::minimum):{//minimum of polygon
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//search for min class
		int minClass=stat.mymax(class_opt);
		for(int iclass=0;iclass<class_opt.size();++iclass){
		  if(polyValues[iclass]>0){
		    if(verbose_opt[0]>1)
		      std::cout << class_opt[iclass] << ": " << polyValues[iclass] << std::endl;
		    if(class_opt[iclass]<minClass)
		      minClass=class_opt[iclass];
		  }
		}
		if(verbose_opt[0]>0)
		  std::cout << "minClass: " << minClass << std::endl;
		writePolygonFeature->SetField(label_opt[0].c_str(),minClass);
		break;
	      }
	      case(rule::maximum):{//maximum of polygon
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//search for max class
		int maxClass=stat.mymin(class_opt);
		for(int iclass=0;iclass<class_opt.size();++iclass){
		  if(polyValues[iclass]>0){
		    if(verbose_opt[0]>1)
		      std::cout << class_opt[iclass] << ": " << polyValues[iclass] << std::endl;
		    if(class_opt[iclass]>maxClass)
		      maxClass=class_opt[iclass];
		  }
		}
		if(verbose_opt[0]>0)
		  std::cout << "maxClass: " << maxClass << std::endl;
		writePolygonFeature->SetField(label_opt[0].c_str(),maxClass);
		break;
	      }
	      case(rule::maxvote):{//maximum votes in polygon
		assert(polygon_opt[0]);//not implemented for points
		if(verbose_opt[0])
		  std::cout << "number of points in polygon: " << nPointPolygon << std::endl;
		//search for max votes
		int maxClass=stat.mymin(class_opt);
		vector<double>::iterator maxit;
		maxit=stat.mymax(polyValues,polyValues.begin(),polyValues.end());
		int maxIndex=distance(polyValues.begin(),maxit);
		maxClass=class_opt[maxIndex];
	      }
	      }
	      if(polygon_opt[0]){
		if(verbose_opt[0]>1)
		  std::cout << "creating polygon feature" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in shapefile";
		    throw(errorString);
		  }
		}
		else{
		  if(writeLayer->CreateFeature( writePolygonFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create polygon feature in shapefile";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writePolygonFeature );
		++ntotalvalid;
		if(verbose_opt[0])
		  std::cout << "ntotalvalid: " << ntotalvalid << std::endl;
	      }
	      else{
		if(verbose_opt[0]>1)
		  std::cout << "creating point feature in centroid" << std::endl;
		if(writeTest){
		  if(writeTestLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in shapefile";
		    throw(errorString);
		  }
		}
		else{
		  //test
		  assert(validFeature);
		  if(writeLayer->CreateFeature( writeCentroidFeature ) != OGRERR_NONE ){
		    std::string errorString="Failed to create point feature in shapefile";
		    throw(errorString);
		  }
		}
		OGRFeature::DestroyFeature( writeCentroidFeature );
		++ntotalvalid;
		if(verbose_opt[0])
		  std::cout << "ntotalvalid: " << ntotalvalid << std::endl;
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
	  progress=static_cast<float>(ifeature+1)/nfeature;
	  pfnProgress(progress,pszMessage,pProgressArg);
	}
	catch(std::string e){
	  std::cout << e << std::endl;
	  continue;
	}
      }//end of getNextFeature
      // if(rbox_opt[0]>0||cbox_opt[0]>0)
      //   boxWriter.close();
      progress=1.0;
      pfnProgress(progress,pszMessage,pProgressArg);
    }
    ogrWriter.close();
    if(test_opt.size())
      ogrTestWriter.close();
  }
  progress=1.0;
  pfnProgress(progress,pszMessage,pProgressArg);
  imgReader.close();
}
  

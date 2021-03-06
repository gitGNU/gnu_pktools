/**********************************************************************
pkextractimg_lib.cc: extract pixel values from raster image using a raster sample
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
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

CPLErr ImgRasterGdal::extractImg(const AppFactory& app){
  // Optionpk<string> image_opt("i", "input", "Raster input dataset containing band information");
  Optionpk<string> sample_opt("s", "sample", "OGR vector dataset with features to be extracted from input data. Output will contain features with input band information included. Sample image can also be GDAL raster dataset.");
  Optionpk<string> output_opt("o", "output", "Output sample dataset");
  Optionpk<int> class_opt("c", "class", "Class(es) to extract from input sample image. Leave empty to extract all valid data pixels from sample dataset");
  Optionpk<float> threshold_opt("t", "threshold", "Probability threshold for selecting samples (randomly). Provide probability in percentage (>0) or absolute (<0). Use a single threshold per vector sample layer. If using raster land cover maps as a sample dataset, you can provide a threshold value for each class (e.g. -t 80 -t 60). Use value 100 to select all pixels for selected class(es)", 100);
  Optionpk<string> ogrformat_opt("f", "f", "Output sample dataset format","SQLite");
  Optionpk<string> ftype_opt("ft", "ftype", "Field type (only Real or Integer)", "Real");
  Optionpk<string> ltype_opt("lt", "ltype", "Label type: In16 or String", "Integer");
  Optionpk<unsigned int> band_opt("b", "band", "Band index(es) to extract (0 based). Leave empty to use all bands");
  Optionpk<unsigned short> bstart_opt("sband", "startband", "Start band sequence number"); 
  Optionpk<unsigned short> bend_opt("eband", "endband", "End band sequence number"); 
  Optionpk<double> srcnodata_opt("srcnodata", "srcnodata", "Invalid value(s) for input image");
  Optionpk<unsigned int> bndnodata_opt("bndnodata", "bndnodata", "Band in input image to check if pixel is valid (used for srcnodata)", 0);
  Optionpk<string> fieldname_opt("bn", "bname", "For single band input data, this extra attribute name will correspond to the raster values. For multi-band input data, multiple attributes with this prefix will be added (e.g. b0, b1, b2, etc.)", "b");
  Optionpk<string> label_opt("cn", "cname", "Name of the class label in the output vector dataset", "label");
  Optionpk<short> down_opt("down", "down", "Down sampling factor", 1);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose mode if > 0", 0,2);

  bstart_opt.setHide(1);
  bend_opt.setHide(1);
  bndnodata_opt.setHide(1);
  srcnodata_opt.setHide(1);
  fieldname_opt.setHide(1);
  label_opt.setHide(1);
  down_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=sample_opt.retrieveOption(app.getArgc(),app.getArgv());
    output_opt.retrieveOption(app.getArgc(),app.getArgv());
    class_opt.retrieveOption(app.getArgc(),app.getArgv());
    threshold_opt.retrieveOption(app.getArgc(),app.getArgv());
    ogrformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    ftype_opt.retrieveOption(app.getArgc(),app.getArgv());
    ltype_opt.retrieveOption(app.getArgc(),app.getArgv());
    band_opt.retrieveOption(app.getArgc(),app.getArgv());
    bstart_opt.retrieveOption(app.getArgc(),app.getArgv());
    bend_opt.retrieveOption(app.getArgc(),app.getArgv());
    bndnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    srcnodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    fieldname_opt.retrieveOption(app.getArgc(),app.getArgv());
    label_opt.retrieveOption(app.getArgc(),app.getArgv());
    down_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
    }

    statfactory::StatFactory stat;
    if(srcnodata_opt.size()){
      while(srcnodata_opt.size()<bndnodata_opt.size())
        srcnodata_opt.push_back(srcnodata_opt[0]);
      stat.setNoDataValues(srcnodata_opt);
    }
    Vector2d<unsigned int> posdata;
    unsigned long int nsample=0;
    unsigned long int ntotalvalid=0;
    unsigned long int ntotalinvalid=0;

    map<int,unsigned long int> nvalid;
    map<int,unsigned long int> ninvalid;
    // vector<unsigned long int> nvalid(class_opt.size());
    // vector<unsigned long int> ninvalid(class_opt.size());
    // if(class_opt.empty()){
    //   nvalid.resize(256);
    //   ninvalid.resize(256);
    // }
    // for(int it=0;it<nvalid.size();++it){
    //   nvalid[it]=0;
    //   ninvalid[it]=0;
    // }

    map <int,short> classmap;//class->index
    for(int iclass=0;iclass<class_opt.size();++iclass){
      nvalid[class_opt[iclass]]=0;
      ninvalid[class_opt[iclass]]=0;
      classmap[class_opt[iclass]]=iclass;
    }

    // ImgRasterGdal imgReader;
    // if(image_opt.empty()){
    //   std::cerr << "No image dataset provided (use option -i). Use --help for help information";
    //     exit(0);
    // }
    if(output_opt.empty()){
      std::cerr << "No output dataset provided (use option -o). Use --help for help information";
      return(CE_Failure);
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
      std::cout << "Number of bands in input image: " << this->nrOfBand() << std::endl;

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

    bool sampleIsRaster=true;

    ImgRasterGdal classReader;
    ImgWriterOgr sampleWriterOgr;

    if(sampleIsRaster){
      if(class_opt.empty()){
        ImgWriterOgr ogrWriter;
        assert(sample_opt.size());
        classReader.open(sample_opt[0]);
        // vector<int> classBuffer(classReader.nrOfCol());
        vector<double> classBuffer(classReader.nrOfCol());
        Vector2d<double> imgBuffer(nband,this->nrOfCol());//[band][col]
        // vector<double> imgBuffer(nband);//[band]
        vector<double> sample(2+nband);//x,y,band values
        Vector2d<double> writeBuffer;
        vector<int> writeBufferClass;
        vector<int> selectedClass;
        Vector2d<double> selectedBuffer;
        unsigned int irow=0;
        unsigned int icol=0;
        if(verbose_opt[0]>1)
          std::cout << "extracting sample from image..." << std::endl;
        progress=0;
        pfnProgress(progress,pszMessage,pProgressArg);
        for(irow=0;irow<classReader.nrOfRow();++irow){
          if(irow%down_opt[0])
            continue;
          classReader.readData(classBuffer,irow);
          double x=0;//geo x coordinate
          double y=0;//geo y coordinate
          double iimg=0;//image x-coordinate in img image
          double jimg=0;//image y-coordinate in img image

          //find col in img
          classReader.image2geo(icol,irow,x,y);
          this->geo2image(x,y,iimg,jimg);
          //nearest neighbour
          if(static_cast<unsigned int>(jimg)<0||static_cast<unsigned int>(jimg)>=this->nrOfRow())
            continue;
          for(unsigned int iband=0;iband<nband;++iband){
            unsigned int theBand=(band_opt.size()) ? band_opt[iband] : iband;
            this->readData(imgBuffer[iband],static_cast<unsigned int>(jimg),theBand);
          }
          for(icol=0;icol<classReader.nrOfCol();++icol){
            if(icol%down_opt[0])
              continue;
            unsigned int theClass=classBuffer[icol];
            unsigned int processClass=0;
            bool valid=false;
            if(class_opt.empty()){
              valid=true;//process every class
              processClass=theClass;
            }
            else{
              for(unsigned int iclass=0;iclass<class_opt.size();++iclass){
                if(classBuffer[icol]==class_opt[iclass]){
                  processClass=iclass;
                  theClass=class_opt[iclass];
                  valid=true;//process this class
                  break;
                }
              }
            }
            classReader.image2geo(icol,irow,x,y);
            sample[0]=x;
            sample[1]=y;
            if(verbose_opt[0]>1){
              std::cout.precision(12);
              std::cout << theClass << " " << x << " " << y << std::endl;
            }
            //find col in img
            this->geo2image(x,y,iimg,jimg);
            //nearest neighbour
            iimg=static_cast<unsigned int>(iimg);
            if(static_cast<unsigned int>(iimg)<0||static_cast<unsigned int>(iimg)>=this->nrOfCol())
              continue;

            for(unsigned int iband=0;iband<nband&&valid;++iband){
              unsigned int theBand=(band_opt.size()) ? band_opt[iband] : iband;
              if(srcnodata_opt.size()&&theBand==bndnodata_opt[0]){
                // vector<unsigned int>::const_iterator bndit=bndnodata_opt.begin();
                for(unsigned int inodata=0;inodata<srcnodata_opt.size()&&valid;++inodata){
                  if(imgBuffer[iband][iimg]==srcnodata_opt[inodata])
                    valid=false;
                }
              }
            }
            // oldimgrow=jimg;

            if(valid){
              for(unsigned int iband=0;iband<imgBuffer.size();++iband){
                sample[iband+2]=imgBuffer[iband][iimg];
              }
              float theThreshold=(threshold_opt.size()>1)?threshold_opt[processClass]:threshold_opt[0];
              if(theThreshold>0){//percentual value
                double p=static_cast<double>(rand())/(RAND_MAX);
                p*=100.0;
                if(p>theThreshold)
                  continue;//do not select for now, go to next column
              }
              // else if(nvalid.size()>processClass){//absolute value
              //   if(nvalid[processClass]>=-theThreshold)
              //     continue;//do not select any more pixels for this class, go to next column to search for other classes
              // }
              writeBuffer.push_back(sample);
              writeBufferClass.push_back(theClass);
              ++ntotalvalid;
              if(nvalid.count(theClass))
                nvalid[theClass]+=1;
              else
                nvalid[theClass]=1;
            }
            else{
              ++ntotalinvalid;
              if(ninvalid.count(theClass))
                ninvalid[theClass]+=1;
              else
                ninvalid[theClass]=1;
            }
          }
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
          ogrWriter.createLayer(layername, this->getProjection(), wkbPoint, papszOptions);
          std::string fieldname="fid";//number of the point
          ogrWriter.createField(fieldname,OFTInteger);
          map<std::string,double> pointAttributes;
          ogrWriter.createField(label_opt[0],labelType);
          for(unsigned int iband=0;iband<nband;++iband){
            ogrWriter.createField(fieldname_opt[iband],fieldType);
          }
          progress=0;
          pfnProgress(progress,pszMessage,pProgressArg);
        
          map<int,short> classDone;
          Vector2d<double> writeBufferTmp;
          vector<int> writeBufferClassTmp;

          if(threshold_opt[0]<0){//absolute threshold
            map<int,unsigned long int>::iterator mapit;
            map<int,unsigned long int> ncopied;
            for(mapit=nvalid.begin();mapit!=nvalid.end();++mapit)
              ncopied[mapit->first]=0;

            cout << "ncopied.size(): " << ncopied.size() << endl;
            while(classDone.size()<nvalid.size()){
              unsigned int index=rand()%writeBufferClass.size();
              unsigned int theClass=writeBufferClass[index];
              float theThreshold=threshold_opt[0];
              if(threshold_opt.size()>1&&class_opt.size())
                theThreshold=threshold_opt[classmap[theClass]];
              theThreshold=-theThreshold;
              if(ncopied[theClass]<theThreshold){
                writeBufferClassTmp.push_back(*(writeBufferClass.begin()+index));
                writeBufferTmp.push_back(*(writeBuffer.begin()+index));
                writeBufferClass.erase(writeBufferClass.begin()+index);
                writeBuffer.erase(writeBuffer.begin()+index);
                ++(ncopied[theClass]);
              }
              else
                classDone[theClass]=1;
              if(ncopied[theClass]>=nvalid[theClass]){
                classDone[theClass]=1;
              }
            }
            writeBuffer=writeBufferTmp;
            writeBufferClass=writeBufferClassTmp;

            //   while(classDone.size()<nvalid.size()){
            //     unsigned int index=rand()%writeBufferClass.size();
            //     unsigned int theClass=writeBufferClass[index];
            //     float theThreshold=threshold_opt[0];
            //     if(threshold_opt.size()>1&&class_opt.size())
            //       theThreshold=threshold_opt[classmap[theClass]];
            //     theThreshold=-theThreshold;
            //     if(nvalid[theClass]>theThreshold){
            //       writeBufferClass.erase(writeBufferClass.begin()+index);
            //       writeBuffer.erase(writeBuffer.begin()+index);
            //       --(nvalid[theClass]);
            //     }
            //     else
            //       classDone[theClass]=1;
            //   }
          }
          for(unsigned int isample=0;isample<writeBuffer.size();++isample){
            if(verbose_opt[0]>1)
              std::cout << "writing sample " << isample << std::endl;
            pointAttributes[label_opt[0]]=writeBufferClass[isample];
            for(unsigned int iband=0;iband<writeBuffer[0].size()-2;++iband){
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
        ImgRasterGdal classReader;
        ImgWriterOgr ogrWriter;
        if(verbose_opt[0]>1){
          std::cout << "reading position from sample dataset " << std::endl;
          std::cout << "class thresholds: " << std::endl;
          for(unsigned int iclass=0;iclass<class_opt.size();++iclass){
            if(threshold_opt.size()>1)
              std::cout << class_opt[iclass] << ": " << threshold_opt[iclass] << std::endl;
            else
              std::cout << class_opt[iclass] << ": " << threshold_opt[0] << std::endl;
          }
        }
        classReader.open(sample_opt[0]);
        vector<int> classBuffer(classReader.nrOfCol());
        // vector<double> classBuffer(classReader.nrOfCol());
        Vector2d<double> imgBuffer(nband,this->nrOfCol());//[band][col]
        // vector<double> imgBuffer(nband);//[band]
        vector<double> sample(2+nband);//x,y,band values
        Vector2d<double> writeBuffer;
        vector<int> writeBufferClass;
        vector<int> selectedClass;
        Vector2d<double> selectedBuffer;
        unsigned int irow=0;
        unsigned int icol=0;
        if(verbose_opt[0]>1)
          std::cout << "extracting sample from image..." << std::endl;
        progress=0;
        pfnProgress(progress,pszMessage,pProgressArg);
        for(irow=0;irow<classReader.nrOfRow();++irow){
          if(irow%down_opt[0])
            continue;
          classReader.readData(classBuffer,irow);
          double x=0;//geo x coordinate
          double y=0;//geo y coordinate
          double iimg=0;//image x-coordinate in img image
          double jimg=0;//image y-coordinate in img image

          //find col in img
          classReader.image2geo(icol,irow,x,y);
          this->geo2image(x,y,iimg,jimg);
          //nearest neighbour
          if(static_cast<unsigned int>(jimg)<0||static_cast<unsigned int>(jimg)>=this->nrOfRow())
            continue;
          for(unsigned int iband=0;iband<nband;++iband){
            unsigned int theBand=(band_opt.size()) ? band_opt[iband] : iband;
            this->readData(imgBuffer[iband],static_cast<unsigned int>(jimg),theBand);
          }

          for(icol=0;icol<classReader.nrOfCol();++icol){
            if(icol%down_opt[0])
              continue;
            unsigned int theClass=0;
            // double theClass=0;
            unsigned int processClass=-1;
            if(class_opt.empty()){//process every class
              if(classBuffer[icol]){
                processClass=0;
                theClass=classBuffer[icol];
              }
            }
            else{
              for(unsigned int iclass=0;iclass<class_opt.size();++iclass){
                if(classBuffer[icol]==class_opt[iclass]){
                  processClass=iclass;
                  theClass=class_opt[iclass];
                }
              }
            }
            if(processClass>=0){
              //         if(classBuffer[icol]==class_opt[0]){
              classReader.image2geo(icol,irow,x,y);
              sample[0]=x;
              sample[1]=y;
              if(verbose_opt[0]>1){
                std::cout.precision(12);
                std::cout << theClass << " " << x << " " << y << std::endl;
              }
              //find col in img
              this->geo2image(x,y,iimg,jimg);
              //nearest neighbour
              iimg=static_cast<unsigned int>(iimg);
              if(static_cast<unsigned int>(iimg)<0||static_cast<unsigned int>(iimg)>=this->nrOfCol())
                continue;
              bool valid=true;

              for(unsigned int iband=0;iband<nband;++iband){
                unsigned int theBand=(band_opt.size()) ? band_opt[iband] : iband;
                if(srcnodata_opt.size()&&theBand==bndnodata_opt[0]){
                  // vector<int>::const_iterator bndit=bndnodata_opt.begin();
                  for(unsigned int inodata=0;inodata<srcnodata_opt.size()&&valid;++inodata){
                    if(imgBuffer[iband][iimg]==srcnodata_opt[inodata])
                      valid=false;
                  }
                }
              }
              if(valid){
                for(unsigned int iband=0;iband<imgBuffer.size();++iband){
                  sample[iband+2]=imgBuffer[iband][iimg];
                }
                float theThreshold=(threshold_opt.size()>1)?threshold_opt[processClass]:threshold_opt[0];
                if(theThreshold>0){//percentual value
                  double p=static_cast<double>(rand())/(RAND_MAX);
                  p*=100.0;
                  if(p>theThreshold)
                    continue;//do not select for now, go to next column
                }
                // else if(nvalid.size()>processClass){//absolute value
                //   if(nvalid[processClass]>=-theThreshold)
                //     continue;//do not select any more pixels for this class, go to next column to search for other classes
                // }
                writeBuffer.push_back(sample);
                writeBufferClass.push_back(theClass);
                ++ntotalvalid;
                if(nvalid.count(theClass))
                  nvalid[theClass]+=1;
                else
                  nvalid[theClass]=1;
              }
              else{
                ++ntotalinvalid;
                if(ninvalid.count(theClass))
                  ninvalid[theClass]+=1;
                else
                  ninvalid[theClass]=1;
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
          char **papszOptions=NULL;
          ostringstream slayer;
          slayer << "training data";
          std::string layername=slayer.str();
          ogrWriter.createLayer(layername, this->getProjection(), wkbPoint, papszOptions);
          std::string fieldname="fid";//number of the point
          ogrWriter.createField(fieldname,OFTInteger);
          map<std::string,double> pointAttributes;
          ogrWriter.createField(label_opt[0],labelType);
          for(unsigned int iband=0;iband<nband;++iband){
            ogrWriter.createField(fieldname_opt[iband],fieldType);
          }
          pfnProgress(progress,pszMessage,pProgressArg);
          progress=0;
          pfnProgress(progress,pszMessage,pProgressArg);

          map<int,short> classDone;
          Vector2d<double> writeBufferTmp;
          vector<int> writeBufferClassTmp;

          if(threshold_opt[0]<0){//absolute threshold
            map<int,unsigned long int>::iterator mapit;
            map<int,unsigned long int> ncopied;
            for(mapit=nvalid.begin();mapit!=nvalid.end();++mapit)
              ncopied[mapit->first]=0;

            while(classDone.size()<nvalid.size()){
              unsigned int index=rand()%writeBufferClass.size();
              unsigned int theClass=writeBufferClass[index];
              float theThreshold=threshold_opt[0];
              if(threshold_opt.size()>1&&class_opt.size())
                theThreshold=threshold_opt[classmap[theClass]];
              theThreshold=-theThreshold;
              if(ncopied[theClass]<theThreshold){
                writeBufferClassTmp.push_back(*(writeBufferClass.begin()+index));
                writeBufferTmp.push_back(*(writeBuffer.begin()+index));
                writeBufferClass.erase(writeBufferClass.begin()+index);
                writeBuffer.erase(writeBuffer.begin()+index);
                ++(ncopied[theClass]);
              }
              else
                classDone[theClass]=1;
              if(ncopied[theClass]>=nvalid[theClass]){
                classDone[theClass]=1;
              }
            }
            writeBuffer=writeBufferTmp;
            writeBufferClass=writeBufferClassTmp;
            // while(classDone.size()<nvalid.size()){
            //   unsigned int index=rand()%writeBufferClass.size();
            //   unsigned int theClass=writeBufferClass[index];
            //   float theThreshold=threshold_opt[0];
            //   if(threshold_opt.size()>1&&class_opt.size())
            //     theThreshold=threshold_opt[classmap[theClass]];
            //   theThreshold=-theThreshold;
            //   if(nvalid[theClass]>theThreshold){
            //     writeBufferClass.erase(writeBufferClass.begin()+index);
            //     writeBuffer.erase(writeBuffer.begin()+index);
            //     --(nvalid[theClass]);
            //   }
            //   else
            //     classDone[theClass]=1;
            // }
          }

          for(unsigned int isample=0;isample<writeBuffer.size();++isample){
            pointAttributes[label_opt[0]]=writeBufferClass[isample];
            for(unsigned int iband=0;iband<writeBuffer[0].size()-2;++iband){
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
            for(unsigned int iclass=0;iclass<class_opt.size();++iclass)
              std::cout << "class " << class_opt[iclass] << " has " << nvalid[iclass] << " samples" << std::endl;
          }
        }
      }
    }
    else{//vector dataset
      cerr << "Error: vector sample not supported, consider using pkextractogr" << endl;
    }//else (vector)
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
  

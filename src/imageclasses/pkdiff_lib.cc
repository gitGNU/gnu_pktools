/**********************************************************************
pkdiff_lib.cc: program to compare two raster image files
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
#include "imageclasses/ImgRaster.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "apps/AppFactory.h"
#include "algorithms/ConfusionMatrix.h"

using namespace std;
using namespace app;

/**
 * @param imgWriter output classified raster dataset
 * @param app application specific option arguments
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRaster::diff(shared_ptr<ImgRaster> imgReference, const AppFactory& app){
  Optionpk<string> mask_opt("m", "mask", "Use the first band of the specified file as a validity mask. Nodata values can be set with the option msknodata.");
  Optionpk<double> msknodata_opt("msknodata", "msknodata", "Mask value(s) where image is invalid. Use negative value for valid data (example: use -t -1: if only -1 is valid value)", 0);
  Optionpk<double> nodata_opt("nodata", "nodata", "No data value(s) in input or reference dataset are ignored");
  Optionpk<unsigned int> band_opt("b", "band", "Input (reference) raster band. Optionally, you can define different bands for input and reference bands respectively: -b 1 -b 0.", 0);
  Optionpk<bool> rmse_opt("rmse", "rmse", "Report root mean squared error", false);
  Optionpk<bool> regression_opt("reg", "reg", "Report linear regression (Input = c0+c1*Reference)", false);
  Optionpk<bool> confusion_opt("cm", "confusion", "Create confusion matrix (to std out)", false);
  Optionpk<string> cmformat_opt("cmf","cmf","Format for confusion matrix (ascii or latex)","ascii");
  Optionpk<string> cmoutput_opt("cmo","cmo","Output file for confusion matrix");
  Optionpk<bool> se95_opt("se95","se95","Report standard error for 95 confidence interval",false);
  Optionpk<string> classname_opt("c", "class", "List of class names.");
  Optionpk<short> classvalue_opt("r", "reclass", "List of class values (use same order as in classname option).");
  Optionpk<string> output_opt("o", "output", "Output dataset (optional)");
  Optionpk<string>  oformat_opt("of", "oformat", "Output image format (see also gdal_translate).","GTiff");
  Optionpk<string> colorTable_opt("ct", "ct", "Color table in ASCII format having 5 columns: id R G B ALFA (0: transparent, 255: solid).");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<short> valueE_opt("\0", "correct", "Value for correct pixels", 0,2);
  Optionpk<short> valueO_opt("\0", "omission", "Value for omission errors: input label > reference label", 1,2);
  Optionpk<short> valueC_opt("\0", "commission", "Value for commission errors: input label < reference label", 2,1);
  Optionpk<short> verbose_opt("v", "verbose", "Verbose level", 0,2);

  output_opt.setHide(1);
  oformat_opt.setHide(1);
  colorTable_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=band_opt.retrieveOption(app.getArgc(),app.getArgv());
    rmse_opt.retrieveOption(app.getArgc(),app.getArgv());
    regression_opt.retrieveOption(app.getArgc(),app.getArgv());
    confusion_opt.retrieveOption(app.getArgc(),app.getArgv());
    classname_opt.retrieveOption(app.getArgc(),app.getArgv());
    classvalue_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    mask_opt.retrieveOption(app.getArgc(),app.getArgv());
    msknodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    output_opt.retrieveOption(app.getArgc(),app.getArgv());
    cmformat_opt.retrieveOption(app.getArgc(),app.getArgv());
    cmoutput_opt.retrieveOption(app.getArgc(),app.getArgv());
    se95_opt.retrieveOption(app.getArgc(),app.getArgv());
    colorTable_opt.retrieveOption(app.getArgc(),app.getArgv());
    valueE_opt.retrieveOption(app.getArgc(),app.getArgv());
    valueO_opt.retrieveOption(app.getArgc(),app.getArgv());
    valueC_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
    if(!doProcess){
      cout << endl;
      std::ostringstream helpStream;
      helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
      throw(helpStream.str());//help was invoked, stop processing
    }
    ImgRaster maskReader;

    if(verbose_opt[0]){
      cout << "no data flag(s) set to";
      for(int iflag=0;iflag<nodata_opt.size();++iflag)
        cout << " " << nodata_opt[iflag];
      cout << endl;
    }

    //band_opt[0] is for input
    //band_opt[1] is for reference
    if(band_opt.size()<2)
      band_opt.push_back(band_opt[0]);

    // if(mask_opt.size())
    //   while(mask_opt.size()<input_opt.size())
    //     mask_opt.push_back(mask_opt[0]);
    vector<short> inputRange;
    vector<short> referenceRange;
    confusionmatrix::ConfusionMatrix cm;
    int nclass=0;
    map<string,short> classValueMap;
    vector<std::string> nameVector(255);//the inverse of the classValueMap
    vector<string> classNames;

    unsigned int ntotalValidation=0;
    unsigned int nflagged=0;
    Vector2d<unsigned int> resultClass;
    vector<float> user;
    vector<float> producer;
    vector<unsigned int> nvalidation;

    if(confusion_opt[0]){
      if(verbose_opt[0])
        cout << "opening input image file" << endl;
      // this->open(input_opt[0],memory_opt[0]);//,imagicX_opt[0],imagicY_opt[0]);
      this->getRange(inputRange,band_opt[0]);
      this->close();

      for(int iflag=0;iflag<nodata_opt.size();++iflag){
        vector<short>::iterator fit;
        fit=find(inputRange.begin(),inputRange.end(),static_cast<short>(nodata_opt[iflag]));
        if(fit!=inputRange.end())
          inputRange.erase(fit);
      }
      nclass=inputRange.size();
      if(verbose_opt[0]){
        cout << "nclass (inputRange.size()): " << nclass << endl;
        cout << "input range: " << endl;
      }
      if(classname_opt.size()){
        assert(classname_opt.size()==classvalue_opt.size());
        for(int iclass=0;iclass<classname_opt.size();++iclass){
          classValueMap[classname_opt[iclass]]=classvalue_opt[iclass];
          assert(classvalue_opt[iclass]<nameVector.size());
          nameVector[classvalue_opt[iclass]]=classname_opt[iclass];
        }
      }
      // nclass=classValueMap.size();
      for(int rc=0;rc<inputRange.size();++rc){
        classNames.push_back(type2string(inputRange[rc]));
        if(verbose_opt[0])
          cout << inputRange[rc] << endl;
      }
      cm.setClassNames(classNames);
      if(verbose_opt[0]){
        cout << "class names: " << endl;
        for(int iclass=0;iclass<cm.nClasses();++iclass)
          cout << iclass << " " << cm.getClass(iclass) << endl;
      }
      resultClass.resize(nclass,nclass);
      user.resize(nclass);
      producer.resize(nclass);
      nvalidation.resize(nclass);
      //initialize
      for(int rc=0;rc<nclass;++rc){
        for(int ic=0;ic<nclass;++ic)
          resultClass[rc][ic]=0;
        nvalidation[rc]=0;
      }
    }

    bool isDifferent=false;
    bool refIsRaster=false;

    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    float progress=0;
    // if(reference_opt[0].find(".shp")!=string::npos){
    ImgRaster gdalWriter;
    // this->open(input_opt[0],memory_opt[0]);
    if(mask_opt.size())
      maskReader.open(mask_opt[0]);
    // maskReader.open(mask_opt[0],memory_opt[0]);
    if(output_opt.size()){
      if(verbose_opt[0])
        cout << "opening output image " << output_opt[0] << endl;
      if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
        string theInterleave="INTERLEAVE=";
        theInterleave+=this->getInterleave();
        option_opt.push_back(theInterleave);
      }
      // gdalWriter.open(output_opt[0],this->nrOfCol(),this->nrOfRow(),1,this->getDataType(),oformat_opt[0],memory_opt[0],option_opt);
      gdalWriter.open(output_opt[0],this->nrOfCol(),this->nrOfRow(),1,this->getDataType(),oformat_opt[0],0,option_opt);
      if(nodata_opt.size())
        gdalWriter.GDALSetNoDataValue(nodata_opt[0]);
      gdalWriter.copyGeoTransform(shared_from_this());
      if(colorTable_opt.size())
        gdalWriter.setColorTable(colorTable_opt[0]);
      else if(this->getColorTable()!=NULL){
        if(verbose_opt[0])
          cout << "set colortable from input image" << endl;
        gdalWriter.setColorTable(this->getColorTable());
      }
    }

    //todo: support different data types!
    vector<double> lineInput(this->nrOfCol());
    vector<double> lineMask(maskReader.nrOfCol());
    vector<double> lineOutput;
    vector<double> bufferInput;//for regression
    vector<double> bufferReference;//for regression
    if(output_opt.size())
      lineOutput.resize(this->nrOfCol());

    unsigned int irow=0;
    unsigned int icol=0;
    double oldreferencerow=-1;
    double oldmaskrow=-1;

    // imgReference->open(reference_opt[0]);//,rmagicX_opt[0],rmagicY_opt[0]);
    if(this->isGeoRef()){
      assert(imgReference->isGeoRef());
      // if(this->getProjection()!=imgReference->getProjection())
      //    << "Warning: projection of input image and reference image are different" << endl;
    }
    vector<double> lineReference(imgReference->nrOfCol());
    if(confusion_opt[0]){
      imgReference->getRange(referenceRange,band_opt[1]);
      for(int iflag=0;iflag<nodata_opt.size();++iflag){
        vector<short>::iterator fit;
        fit=find(referenceRange.begin(),referenceRange.end(),static_cast<unsigned short>(nodata_opt[iflag]));
        if(fit!=referenceRange.end())
          referenceRange.erase(fit);
      }
      if(verbose_opt[0]){
        cout << "reference range: " << endl;
        for(int rc=0;rc<referenceRange.size();++rc)
          cout << referenceRange[rc] << endl;
      }
      if(referenceRange.size()!=inputRange.size()){
        if(confusion_opt[0]||output_opt.size()){
          cout << "reference range is not equal to input range!" << endl;
          cout << "Kappa: " << 0 << endl;
          cout << "total weighted: " << 0 << endl;
        }
        else{
          std::ostringstream errorStream;
          errorStream << "reference range is not equal to input range!" << endl;
          errorStream << "Images are different" << endl;
          throw(errorStream.str());
        }
      }
    }
    double rmse=0;
    // for(irow=0;irow<this->nrOfRow()&&!isDifferent;++irow){
    for(irow=0;irow<this->nrOfRow();++irow){
      //read line in lineInput, lineReference and lineMask
      this->readData(lineInput,irow,band_opt[0]);
      double x,y;//geo coordinates
      double ireference,jreference;//image coordinates in reference image
      double imask,jmask;//image coordinates in mask image
      for(icol=0;icol<this->nrOfCol();++icol){
        //find col in reference
        this->image2geo(icol,irow,x,y);
        imgReference->geo2image(x,y,ireference,jreference);
        if(ireference<0||ireference>=imgReference->nrOfCol()){
          if(rmse_opt[0]||regression_opt[0])
            continue;
          else{
            std::ostringstream errorStream;
             errorStream << ireference << " out of reference range!" << endl;
             errorStream << x << " " << y << " " << icol << " " << irow << endl;
             errorStream << x << " " << y << " " << ireference << " " << jreference << endl;
            throw(errorStream.str());
          }
        }
        if(jreference!=oldreferencerow){
          if(jreference<0||jreference>=imgReference->nrOfRow()){
            if(rmse_opt[0]||regression_opt[0])
              continue;
            else{
              std::ostringstream errorStream;
              errorStream << jreference << " out of reference range!" << endl;
              errorStream << x << " " << y << " " << icol << " " << irow << endl;
              errorStream << x << " " << y << " " << ireference << " " << jreference << endl;
              throw(errorStream.str());
            }
          }
          else{
            imgReference->readData(lineReference,static_cast<unsigned int>(jreference),band_opt[1]);
            oldreferencerow=jreference;
          }
        }
        bool flagged=false;
        for(int iflag=0;iflag<nodata_opt.size();++iflag){
          if((lineInput[icol]==nodata_opt[iflag])||(lineReference[ireference]==nodata_opt[iflag])){
            if(output_opt.size())
              lineOutput[icol]=nodata_opt[iflag];
            flagged=true;
            break;
          }
        }
        if(mask_opt.size()){
          maskReader.geo2image(x,y,imask,jmask);
          if(jmask>=0&&jmask<maskReader.nrOfRow()){
            if(jmask!=oldmaskrow)
              maskReader.readData(lineMask,jmask);
            for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
              if(lineMask[icol]==msknodata_opt[ivalue]){
                flagged=true;
                break;
              }
            }
          }
        }
        if(!flagged){
          if(rmse_opt[0]){//divide by image size to prevent overflow. At the end we need to take care about flagged pixels by normalizing...
            rmse+=static_cast<double>(lineInput[icol]-lineReference[ireference])*(lineInput[icol]-lineReference[ireference])/this->nrOfCol()/this->nrOfRow();
          }
          else if(regression_opt[0]){
            bufferInput.push_back(lineInput[icol]);
            bufferReference.push_back(lineReference[ireference]);
          }
          if(confusion_opt[0]){
            ++ntotalValidation;
            int rc=distance(referenceRange.begin(),find(referenceRange.begin(),referenceRange.end(),lineReference[ireference]));
            int ic=distance(inputRange.begin(),find(inputRange.begin(),inputRange.end(),lineInput[icol]));
            assert(rc<nclass);
            assert(ic<nclass);
            ++nvalidation[rc];
            ++resultClass[rc][ic];
            if(verbose_opt[0]>1)
              cout << "increment: " << rc << " " << referenceRange[rc] << " " << ic << " " << inputRange[ic] << endl;
            cm.incrementResult(cm.getClass(rc),cm.getClass(ic),1);
          }
          if(lineInput[icol]==lineReference[ireference]){//correct
            if(output_opt.size()){
              lineOutput[icol]=valueE_opt[0];
              if(nodata_opt.size()){
                if(valueE_opt[0]==nodata_opt[0])
                  lineOutput[icol]=lineInput[icol];
              }
            }
          }
          else{//error
            if(output_opt.empty()&&!confusion_opt[0]&&!rmse_opt[0]&&!regression_opt[0]){
              isDifferent=true;
              break;
            }
            if(output_opt.size()){
              if(lineInput[icol]>lineReference[ireference])
                lineOutput[icol]=valueO_opt[0];//omission error
              else
                lineOutput[icol]=valueC_opt[0];//commission error
            }
          }
        }
        else{
          ++nflagged;
          if(output_opt.size()){
            if(nodata_opt.size())
              lineOutput[icol]=nodata_opt[0];
            else //should never occur?
              lineOutput[icol]=0;
          }
        }
      }
      if(output_opt.size()){
        gdalWriter.writeData(lineOutput,irow);
      }
      else if(isDifferent&&!confusion_opt[0]&&!rmse_opt[0]&&!regression_opt[0]){//we can break off here, files are different...
        if(!verbose_opt[0])
          pfnProgress(1.0,pszMessage,pProgressArg);
        break;
      }
      progress=static_cast<float>(irow+1.0)/this->nrOfRow();
      if(!verbose_opt[0])
        pfnProgress(progress,pszMessage,pProgressArg);
    }
    if(output_opt.size())
      gdalWriter.close();
    else if(!confusion_opt[0]){
      if(rmse_opt[0]){
        double normalization=1.0*(this->nrOfCol())*(this->nrOfRow())/(this->nrOfCol()*(this->nrOfRow())-nflagged);
        if(verbose_opt[0]){
          cout << "normalization: " << normalization << endl;
          cout << "rmse before sqrt and normalization: " << rmse << endl;
        }
        cout << "--rmse " << sqrt(rmse/normalization) << endl;
      }
      else if(regression_opt[0]){
        double err=0;
        double c0=0;
        double c1=1;
        statfactory::StatFactory stat;
        if(bufferInput.size()&&bufferReference.size()){
          err=stat.linear_regression_err(bufferInput,bufferReference,c0,c1);
        }
        if(verbose_opt[0]){
          cout << "bufferInput.size(): " << bufferInput.size() << endl;
          cout << "bufferReference.size(): " << bufferReference.size() << endl;
          double theMin=0;
          double theMax=0;
          stat.minmax(bufferInput,bufferInput.begin(),bufferInput.end(),theMin,theMax);
          cout << "min,  max input: " << theMin << ", " << theMax << endl;
          theMin=0;
          theMax=0;
          stat.minmax(bufferReference,bufferReference.begin(),bufferReference.end(),theMin,theMax);
          cout << "min,  max reference: " << theMin << ", " << theMax << endl;
        }
        cout << "--c0 " << c0 << "--c1 " << c1 << " --rmse: " << err << endl;

      }
      else if(isDifferent){
        std::ostringstream outputStream;
        outputStream << "Error: no input file provided (use option -i). Use --help for help information" << std::endl;
        outputStream << "Images are different" << endl;
        throw(outputStream.str());
      }
      else{
        std::ostringstream outputStream;
        outputStream << "Images are identical" << endl;
        throw(outputStream.str());
      }
    }
    // imgReference->close();
    // this->close();
    if(mask_opt.size())
      maskReader.close();

    if(confusion_opt[0]){
      cm.setFormat(cmformat_opt[0]);
      cm.reportSE95(se95_opt[0]);
      ofstream outputFile;
      if(cmoutput_opt.size()){
        outputFile.open(cmoutput_opt[0].c_str(),ios::out);
        outputFile << cm << endl;
      }
      else
        cout << cm << endl;
    }
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
}

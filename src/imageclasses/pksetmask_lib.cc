/**********************************************************************
pksetmask_lib.cc: program to apply mask image (set invalid values) to raster image
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

#include "imageclasses/ImgRasterGdal.h"
#include "base/Optionpk.h"
#include "apps/AppFactory.h"

using namespace std;
using namespace app;

/**
 * @param app application specific option arguments
 * @return output image
 **/
shared_ptr<ImgRasterGdal> ImgRasterGdal::setMask(const AppFactory& app){
  shared_ptr<ImgRasterGdal> imgWriter=ImgRasterGdal::createImg();
  setMask(*imgWriter, app);
  return(imgWriter);
}

/**
 * @param imgWriter output raster setmask dataset
 * @return CE_None if successful, CE_Failure if failed
 **/
CPLErr ImgRasterGdal::setMask(ImgRasterGdal& imgWriter, const AppFactory& app){
  //command line options
  Optionpk<string> mask_opt("m", "mask", "Mask image(s)");
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<int> msknodata_opt("msknodata", "msknodata", "Mask value(s) where image has nodata. Use one value for each mask, or multiple values for a single mask.", 1);
  Optionpk<short> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed). Provide band for each mask.", 0);
  Optionpk<char> operator_opt("p", "operator", "Operator: < = > !. Use operator for each msknodata option", '=');
  Optionpk<int> nodata_opt("nodata", "nodata", "nodata value to put in image if not valid", 0);
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0,2);

  otype_opt.setHide(1);
  mskband_opt.setHide(1);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=mask_opt.retrieveOption(app.getArgc(),app.getArgv());
    msknodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    mskband_opt.retrieveOption(app.getArgc(),app.getArgv());
    nodata_opt.retrieveOption(app.getArgc(),app.getArgv());
    operator_opt.retrieveOption(app.getArgc(),app.getArgv());
    otype_opt.retrieveOption(app.getArgc(),app.getArgv());
    verbose_opt.retrieveOption(app.getArgc(),app.getArgv());
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    return(CE_Failure);
  }
  if(!doProcess){
    cout << endl;
    std::ostringstream helpStream;
    helpStream << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    throw(helpStream.str());//help was invoked, stop processing
  }

  if(verbose_opt[0])
    cout << "number of mask images: " << mask_opt.size() << endl;

  //duplicate band used for mask if not explicitly provided
  while(mskband_opt.size()<mask_opt.size())
    mskband_opt.push_back(mskband_opt[0]);

  vector<ImgRasterGdal> maskReader(mask_opt.size());
  for(int imask=0;imask<mask_opt.size();++imask){
    if(verbose_opt[0])
      cout << "opening mask image file " << mask_opt[imask] << endl;
    maskReader[imask].open(mask_opt[imask]);
  }

  GDALDataType theType=GDT_Unknown;
  if(otype_opt.size()){
    theType=getGDALDataType(otype_opt[0]);
    if(theType==GDT_Unknown)
      std::cout << "Warning: unknown output pixel type: " << otype_opt[0] << ", using input type as default" << std::endl;
  }
  //if output type not set, get type from input image
  if(theType==GDT_Unknown){
    theType=getDataType();
    if(verbose_opt[0])
      cout << "Using data type from input image: " << GDALGetDataTypeName(theType) << endl;
  }
  if(verbose_opt[0])
    cout << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
  // ImgRasterGdal imgWriter;
  try{
    imgWriter.open(nrOfCol(),nrOfRow(),nrOfBand(),theType);
    for(unsigned int iband=0;iband<nrOfBand();++iband)
      imgWriter.GDALSetNoDataValue(nodata_opt[0],iband);
    imgWriter.setProjection(getProjection());
    imgWriter.copyGeoTransform(*this);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    return(CE_Failure);
  }
  // if(verbose_opt[0])
  //   cout << "opening output image file " << output_opt[0] << endl;
  // imgWriter.open(output_opt[0],;
  if (getColorTable()!=NULL)//copy colorTable from input image
    imgWriter.setColorTable(getColorTable());
  while(nodata_opt.size()<msknodata_opt.size())
    nodata_opt.push_back(nodata_opt.back());
  if(operator_opt.size()!=msknodata_opt.size()&&operator_opt.size()!=1){
    std::string errorString="Error: number of operators and masks do not match";
    throw(errorString);
  }
  if(verbose_opt[0]){
    cout << " mask files selected: " << mask_opt.size() << endl;
    for(int iv=0;iv<msknodata_opt.size();++iv){
      char op=(operator_opt.size()==msknodata_opt.size())?operator_opt[iv]:operator_opt[0];
      cout << op << " " << msknodata_opt[iv] << "->" << nodata_opt[iv] << endl;
    }
  }

  Vector2d<double> lineInput(nrOfBand(),nrOfCol());
  Vector2d<double> lineOutput(imgWriter.nrOfBand(),imgWriter.nrOfCol());
  assert(lineOutput.size()==lineInput.size());
  assert(nrOfCol()==imgWriter.nrOfCol());
  Vector2d<double> lineMask(mask_opt.size());
  for(int imask=0;imask<mask_opt.size();++imask){
    if(verbose_opt[0])
      cout << "mask " << imask << " has " << maskReader[imask].nrOfCol() << " columns and " << maskReader[imask].nrOfRow() << " rows" << endl;
    lineMask[imask].resize(maskReader[imask].nrOfCol());
  }
  unsigned int irow=0;
  unsigned int icol=0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  float progress=0;
  if(!verbose_opt[0])
    pfnProgress(progress,pszMessage,pProgressArg);
  vector<double> oldRowMask(mask_opt.size());
  for(int imask=0;imask<mask_opt.size();++imask)
    oldRowMask[imask]=-1;
  for(irow=0;irow<nrOfRow();++irow){
    //read line in lineInput buffer
    for(unsigned int iband=0;iband<nrOfBand();++iband){
      try{
        readData(lineInput[iband],irow,iband);
      }
      catch(string errorstring){
        cerr << errorstring << endl;
        exit(1);
      }
    }
    double x,y;//geo coordinates
    double colMask,rowMask;//image coordinates in mask image
    for(icol=0;icol<nrOfCol();++icol){
      if(mask_opt.size()>1){//multiple masks
        for(int imask=0;imask<mask_opt.size();++imask){
          image2geo(icol,irow,x,y);
          maskReader[imask].geo2image(x,y,colMask,rowMask);
          colMask=static_cast<unsigned int>(colMask);
          rowMask=static_cast<unsigned int>(rowMask);
          bool masked=false;
          if(rowMask>=0&&rowMask<maskReader[imask].nrOfRow()&&colMask>=0&&colMask<maskReader[imask].nrOfCol()){
            if(static_cast<unsigned int>(rowMask)!=static_cast<unsigned int>(oldRowMask[imask])){
              assert(rowMask>=0&&rowMask<maskReader[imask].nrOfRow());
              try{
                maskReader[imask].readData(lineMask[imask],static_cast<unsigned int>(rowMask),mskband_opt[imask]);
              }
              catch(string errorstring){
                cerr << errorstring << endl;
                exit(1);
              }
              oldRowMask[imask]=rowMask;
            }
          }
          else
            continue;//no coverage in this mask
          int ivalue=0;
          if(mask_opt.size()==msknodata_opt.size())//one invalid value for each mask
            ivalue=msknodata_opt[imask];
          else//use same invalid value for each mask
            ivalue=msknodata_opt[0];
          char op=(operator_opt.size()==mask_opt.size())?operator_opt[imask]:operator_opt[0];
          switch(op){
          case('='):
          default:
            if(lineMask[imask][colMask]==ivalue)
              masked=true;
          break;
          case('<'):
            if(lineMask[imask][colMask]<ivalue)
              masked=true;
            break;
          case('>'):
            if(lineMask[imask][colMask]>ivalue)
              masked=true;
            break;
          case('!'):
            if(lineMask[imask][colMask]!=ivalue)
              masked=true;
            break;
          }
          if(masked){
            if(verbose_opt[0]>1)
              cout << "image masked at (col=" << icol << ",row=" << irow <<") with mask " << mask_opt[imask] << " and value " << ivalue << endl;
            for(unsigned int iband=0;iband<nrOfBand();++iband){
              if(mask_opt.size()==nodata_opt.size())//one flag value for each mask
                lineInput[iband][icol]=nodata_opt[imask];
              else
                lineInput[iband][icol]=nodata_opt[0];
            }
            masked=false;
            break;
          }
        }
      }
      else{//potentially more invalid values for single mask
        image2geo(icol,irow,x,y);
        maskReader[0].geo2image(x,y,colMask,rowMask);
        colMask=static_cast<unsigned int>(colMask);
        rowMask=static_cast<unsigned int>(rowMask);
        bool masked=false;
        if(rowMask>=0&&rowMask<maskReader[0].nrOfRow()&&colMask>=0&&colMask<maskReader[0].nrOfCol()){
          if(static_cast<unsigned int>(rowMask)!=static_cast<unsigned int>(oldRowMask[0])){
            assert(rowMask>=0&&rowMask<maskReader[0].nrOfRow());
            try{
              // maskReader[0].readData(lineMask[0],static_cast<unsigned int>(rowMask));
              maskReader[0].readData(lineMask[0],static_cast<unsigned int>(rowMask),mskband_opt[0]);
            }
            catch(string errorstring){
              cerr << errorstring << endl;
              exit(1);
            }
            oldRowMask[0]=rowMask;
          }
          for(int ivalue=0;ivalue<msknodata_opt.size();++ivalue){
            assert(msknodata_opt.size()==nodata_opt.size());
            char op=(operator_opt.size()==msknodata_opt.size())?operator_opt[ivalue]:operator_opt[0];
            switch(op){
            case('='):
            default:
              if(lineMask[0][colMask]==msknodata_opt[ivalue])
                masked=true;
            break;
            case('<'):
              if(lineMask[0][colMask]<msknodata_opt[ivalue])
                masked=true;
              break;
            case('>'):
              if(lineMask[0][colMask]>msknodata_opt[ivalue])
                masked=true;
              break;
            case('!'):
              if(lineMask[0][colMask]!=msknodata_opt[ivalue])
                masked=true;
              break;
            }
            if(masked){
              for(unsigned int iband=0;iband<nrOfBand();++iband)
                lineInput[iband][icol]=nodata_opt[ivalue];
              masked=false;
              break;
            }
          }
        }
      }
      for(unsigned int iband=0;iband<lineOutput.size();++iband)
        lineOutput[iband][icol]=lineInput[iband][icol];
    }
    //write buffer lineOutput to output file
    for(unsigned int iband=0;iband<imgWriter.nrOfBand();++iband){
      try{
        imgWriter.writeData(lineOutput[iband],irow,iband);
      }
      catch(string errorstring){
        cerr << errorstring << endl;
        exit(1);
      }
    }
    //progress bar
    progress=static_cast<float>(irow+1.0)/imgWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  for(int imask=0;imask<mask_opt.size();++imask)
    maskReader[imask].close();
  return(CE_None);
}

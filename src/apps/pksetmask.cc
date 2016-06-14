/**********************************************************************
pksetmask.cc: program to apply mask image (set invalid values) to raster image
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

#include "imageclasses/ImgRaster.h"
#include "base/Optionpk.h"
/******************************************************************************/
/*! \page pksetmask pksetmask
 program to apply mask image (set invalid values) to raster image
## SYNOPSIS

<code>
  Usage: pksetmask -i input -m mask [-m mask]* -o output
</code>

<code>
  
  Options: [-mskband value]* [-msknodata value -nodata value]*
 
  Advanced options: [--operator '<'|'='|'<'] [-ot type] [-of format] [-co option]* [-ct table] 

</code>

\section pksetmask_description Description

The utility pksetmask sets a mask provided with option -m to an input raster dataset. The default operator is '='. Values in the input raster data where the mask has a nodata value (set with the option -msknodata) will then be set to nodata (set with -nodata). Other operators are less than (--operator '<') and larger than (--operator '>').\section pksetmask_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image | 
 | m      | mask                 | std::string |       |Mask image(s) | 
 | msknodata | msknodata            | int  | 1     |Mask value(s) where image has nodata. Use one value for each mask, or multiple values for a single mask. | 
 | mskband | mskband              | short | 0     |Mask band to read (0 indexed). Provide band for each mask. | 
 | o      | output               | std::string |       |Output mask file | 
 | nodata | nodata               | int  | 0     |nodata value to put in image if not valid | 
 | p      | operator             | char | =     |Operator: < = > !. Use operator for each msknodata option | 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | of     | oformat              | std::string | GTiff |Output image format (see also gdal_translate)| 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 
 | mem    | mem                  | unsigned long int | 0 |Buffer size (in MB) to read image data blocks in memory | 

Usage: pksetmask -i input -m mask [-m mask]* [-msknodata value -nodata value]* -o output

Examples
========
Some examples how to use pksetmask can be found \ref examples_pksetmask "here"
FAQ
========
Frequently asked questions on pksetmask can be found \ref faq_pksetmask "here"
**/

using namespace std;

int main(int argc, char *argv[])
{
  //command line options
  Optionpk<string> input_opt("i", "input", "Input image");
  Optionpk<string> mask_opt("m", "mask", "Mask image(s)");
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<string> otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string> oformat_opt("of", "oformat", "Output image format (see also gdal_translate)","GTiff");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<int> msknodata_opt("msknodata", "msknodata", "Mask value(s) where image has nodata. Use one value for each mask, or multiple values for a single mask.", 1);
  Optionpk<short> mskband_opt("mskband", "mskband", "Mask band to read (0 indexed). Provide band for each mask.", 0);
  Optionpk<char> operator_opt("p", "operator", "Operator: < = > !. Use operator for each msknodata option", '=');
  Optionpk<int> nodata_opt("nodata", "nodata", "nodata value to put in image if not valid", 0);
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<unsigned long int>  memory_opt("mem", "mem", "Buffer size (in MB) to read image data blocks in memory",0,1);
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0,2);

  otype_opt.setHide(1);
  oformat_opt.setHide(1);
  option_opt.setHide(1);
  colorTable_opt.setHide(1);
  mskband_opt.setHide(1);
  memory_opt.setHide(1);
  
  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    msknodata_opt.retrieveOption(argc,argv);
    mskband_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    operator_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    oformat_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    memory_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pksetmask -i input -m mask [-m mask]* [-msknodata value -nodata value]* -o output" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  if(verbose_opt[0])
     cout << "number of mask images: " << mask_opt.size() << endl;

  //duplicate band used for mask if not explicitly provided
  while(mskband_opt.size()<mask_opt.size())
    mskband_opt.push_back(mskband_opt[0]);

  vector<ImgRaster> maskReader(mask_opt.size()); 
  for(int imask=0;imask<mask_opt.size();++imask){
    if(verbose_opt[0])
      cout << "opening mask image file " << mask_opt[imask] << endl;
    maskReader[imask].open(mask_opt[imask],memory_opt[0]);
  }
  assert(input_opt.size());
  if(verbose_opt[0])
    cout << "opening input image file " << input_opt[0] << endl;
  ImgRaster inputReader;
  inputReader.open(input_opt[0],memory_opt[0]);
  string imageType;//=inputReader.getImageType();
  if(oformat_opt.size())//default
    imageType=oformat_opt[0];
  GDALDataType theType=GDT_Unknown;
  if(verbose_opt[0]){
    std::cout << "Image type: " << imageType << std::endl;
    std::cout << "possible output data types: ";
  }
  for(int iType = 0; iType < GDT_TypeCount; ++iType){
    if(verbose_opt[0])
      cout << " " << GDALGetDataTypeName((GDALDataType)iType);
    if( GDALGetDataTypeName((GDALDataType)iType) != NULL
        && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                 otype_opt[0].c_str()))
      theType=(GDALDataType) iType;
  }
  if(theType==GDT_Unknown)
    theType=inputReader.getDataType();

  assert(output_opt.size());
  if(verbose_opt[0]){
    std::cout << std::endl << "Output data type:  " << GDALGetDataTypeName(theType) << std::endl;
    std::cout << "opening output image for writing: " << output_opt[0] << std::endl;
  }
  ImgRaster outputWriter;
  try{
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=inputReader.getInterleave();
      option_opt.push_back(theInterleave);
    }
    outputWriter.open(output_opt[0],inputReader.nrOfCol(),inputReader.nrOfRow(),inputReader.nrOfBand(),theType,imageType,memory_opt[0],option_opt);
    for(unsigned int iband=0;iband<inputReader.nrOfBand();++iband)
      outputWriter.GDALSetNoDataValue(nodata_opt[0],iband);
    outputWriter.setProjection(inputReader.getProjection());
    outputWriter.copyGeoTransform(inputReader);
  }
  catch(string errorstring){
    cout << errorstring << endl;
    exit(1);
  }
  // if(verbose_opt[0])
  //   cout << "opening output image file " << output_opt[0] << endl;
  // outputWriter.open(output_opt[0],inputReader);
  if(colorTable_opt.size()){
    if(colorTable_opt[0]!="none")
      outputWriter.setColorTable(colorTable_opt[0]);
  }
  else if (inputReader.getColorTable()!=NULL)//copy colorTable from input image
    outputWriter.setColorTable(inputReader.getColorTable());
  if(inputReader.isGeoRef()){
    for(int imask=0;imask<mask_opt.size();++imask)
      assert(maskReader[imask].isGeoRef());
  }
  assert(nodata_opt.size()==msknodata_opt.size());
  assert(operator_opt.size()==msknodata_opt.size()||operator_opt.size()==1);
  if(verbose_opt[0]){
    cout << " mask files selected: " << mask_opt.size() << endl;
    for(int iv=0;iv<msknodata_opt.size();++iv){
      char op=(operator_opt.size()==msknodata_opt.size())?operator_opt[iv]:operator_opt[0];
      cout << op << " " << msknodata_opt[iv] << "->" << nodata_opt[iv] << endl;
    }
  }
  
  Vector2d<double> lineInput(inputReader.nrOfBand(),inputReader.nrOfCol());
  Vector2d<double> lineOutput(outputWriter.nrOfBand(),outputWriter.nrOfCol());
  assert(lineOutput.size()==lineInput.size());
  assert(inputReader.nrOfCol()==outputWriter.nrOfCol());
  // Vector2d<int> lineMask(mask_opt.size());
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
  // double oldRowMask=-1;
  vector<double> oldRowMask(mask_opt.size());
  for(int imask=0;imask<mask_opt.size();++imask)
    oldRowMask[imask]=-1;
  for(irow=0;irow<inputReader.nrOfRow();++irow){
    //read line in lineInput buffer
    for(unsigned int iband=0;iband<inputReader.nrOfBand();++iband){
      try{
        inputReader.readData(lineInput[iband],irow,iband);
      }
      catch(string errorstring){
        cerr << errorstring << endl;
        exit(1);
      }
    }
    double x,y;//geo coordinates
    double colMask,rowMask;//image coordinates in mask image
    for(icol=0;icol<inputReader.nrOfCol();++icol){
      if(mask_opt.size()>1){//multiple masks
        for(int imask=0;imask<mask_opt.size();++imask){
	  inputReader.image2geo(icol,irow,x,y);
	  maskReader[imask].geo2image(x,y,colMask,rowMask);
	  colMask=static_cast<unsigned int>(colMask);
	  rowMask=static_cast<unsigned int>(rowMask);
          bool masked=false;
          if(rowMask>=0&&rowMask<maskReader[imask].nrOfRow()&&colMask>=0&&colMask<maskReader[imask].nrOfCol()){
	    if(static_cast<unsigned int>(rowMask)!=static_cast<unsigned int>(oldRowMask[imask])){
	      assert(rowMask>=0&&rowMask<maskReader[imask].nrOfRow());
	      try{
		// maskReader[imask].readData(lineMask[imask],static_cast<unsigned int>(rowMask));
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
	    for(unsigned int iband=0;iband<inputReader.nrOfBand();++iband){
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
	inputReader.image2geo(icol,irow,x,y);
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
              for(unsigned int iband=0;iband<inputReader.nrOfBand();++iband)
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
    for(unsigned int iband=0;iband<outputWriter.nrOfBand();++iband){
      try{
        outputWriter.writeData(lineOutput[iband],irow,iband);
      }
      catch(string errorstring){
        cerr << errorstring << endl;
        exit(1);
      }
    }
    //progress bar
    progress=static_cast<float>(irow+1.0)/outputWriter.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
  inputReader.close();
  for(int imask=0;imask<mask_opt.size();++imask)
    maskReader[imask].close();
  outputWriter.close();
}

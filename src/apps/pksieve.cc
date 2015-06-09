/**********************************************************************
pksieve.cc: program to sieve filter raster image
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
#include "cpl_string.h"
#include "gdal_priv.h"
#include "gdal.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
// #include "imageclasses/ImgWriterOgr.h"
#include "base/Optionpk.h"
#include "ogrsf_frmts.h"
extern "C" {
#include "gdal_alg.h"
#include "ogr_api.h"
}
/******************************************************************************/
/*! \page pksieve pksieve
 program to sieve filter raster image
## SYNOPSIS

<code>
  Usage: pksieve -i input [-s size] -o output
</code>

<code>
  
  Options: [-c 4|8] [-b band] [-m mask] [-ot type] [-of format] [-co option]* [-ct table] 

</code>

\section pksieve_description Description

The utility pksieve filters small objects (maximum size defined with the option -s) in a raster by replacing them to the largest neighbor object. In this context, objects are defined as pixels of the same value that are also connected. The connection can be defined in four directions (N-S and W-E: set option -c 4) or eight directions (N-S, W-E and diagonals NW-SE, NE-SW: set option -c 8).\section pksieve_options Options
 - use either `-short` or `--long` options (both `--long=value` and `--long value` are supported)
 - short option `-h` shows basic options only, long option `--help` shows all options
|short|long|type|default|description|
|-----|----|----|-------|-----------|
 | i      | input                | std::string |       |Input image file | 
 | s      | size                 | int  | 0     |raster polygons with sizes smaller than this will be merged into their largest neighbour. No sieve is performed if size = 0 | 
 | o      | output               | std::string |       |Output image file | 
 | c      | connect              | int  | 8     |the connectedness: 4 directions or 8 directions | 
 | b      | band                 | int  | 0     |the band to be used from input file | 
 | m      | mask                 | std::string |       |Use the first band of the specified file as a validity mask (zero is invalid, non-zero is valid). | 
 | ot     | otype                | std::string |       |Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image | 
 | co     | co                   | std::string |       |Creation option for output file. Multiple options can be specified. | 
 | ct     | ct                   | std::string |       |color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid) | 

Usage: pksieve -i input [-s size] -o output


Examples
========
Some examples how to use pksieve can be found \ref examples_pksieve "here"
**/

using namespace std;

int main(int argc,char **argv) {
  Optionpk<string> input_opt("i", "input", "Input image file");
  Optionpk<string> mask_opt("m", "mask", "Use the first band of the specified file as a validity mask (zero is invalid, non-zero is valid).");
  Optionpk<string> output_opt("o", "output", "Output image file");
  Optionpk<int> band_opt("b", "band", "the band to be used from input file", 0);
  Optionpk<int> connect_opt("c", "connect", "the connectedness: 4 directions or 8 directions", 8);
  Optionpk<int> size_opt("s", "size", "raster polygons with sizes smaller than this will be merged into their largest neighbour. No sieve is performed if size = 0", 0);
  Optionpk<string>  otype_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode if > 0", 0,2);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    size_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    connect_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    otype_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pksieve -i input [-s size] -o output" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  GDALAllRegister();

  double dfComplete=0.0;
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  pfnProgress(dfComplete,pszMessage,pProgressArg);
  
  ImgReaderGdal maskReader;
  GDALRasterBand *maskBand=NULL;
  if(mask_opt.size()){
    if(verbose_opt[0])
      cout << "opening mask file " << mask_opt[0] << endl;
    maskReader.open(mask_opt[0]);
    maskBand = maskReader.getRasterBand(0);
  }

  assert(input_opt.size());
  assert(output_opt.size());
  ImgReaderGdal inputReader(input_opt[0]);
  GDALRasterBand  *inputBand;
  inputBand=inputReader.getRasterBand(band_opt[0]);

  ImgWriterGdal outputWriter;
  GDALRasterBand *outputBand=NULL;
  if(verbose_opt[0])
    cout << "opening output file " << output_opt[0] << endl;
  outputWriter.open(output_opt[0],inputReader);
  if(colorTable_opt.size()){
    if(colorTable_opt[0]!="none")
      outputWriter.setColorTable(colorTable_opt[0]);
  }
  else if (inputReader.getColorTable()!=NULL)//copy colorTable from input image
    outputWriter.setColorTable(inputReader.getColorTable());
  outputBand = outputWriter.getRasterBand(0);
  //sieve filter to remove small raster elements (overwrite input band)
  if(size_opt[0]){
    if(GDALSieveFilter((GDALRasterBandH)inputBand, (GDALRasterBandH)maskBand, (GDALRasterBandH)outputBand, size_opt[0], connect_opt[0],NULL,pfnProgress,pProgressArg)!=CE_None)
      cerr << CPLGetLastErrorMsg() << endl;
    else{
      dfComplete=1.0;
      pfnProgress(dfComplete,pszMessage,pProgressArg);
    }
  }
  inputReader.close();
  if(mask_opt.size())
    maskReader.close();
  outputWriter.close();
}


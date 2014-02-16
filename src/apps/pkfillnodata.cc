/**********************************************************************
pkfillnodata.cc: program to fill holes in raster image
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
#include "cpl_string.h"
#include "gdal_priv.h"
#include "gdal.h"
extern "C" {
#include "gdal_alg.h"
}
#include <string>
#include "base/Optionpk.h"

int main(int argc,char **argv) {
  Optionpk<std::string> input_opt("i", "input", "Input image file");
  Optionpk<int> band_opt("b", "band", "band(s) to process (Default is -1: process all bands)");
  Optionpk<std::string> mask_opt("m", "mask", "Mask band indicating pixels to be interpolated (zero valued) ");
  Optionpk<std::string> output_opt("o", "output", "Output image file");
  Optionpk<double> distance_opt("d", "distance", "Maximum number of pixels to search in all directions to find values to interpolate from", 0);
  Optionpk<int> iteration_opt("it", "iteration", "Number of 3x3 smoothing filter passes to run (default 0)", 0);
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    distance_opt.retrieveOption(argc,argv);
    iteration_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(std::string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  assert(input_opt.size());
  assert(mask_opt.size());
  assert(output_opt.size());
  GDALAllRegister();
  GDALDataset *gds_input;
  if(verbose_opt[0])
    std::cout << "opening input file " << input_opt[0] << std::endl;
  gds_input = (GDALDataset *) GDALOpen(input_opt[0].c_str(), GA_ReadOnly);
  if(gds_input == NULL){
    std::string errorString="FileOpenError";
    throw(errorString);
  }

  GDALDataset *gds_mask;
  if(verbose_opt[0])
    std::cout << "opening mask file " << mask_opt[0] << std::endl;
  gds_mask = (GDALDataset *) GDALOpen(mask_opt[0].c_str(), GA_ReadOnly );
  if(gds_mask == NULL){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
  GDALRasterBand *maskBand;
  if(verbose_opt[0])
    std::cout << "get mask raster band" << std::endl;
  maskBand=gds_mask->GetRasterBand(1);
  

  GDALDriver *poDriver;
  poDriver = GetGDALDriverManager()->GetDriverByName(gds_input->GetDriver()->GetDescription());
  if( poDriver == NULL ){
    std::string errorString="FileOpenError";
    throw(errorString);
  }
  if(verbose_opt[0])
    std::cout << "copying input file to " << output_opt[0] << std::endl;
  poDriver->CopyFiles(output_opt[0].c_str(),input_opt[0].c_str());
  GDALDataset *gds_out;
  gds_out=(GDALDataset *) GDALOpen(output_opt[0].c_str(), GA_Update);

  if(band_opt.empty()){
    band_opt.clear();
    for(int iband=0;iband<gds_input->GetRasterCount();++iband)
      band_opt.push_back(iband);
  }
  GDALRasterBand *targetBand;
  for(unsigned short iband=0;iband<band_opt.size();++iband){
    targetBand=gds_out->GetRasterBand(band_opt[iband]+1);
    if(verbose_opt[0])
      std::cout << "copying input file to " << output_opt[0] << std::endl;
    double dfComplete=0.0;
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    pfnProgress(dfComplete,pszMessage,pProgressArg);
    if(GDALFillNodata(targetBand,maskBand,distance_opt[0],0,iteration_opt[0],NULL,pfnProgress,pProgressArg)!=CE_None){
      std::cerr << CPLGetLastErrorMsg() << std::endl;
      exit(1);
    }
    else{
      dfComplete=1.0;
      pfnProgress(dfComplete,pszMessage,pProgressArg);
    }

    //   gds_out=poDriver->CreateCopy(output_opt[0].c_str(),gds_input, FALSE,NULL,NULL,NULL);
    //   char **papszParmList=NULL;
    //   gds_out=poDriver->Create(output_opt[0].c_str(),targetBand->GetXSize(),targetBand->GetYSize(),1,targetBand->GetRasterDataType(),papszParmList);
    //   char **papszMetadata;
    //   papszMetadata = poDriver->GetMetadata();
    //   assert( CSLFetchBoolean( papszMetadata, GDAL_DCAP_CREATE, FALSE ));
    //   ostringstream compressList;
    //   gds_out->SetMetadataItem("INTERLEAVE",gds_input->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE"),"IMAGE_STRUCTURE");
    //   gds_out->SetMetadataItem("COMPRESSION",gds_input->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE"),"IMAGE_STRUCTURE");
    //   if(gds_input->GetProjectionRef()!=NULL){
    //     gds_out->SetProjection(gds_input->GetProjectionRef());
    //     double adfGeoTransform[6];
    //     gds_input->GetGeoTransform(adfGeoTransform);
    //     gds_out->SetGeoTransform(adfGeoTransform);
    //   }
  }
  GDALClose(gds_input);
  GDALClose(gds_mask);
  GDALClose(gds_out);
  GDALDumpOpenDatasets(stderr);
  GDALDestroyDriverManager();
}


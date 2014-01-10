/**********************************************************************
pkreclass.cc: program to replace pixel values in raster image
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
#include <map>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "imageclasses/ImgWriterOgr.h"
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"

using namespace std;

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input image");
  Optionpk<string> mask_opt("m", "mask", "Mask image(s)");
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<unsigned short> masknodata_opt("msknodata", "msknodata", "Mask value(s) where image has nodata. Use one value for each mask, or multiple values for a single mask.", 1);
  Optionpk<int> nodata_opt("nodata", "nodata", "nodata value to put in image if not valid (0)", 0);
  Optionpk<string> colorTable_opt("ct", "ct", "color table (file with 5 columns: id R G B ALFA (0: transparent, 255: solid)");
  Optionpk<unsigned short>  band_opt("b", "band", "band index to replace (other bands are copied to output)", 0);
  Optionpk<string> type_opt("ot", "otype", "Data type for output image ({Byte/Int16/UInt16/UInt32/Int32/Float32/Float64/CInt16/CInt32/CFloat32/CFloat64}). Empty string: inherit type from input image", "");
  Optionpk<string> code_opt("code", "code", "Recode text file (2 colums: from to)");
  Optionpk<string> class_opt("c", "class", "list of classes to reclass (in combination with reclass option)");
  Optionpk<string> reclass_opt("r", "reclass", "list of recoded class(es) (in combination with class option)");
  Optionpk<string> fieldname_opt("n", "fname", "field name of the shape file to be replaced", "label");
  Optionpk<string> option_opt("co", "co", "Creation option for output file. Multiple options can be specified.");
  Optionpk<string> description_opt("d", "description", "Set image description");
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    mask_opt.retrieveOption(argc,argv);
    masknodata_opt.retrieveOption(argc,argv);
    nodata_opt.retrieveOption(argc,argv);
    code_opt.retrieveOption(argc,argv);
    class_opt.retrieveOption(argc,argv);
    reclass_opt.retrieveOption(argc,argv);
    colorTable_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    type_opt.retrieveOption(argc,argv);
    band_opt.retrieveOption(argc,argv);
    fieldname_opt.retrieveOption(argc,argv);
    option_opt.retrieveOption(argc,argv);
    description_opt.retrieveOption(argc,argv);
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
  if(input_opt.empty()){
    std::cerr << "No input file provided (use option -i). Use --help for help information" << std::endl;
    exit(0);
  }
  if(output_opt.empty()){
    std::cerr << "No output file provided (use option -o). Use --help for help information" << std::endl;
    exit(0);
  }
    
  // vector<short> bandVector;
  // for(int iband=0;iband<band_opt.size();++iband)
  //   bandVector.push_back(band_opt[iband]);
  map<string,string> codemapString;//map with codes: codemapString[theKey(from)]=theValue(to)
  map<double,double> codemap;//map with codes: codemap[theKey(from)]=theValue(to)
  if(code_opt.size()){
    if(verbose_opt[0])
      cout << "opening code text file " << code_opt[0] << endl;
    ifstream codefile;
    codefile.open(code_opt[0].c_str());
    string theKey;
    string theValue;
    while(codefile>>theKey){
      codefile >> theValue;
      codemapString[theKey]=theValue;
      codemap[string2type<double>(theKey)]=string2type<double>(theValue);
    }
    codefile.close();
  }
  else{//use combination of class_opt and reclass_opt
    assert(class_opt.size()==reclass_opt.size());
    for(int iclass=0;iclass<class_opt.size();++iclass){
      codemapString[class_opt[iclass]]=reclass_opt[iclass];
      codemap[string2type<double>(class_opt[iclass])]=string2type<double>(reclass_opt[iclass]);
    }
  }
  assert(codemapString.size());
  assert(codemap.size());
  //if verbose true, print the codes to screen
  if(verbose_opt[0]){
    map<string,string>::iterator mit;
    cout << codemapString.size() << " codes used: " << endl;
    for(mit=codemapString.begin();mit!=codemapString.end();++mit)
      cout << (*mit).first << " " << (*mit).second << endl;
  }
  bool refIsRaster=false;
  ImgReaderOgr ogrReader;
  try{
    ogrReader.open(input_opt[0]);
  }
  catch(string errorString){
    refIsRaster=true;
  }
  // if(input_opt[0].find(".shp")!=string::npos){//shape file
  if(!refIsRaster){
    if(verbose_opt[0])
      cout << "opening " << input_opt[0] << " for reading " << endl;
    // ImgReaderOgr ogrReader(input_opt[0]);
    if(verbose_opt[0])
      cout << "opening " << output_opt[0] << " for writing " << endl;
    ImgWriterOgr ogrWriter(output_opt[0],ogrReader);
    if(verbose_opt[0])
      cout << "copied layer from " << input_opt[0] << endl << flush;
    OGRFeatureDefn *poFDefn = ogrWriter.getLayer()->GetLayerDefn();
    //start reading features from the layer
    if(verbose_opt[0])
      cout << "reset reading" << endl;
    ogrReader.getLayer()->ResetReading();
    unsigned long int ifeature=0;
    if(verbose_opt[0])
      cout << "going through features" << endl << flush;
    while(true){
//     while( (poFeature = ogrWriter.getLayer()->GetNextFeature()) != NULL ){
      OGRFeature *poFeature;
      poFeature=ogrReader.getLayer()->GetNextFeature();
      if(poFeature== NULL)
        break;
      OGRFeatureDefn *poFDefn = ogrWriter.getLayer()->GetLayerDefn();
      string featurename;
      for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
        OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
        string fieldname=poFieldDefn->GetNameRef();
        if(fieldname==fieldname_opt[0]){
          string fromClass=poFeature->GetFieldAsString(iField);
          string toClass=fromClass;
          if(codemapString.find(fromClass)!=codemapString.end())
            toClass=codemapString[fromClass];
          poFeature->SetField(iField,toClass.c_str());
          if(verbose_opt[0])
            cout << "feature " << ifeature << ": " << fromClass << "->" << poFeature->GetFieldAsInteger(iField) << endl << flush;
//             cout << "feature " << ifeature << ": " << fromClass << "->" << toClass << endl << flush;
        }
      }
      //do not forget to actually write feature to file!!!
      ogrWriter.createFeature(poFeature);
      OGRFeature::DestroyFeature( poFeature );
      ++ifeature;
    }
    if(verbose_opt[0])
       cout << "replaced " << ifeature << " features" << endl;
    ogrReader.close();
    ogrWriter.close();
  }
  else{//image file
    ImgReaderGdal inputReader;
    vector<ImgReaderGdal> maskReader(mask_opt.size()); 
    ImgWriterGdal outputWriter;
    if(verbose_opt[0])
      cout << "opening input image file " << input_opt[0] << endl;
    inputReader.open(input_opt[0]);
    for(int imask=0;imask<mask_opt.size();++imask){
      if(verbose_opt[0])
        cout << "opening mask image file " << mask_opt[imask] << endl;
      maskReader[imask].open(mask_opt[imask]);
    }
    if(verbose_opt[0]){
      cout << "opening output image file " << output_opt[0] << endl;
      cout << "data type: " << type_opt[0] << endl;
    }
    //create output image with user defined data type 
    GDALDataType theType=GDT_Unknown;
    if(verbose_opt[0])
      cout << "possible output data types: ";
    for(int iType = 0; iType < GDT_TypeCount; ++iType){
      if(verbose_opt[0])
        cout << " " << GDALGetDataTypeName((GDALDataType)iType);
      if( GDALGetDataTypeName((GDALDataType)iType) != NULL
          && EQUAL(GDALGetDataTypeName((GDALDataType)iType),
                   type_opt[0].c_str()))
        theType=(GDALDataType) iType;
    }
    if(theType==GDT_Unknown)
      theType=inputReader.getDataType();
    if(verbose_opt[0])
      cout << endl << "Output pixel type:  " << GDALGetDataTypeName(theType) << endl;
    if(option_opt.findSubstring("INTERLEAVE=")==option_opt.end()){
      string theInterleave="INTERLEAVE=";
      theInterleave+=inputReader.getInterleave();
      option_opt.push_back(theInterleave);
    }
    outputWriter.open(output_opt[0],inputReader.nrOfCol(),inputReader.nrOfRow(),inputReader.nrOfBand(),theType,inputReader.getImageType(),option_opt);
    for(int iband=0;iband<inputReader.nrOfBand();++iband)
      outputWriter.GDALSetNoDataValue(nodata_opt[0],iband);
    if(description_opt.size())
      outputWriter.setImageDescription(description_opt[0]);

    if(colorTable_opt.size()){
      if(colorTable_opt[0]!="none")
        outputWriter.setColorTable(colorTable_opt[0]);
    }
    else if (inputReader.getColorTable()!=NULL)//copy colorTable from input image
      outputWriter.setColorTable(inputReader.getColorTable());
    
    //if input image is georeferenced, copy projection info to output image
    if(inputReader.isGeoRef()){
      for(int imask=0;imask<mask_opt.size();++imask)
        assert(maskReader[imask].isGeoRef());
    }
    outputWriter.setProjection(inputReader.getProjection());
    double ulx,uly,lrx,lry;
    inputReader.getBoundingBox(ulx,uly,lrx,lry);
    outputWriter.copyGeoTransform(inputReader);
    assert(nodata_opt.size()==masknodata_opt.size());
    if(verbose_opt[0]&&mask_opt.size()){
      for(int iv=0;iv<masknodata_opt.size();++iv)
        cout << masknodata_opt[iv] << "->" << nodata_opt[iv] << endl;
    }

    assert(outputWriter.nrOfCol()==inputReader.nrOfCol());
    // Vector2d<int> lineInput(inputReader.nrOfBand(),inputReader.nrOfCol());
    Vector2d<double> lineInput(inputReader.nrOfBand(),inputReader.nrOfCol());
    Vector2d<short> lineMask(mask_opt.size());
    for(int imask=0;imask<mask_opt.size();++imask)
      lineMask[imask].resize(maskReader[imask].nrOfCol());
    Vector2d<double> lineOutput(outputWriter.nrOfBand(),outputWriter.nrOfCol());
    int irow=0;
    int icol=0;
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    pfnProgress(progress,pszMessage,pProgressArg);
    double oldRowMask=-1;
    for(irow=0;irow<inputReader.nrOfRow();++irow){
      //read line in lineInput buffer
      for(int iband=0;iband<inputReader.nrOfBand();++iband){
        try{
          // inputReader.readData(lineInput[iband],GDT_Int32,irow,iband);
          inputReader.readData(lineInput[iband],GDT_Float64,irow,iband);
        }
        catch(string errorstring){
          cerr << errorstring << endl;
          exit(1);
        }
      }
      double x,y;//geo coordinates
      double colMask,rowMask;//image coordinates in mask image
      for(icol=0;icol<inputReader.nrOfCol();++icol){
        bool masked=false;
        if(mask_opt.size()>1){//multiple masks
          for(int imask=0;imask<mask_opt.size();++imask){
	    inputReader.image2geo(icol,irow,x,y);
	    maskReader[imask].geo2image(x,y,colMask,rowMask);
            if(static_cast<int>(rowMask)!=static_cast<int>(oldRowMask)){
              assert(rowMask>=0&&rowMask<maskReader[imask].nrOfRow());
              try{
                maskReader[imask].readData(lineMask[imask],GDT_Int16,static_cast<int>(rowMask));
              }
              catch(string errorstring){
                cerr << errorstring << endl;
                exit(1);
              }
              oldRowMask=rowMask;
            }
            short ivalue=0;
            if(mask_opt.size()==masknodata_opt.size())//one invalid value for each mask
              ivalue=masknodata_opt[imask];
            else//use same invalid value for each mask
              ivalue=masknodata_opt[0];
            if(lineMask[imask][colMask]==ivalue){
              for(int iband=0;iband<inputReader.nrOfBand();++iband)
                lineInput[iband][icol]=nodata_opt[imask];
              masked=true;
              break;
            }
          }
        }
        else if(mask_opt.size()){//potentially more invalid values for single mask
	  inputReader.image2geo(icol,irow,x,y);
	  maskReader[0].geo2image(x,y,colMask,rowMask);
          if(static_cast<int>(rowMask)!=static_cast<int>(oldRowMask)){
            assert(rowMask>=0&&rowMask<maskReader[0].nrOfRow());
            try{
              maskReader[0].readData(lineMask[0],GDT_Int16,static_cast<int>(rowMask));
            }
            catch(string errorstring){
              cerr << errorstring << endl;
              exit(1);
            }
            oldRowMask=rowMask;
          }
          for(int ivalue=0;ivalue<masknodata_opt.size();++ivalue){
            assert(masknodata_opt.size()==nodata_opt.size());
            if(lineMask[0][colMask]==masknodata_opt[ivalue]){
              for(int iband=0;iband<inputReader.nrOfBand();++iband)
                lineInput[iband][icol]=nodata_opt[ivalue];
              masked=true;
              break;
            }
          }
        }
        for(int iband=0;iband<lineOutput.size();++iband){
          lineOutput[iband][icol]=lineInput[iband][icol];
          if(find(band_opt.begin(),band_opt.end(),iband)!=band_opt.end()){
            if(!masked && codemap.find(lineInput[iband][icol])!=codemap.end()){
              double toValue=codemap[lineInput[iband][icol]];
	      lineOutput[iband][icol]=toValue;
	    }
	  }
        }
      }
      //write buffer lineOutput to output file
      try{
        for(int iband=0;iband<outputWriter.nrOfBand();++iband)
          outputWriter.writeData(lineOutput[iband],GDT_Float64,irow,iband);
      }
      catch(string errorstring){
        cerr << errorstring << endl;
        exit(1);
      }
      //progress bar
      progress=static_cast<float>((irow+1.0)/outputWriter.nrOfRow());
      pfnProgress(progress,pszMessage,pProgressArg);
    }
    inputReader.close();
    outputWriter.close();
  }
}

/**********************************************************************
pkeditogr.cc: program to edit (rename fields) ogr fil
Copyright (C) 2008-2013 Pieter Kempeneers

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
  Optionpk<string> input_opt("i", "input", "Input vector file");
  Optionpk<string> output_opt("o", "output", "Output vector file");
  Optionpk<string> ogrformat_opt("f", "f", "Output OGR file format","ESRI Shapefile");
  Optionpk<string> selectField_opt("select", "select", "select field (combined with like opt)");
  //selectField can also be done via ogr2ogr using the -select option
  Optionpk<string> like_opt("like", "like", "substring(s) to be found in select field. If multiple substrings are provided, feature will be selected if one of them is found (stringent option must be false)");
  Optionpk<bool> stringent_opt("st", "stringent", "string in like option must exactly match to select feature)",false);
  // Optionpk<string> field_opt("as_field", "_field", "output field names (number must exactly match input fields)");
  //renaming fields can also be done via ogr2ogr using the -sql option:
  //ogr2ogr outdataset indataset -sql "SELECT src_field1 AS dst_field1, src_field2 AS dst_field2 FROM sourcelayer"
  Optionpk<long int> setfeature_opt("sf", "sf", "id of feature(s) to set (start from 0)");
  Optionpk<string> setname_opt("sn", "sn", "name(s) of field(s) to set");
  Optionpk<string> setvalue_opt("sv", "sv", "value(s) of field(s) to set");
  Optionpk<string> addname_opt("an", "an", "name(s) of field(s) to add (number must exactly match field types)");
  Optionpk<string> addtype_opt("at", "at", "type(s) of field(s) to add (number must exactly match fieldnames to add", "Real");
  Optionpk<string> addvalue_opt("av", "av", "value(s) of field(s) to add");
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    ogrformat_opt.retrieveOption(argc,argv);
    selectField_opt.retrieveOption(argc,argv);
    like_opt.retrieveOption(argc,argv);
    stringent_opt.retrieveOption(argc,argv);
    // field_opt.retrieveOption(argc,argv);
    addname_opt.retrieveOption(argc,argv);
    addtype_opt.retrieveOption(argc,argv);
    addvalue_opt.retrieveOption(argc,argv);
    setfeature_opt.retrieveOption(argc,argv);
    setname_opt.retrieveOption(argc,argv);
    setvalue_opt.retrieveOption(argc,argv);
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
    std::cerr << "No input file provided (use option -i). Use --help for help information";
      exit(0);//help was invoked, stop processing
  }
  if(output_opt.empty()){
    std::cerr << "No output file provided (use option -o). Use --help for help information";
      exit(0);//help was invoked, stop processing
  }
  if(verbose_opt[0])
    cout << "opening " << input_opt[0] << " for reading " << endl;
  ImgReaderOgr ogrReader(input_opt[0]);
  if(verbose_opt[0])
    cout << "opening " << output_opt[0] << " for writing " << endl;

  OGRFieldType fieldType[addtype_opt.size()];
  int ogr_typecount=11;//hard coded for now!
  for(int it=0;it<addtype_opt.size();++it){
    for(int iType = 0; iType < ogr_typecount; ++iType){
      if( OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType) != NULL
          && EQUAL(OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType),
                   addtype_opt[it].c_str()))
        fieldType[it]=(OGRFieldType) iType;
    }
    if(verbose_opt[0]>1)
      std::cout << std::endl << "field type is: " << OGRFieldDefn::GetFieldTypeName(fieldType[it]) << std::endl;
  }

  //start reading features from the layer
  // if(field_opt.size())
  //   assert(field_opt.size()==ogrReader.getFieldCount());
  unsigned long int ifeature=0;
  if(verbose_opt[0])
    cout << "going through features" << endl << flush;

  ImgWriterOgr ogrWriter(output_opt[0],ogrformat_opt[0]);

  //support multiple layers
  int nlayer=ogrReader.getLayerCount();
  if(verbose_opt[0])
    std::cout << "number of layers: " << nlayer << endl;
      
  for(int ilayer=0;ilayer<nlayer;++ilayer){
    OGRLayer *readLayer=ogrReader.getLayer(ilayer);
    if(verbose_opt[0])
      cout << "reset reading" << endl;
    readLayer->ResetReading();

    OGRLayer *writeLayer=ogrWriter.createLayer(output_opt[0],ogrReader.getProjection(),ogrReader.getGeometryType(ilayer),NULL);
    std::vector<OGRFieldDefn*> readFields;
    std::vector<OGRFieldDefn*> writeFields;
    ogrReader.getFields(readFields,ilayer);
    writeFields=readFields;
    try{
      for(int ifield=0;ifield<readFields.size();++ifield){
	// if(field_opt.size()>ifield)
	//   writeFields[ifield]->SetName(field_opt[ifield].c_str());
	if(verbose_opt[0])
	  std::cout << readFields[ifield]->GetNameRef() << " -> " << writeFields[ifield]->GetNameRef() << std::endl;
	if(writeLayer->CreateField(writeFields[ifield]) != OGRERR_NONE ){
	  ostringstream es;
	  // if(field_opt.size()>ifield)
	  //   es << "Creating field " << field_opt[ifield] << " failed";
	  // else
	  es << "Creating field " << readFields[ifield] << " failed";
	  string errorString=es.str();
	  throw(errorString);
	}
      }
    }
    catch(string errorString){
      std::cerr << errorString << std::endl;
      exit(1);
    }
    if(verbose_opt[0])
      std::cout << "add " << addname_opt.size() << " fields" << std::endl;
    if(addname_opt.size()){
      assert(addname_opt.size()==addtype_opt.size());
      while(addvalue_opt.size()<addname_opt.size())
	addvalue_opt.push_back(addvalue_opt.back());
    }
    for(int iname=0;iname<addname_opt.size();++iname){
      if(verbose_opt[0])
	std::cout << addname_opt[iname] << " " << std::endl;
      ogrWriter.createField(addname_opt[iname],fieldType[iname]);
    }
    if(verbose_opt[0]){
      std::cout << std::endl;
      std::cout << addname_opt.size() << " fields created" << std::endl;
    }
    const char* pszMessage;
    void* pProgressArg=NULL;
    GDALProgressFunc pfnProgress=GDALTermProgress;
    double progress=0;
    OGRFeature *poFeature;
    unsigned long int nfeature=ogrReader.getFeatureCount(ilayer);
    while((poFeature = ogrReader.getLayer(ilayer)->GetNextFeature()) != NULL ){
      if(verbose_opt[0])
	std::cout << "feature " << ifeature << std::endl;
      ++ifeature;
      bool doSelect;
      if(like_opt.empty())
	doSelect=true;
      else{
	assert(selectField_opt.size());
	int fieldIndex=poFeature->GetFieldIndex(selectField_opt[0].c_str());
	string fieldValue=poFeature->GetFieldAsString(fieldIndex);
	if(stringent_opt[0]){
	  if(fieldValue==like_opt[0])
	    doSelect=true;
	  else
	    doSelect=false;
	}
	else{
	  doSelect=false;
	  for(int ilike=0;ilike<like_opt.size();++ilike){
	    if(fieldValue.find(like_opt[ilike])!=std::string::npos){
	      if(verbose_opt[0])
		std::cout << "found " << like_opt[ilike] << " in " << fieldValue << std::endl;
	      doSelect=true;
	      break;
	    }
	  }
	}
      }
      if(!doSelect){
	if(verbose_opt[0])
	  std::cout << "string not found in feature " << ifeature << std::endl;
	continue;
      }
      OGRFeature *poDstFeature = NULL;
      poDstFeature=ogrWriter.createFeature(ilayer);
      if( poDstFeature->SetFrom( poFeature, TRUE ) != OGRERR_NONE ){
	const char* fmt;
	string errorString="Unable to translate feature %d from layer %s.\n";
	fmt=errorString.c_str();
	CPLError( CE_Failure, CPLE_AppDefined,
		  fmt,
		  poFeature->GetFID(), ogrWriter.getLayerName().c_str() );
	OGRFeature::DestroyFeature( poFeature );
	OGRFeature::DestroyFeature( poDstFeature );
      }
      long int fid=poFeature->GetFID();
      poDstFeature->SetFID( poFeature->GetFID() );
      for(int ifeature=0;ifeature<setfeature_opt.size();++ifeature){
	if(fid==setfeature_opt[ifeature]){
	  switch(poDstFeature->GetFieldDefnRef(fid)->GetType()){
	  case(OFTReal):
	    poDstFeature->SetField(setname_opt[ifeature].c_str(),string2type<float>(setvalue_opt[ifeature]));
	    break;
	  case(OFTInteger):
	    poDstFeature->SetField(setname_opt[ifeature].c_str(),string2type<int>(setvalue_opt[ifeature]));
	    break;
	  case(OFTString):
	    poDstFeature->SetField(setname_opt[ifeature].c_str(),setvalue_opt[ifeature].c_str());
	    break;
	  default:
	    std::cerr << "Error: field type not supported" << std::endl;
	    exit(1);
	    break;
	  }
	}
      }

      //set default values for new fields
      if(verbose_opt[0])
	std::cout << "set default values for new fields in feature " << ifeature << std::endl;
      for(int iname=0;iname<addname_opt.size();++iname){
	switch(fieldType[iname]){
	case(OFTReal):
	  if(verbose_opt[0])
	    std::cout << "set field " << addname_opt[iname] << " to default " << string2type<float>(addvalue_opt[iname]) << std::endl;
	  poDstFeature->SetField(addname_opt[iname].c_str(),string2type<float>(addvalue_opt[iname]));
	  break;
	case(OFTInteger):
	  if(verbose_opt[0])
	    std::cout << "set field " << addname_opt[iname] << " to default " << string2type<int>(addvalue_opt[iname]) << std::endl;
	  poDstFeature->SetField(addname_opt[iname].c_str(),string2type<int>(addvalue_opt[iname]));
	  break;
	case(OFTString):
	  if(verbose_opt[0])
	    std::cout << "set field " << addname_opt[iname] << " to default " << addvalue_opt[iname] << std::endl;
	  poDstFeature->SetField(addname_opt[iname].c_str(),addvalue_opt[iname].c_str());
	  break;
	default:
	  std::cerr << "Error: field type not supported" << std::endl;
	  exit(1);
	  break;
	}
      }
      CPLErrorReset();
      if(verbose_opt[0])
	std::cout << "create feature" << std::endl;
      if(ogrWriter.createFeature( poDstFeature,ilayer ) != OGRERR_NONE){
	const char* fmt;
	string errorString="Unable to translate feature %d from layer %s.\n";
	fmt=errorString.c_str();
	CPLError( CE_Failure, CPLE_AppDefined,
		  fmt,
		  poFeature->GetFID(), ogrWriter.getLayerName().c_str() );
	OGRFeature::DestroyFeature( poDstFeature );
      }
      OGRFeature::DestroyFeature( poFeature );
      OGRFeature::DestroyFeature( poDstFeature );
      progress=static_cast<float>(ifeature+1)/nfeature;
      pfnProgress(progress,pszMessage,pProgressArg);
    }
  }
  if(verbose_opt[0])
    std::cout << "replaced " << ifeature << " features" << std::endl;
  ogrReader.close();
  ogrWriter.close();
}

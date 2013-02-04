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

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input image");
  Optionpk<string> output_opt("o", "output", "Output mask file");
  Optionpk<string> field_opt("f", "field", "output field names (number must exactly match input fields)");
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    field_opt.retrieveOption(argc,argv);
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
  if(verbose_opt[0])
    cout << "opening " << input_opt[0] << " for reading " << endl;
  ImgReaderOgr ogrReader(input_opt[0]);
  if(verbose_opt[0])
    cout << "opening " << output_opt[0] << " for writing " << endl;
  //start reading features from the layer
  if(verbose_opt[0])
    cout << "reset reading" << endl;
  ogrReader.getLayer()->ResetReading();
  assert(field_opt.size()==ogrReader.getFieldCount());
  unsigned long int ifeature=0;
  if(verbose_opt[0])
    cout << "going through features" << endl << flush;

  ImgWriterOgr ogrWriter(output_opt[0]);
  OGRLayer* writeLayer=ogrWriter.createLayer(output_opt[0],ogrReader.getProjection(),ogrReader.getGeometryType(),NULL);
  std::vector<OGRFieldDefn*> readFields;
  std::vector<OGRFieldDefn*> writeFields;
  ogrReader.getFields(readFields);
  writeFields=readFields;
  try{
    for(int ifield=0;ifield<readFields.size();++ifield){
      writeFields[ifield]->SetName(field_opt[ifield].c_str());
      if(writeLayer->CreateField(writeFields[ifield]) != OGRERR_NONE ){
        ostringstream es;
        es << "Creating field " << field_opt[ifield] << " failed";
        string errorString=es.str();
        throw(errorString);
      }
    }
  }
  catch(string errorString){
    std::cerr << errorString << std::endl;
    exit(1);
  }
  OGRFeature *poFeature;
  while(true){// (poFeature = imgReaderOgr.getLayer()->GetNextFeature()) != NULL ){
    poFeature = ogrReader.getLayer()->GetNextFeature();
    if( poFeature == NULL )
      break;
    OGRFeature *poDstFeature = NULL;
    poDstFeature=ogrWriter.createFeature();
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
    poDstFeature->SetFID( poFeature->GetFID() );
    OGRFeature::DestroyFeature( poFeature );

    CPLErrorReset();
    if(ogrWriter.createFeature( poDstFeature ) != OGRERR_NONE){
      const char* fmt;
      string errorString="Unable to translate feature %d from layer %s.\n";
      fmt=errorString.c_str();
      CPLError( CE_Failure, CPLE_AppDefined,
                fmt,
                poFeature->GetFID(), ogrWriter.getLayerName().c_str() );
      OGRFeature::DestroyFeature( poDstFeature );
    }
    OGRFeature::DestroyFeature( poDstFeature );
  }
  if(verbose_opt[0])
    cout << "replaced " << ifeature << " features" << endl;
  ogrReader.close();
  ogrWriter.close();
}

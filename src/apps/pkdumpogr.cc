/**********************************************************************
pkdumpogr.cc: dump ogr file to text file or standard output
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
#include <math.h>
#include <string>
#include <fstream>
#include <assert.h>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "pkdumpogr.h"

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i", "input", "Input shape file");
  Optionpk<string> layer_opt("ln", "lname", "Layer name(s) in sample (leave empty to select all)");
  Optionpk<string> output_opt("o", "output", "Output ASCII file");
  Optionpk<string> attribute_opt("n", "name", "names of the attributes to select. Each attribute is stored in a separate band. Default is ALL: write all attributes", "ALL");
  Optionpk<bool> pos_opt("pos","pos","include position (x and y)",false);
  Optionpk<bool> transpose_opt("t","transpose","transpose output (does not work for -n ALL ",false);
  Optionpk<char> fs_opt("fs","fs","field separator.",' ');
  Optionpk<short> verbose_opt("v", "verbose", "verbose (Default: 0)", 0,2);

  fs_opt.setHide(1);
  verbose_opt.setHide(2);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    layer_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    attribute_opt.retrieveOption(argc,argv);
    pos_opt.retrieveOption(argc,argv);
    transpose_opt.retrieveOption(argc,argv);
    fs_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    cout << endl;
    cout << "Usage: pkdumpogr -i input.txt [-o output]" << endl;
    cout << endl;
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  assert(input_opt.size());
  ImgReaderOgr inputReader;
  try{
    inputReader.open(input_opt[0]);
  }
  catch(string errorstring){
    cerr << errorstring << endl;
  }
  ofstream outputFile;
  if(output_opt.size())
    outputFile.open(output_opt[0].c_str(),ios::out);

  inputReader.setFieldSeparator(fs_opt[0]);

  //support multiple layers
  int nlayerRead=inputReader.getDataSource()->GetLayerCount();
  if(verbose_opt[0])
    cout << "number of layers: " << nlayerRead << endl;
      
  for(int ilayer=0;ilayer<nlayerRead;++ilayer){
    OGRLayer *readLayer=inputReader.getLayer(ilayer);
    string currentLayername=readLayer->GetName();
    if(layer_opt.size())
      if(find(layer_opt.begin(),layer_opt.end(),currentLayername)==layer_opt.end())
	continue;
    if(verbose_opt[0])
      cout << "processing layer " << currentLayername << endl;
    // if(layer_opt.size())
    //   cout << " --lname " << currentLayername << endl;
      
    if(attribute_opt[0]=="ALL"){
      attribute_opt.clear();
      OGRFeatureDefn *poFDefn = readLayer->GetLayerDefn();
      for(int iField=0;iField<poFDefn->GetFieldCount();++iField){
	OGRFieldDefn *poFieldDefn = poFDefn->GetFieldDefn(iField);
	std::string fieldname=poFieldDefn->GetNameRef();
	attribute_opt.push_back(fieldname);
      }
    }
    // if(attribute_opt[0]!="ALL"){
      vector<double> xvector;
      vector<double> yvector;
      if(inputReader.getGeometryType()==wkbPoint)
	inputReader.readXY(xvector,yvector);
      Vector2d<std::string> theData(attribute_opt.size());
      for(int ifield=0;ifield<attribute_opt.size();++ifield){
	if(verbose_opt[0])
	  cout << "field: " << ifield << endl;
	theData[ifield].clear();
	inputReader.readData(theData[ifield],OFTReal,attribute_opt[ifield],ilayer,verbose_opt[0]);
      }
      if(verbose_opt[0]){
	std::cout << "number of fields: " << theData.size() << std::endl;
	std::cout << "number of samples: " << theData[0].size() << std::endl;
      }
      if(transpose_opt[0]){
	if(pos_opt[0]&&(inputReader.getGeometryType()==wkbPoint)){
	  if(output_opt.size()){
	    outputFile << "X" << fs_opt[0];
	    for(int isample=0;isample<xvector.size();++isample){
	      outputFile << xvector[isample];
	      if(isample<xvector.size()-1)
		outputFile << fs_opt[0];
	      else
		outputFile << std::endl;
	    }
	    outputFile << "Y" << fs_opt[0];
	    for(int isample=0;isample<yvector.size();++isample){
	      outputFile << yvector[isample];
	      if(isample<yvector.size()-1)
		outputFile << fs_opt[0];
	      else
		outputFile << std::endl;
	    }
	  }
	  else{
	    std::cout << "X" << fs_opt[0];
	    for(int isample=0;isample<xvector.size();++isample){
	      std::cout << xvector[isample];
	      if(isample<xvector.size()-1)
		std::cout << fs_opt[0];
	      else
		std::cout << std::endl;
	    }
	    std::cout << "Y" << fs_opt[0];
	    for(int isample=0;isample<yvector.size();++isample){
	      std::cout << yvector[isample];
	      if(isample<yvector.size()-1)
		std::cout << fs_opt[0];
	      else
		std::cout << std::endl;
	    }
	  }
	}
	for(int ifield=0;ifield<theData.size();++ifield){
	  if(output_opt.size()){
	    outputFile << ifield << fs_opt[0];
	    for(int isample=0;isample<theData[0].size();++isample){
	      outputFile << theData[ifield][isample];
	      if(isample<theData[0].size()-1)
		outputFile << fs_opt[0];
	      else
		outputFile << std::endl;
	    }
	  }
	  else{
	    std::cout << ifield << fs_opt[0];
	    for(int isample=0;isample<theData[0].size();++isample){
	      std::cout << theData[ifield][isample];
	      if(isample<theData[0].size()-1)
		std::cout << fs_opt[0];
	      else
		std::cout << std::endl;
	    }
	  }
	}
      }
      else{
	for(int isample=0;isample<theData[0].size();++isample){
	  if(output_opt.size()){
	    outputFile << isample << fs_opt[0];
	    if(pos_opt[0])
	      outputFile << xvector[isample] << fs_opt[0] << yvector[isample] << fs_opt[0];
	    for(int ifield=0;ifield<theData.size();++ifield){
	      outputFile << theData[ifield][isample];
	      if(ifield<theData.size()-1)
		outputFile << fs_opt[0];
	      else
		outputFile << std::endl;
	    }
	  }
	  else{
	    std::cout << isample << fs_opt[0];
	    if(pos_opt[0])
	      std::cout  << xvector[isample] << fs_opt[0] << yvector[isample] << fs_opt[0];
	    for(int ifield=0;ifield<theData.size();++ifield){
	      std::cout << theData[ifield][isample];
	      if(ifield<theData.size()-1)
		std::cout << fs_opt[0];
	      else
		std::cout << std::endl;
	    }
	  }
	}
      }
      if(output_opt.size())
	outputFile.close();
    // }
    // else{
    //   if(output_opt.size()){
    // 	ofstream outputFile(output_opt[0].c_str(),ios::out);
    // 	outputFile << inputReader;
    // 	outputFile.close();
    //   }
    //   else
    // 	std::cout << inputReader;
    // }
  }
  inputReader.close();
}


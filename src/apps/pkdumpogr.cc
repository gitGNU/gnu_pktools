/**********************************************************************
pkdumpogr.cc: dump ogr file to text file or standard output
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
#include <math.h>
#include <string>
#include <fstream>
#include <assert.h>
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderOgr.h"
#include "algorithms/Histogram.h"
#include "pkdumpogr.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

int main(int argc, char *argv[])
{
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string> input_opt("i", "input", "Input shape file", "");
  Optionpk<string> output_opt("o", "output", "Output ASCII file", "");
  Optionpk<string> attribute_opt("n", "name", "names of the attributes to select. Each attribute is stored in a separate band. Default is ALL: write all attributes", "ALL");
  Optionpk<bool> pos_opt("pos","pos","include position (x and y)",false);
  Optionpk<bool> transpose_opt("t","transpose","transpose output (does not work for -n ALL ",false);
  Optionpk<short> verbose_opt("v", "verbose", "verbose (Default: 0)", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  attribute_opt.retrieveOption(argc,argv);
  pos_opt.retrieveOption(argc,argv);
  transpose_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

  if(version_opt[0]||todo_opt[0]){
    cout << version_opt.getHelp() << endl;
    cout << "todo: " << todo_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }
  if(help_opt[0]){
    cout << "usage: pkdumpogr -i input [OPTIONS]" << endl;
    exit(0);
  }

  ImgReaderOgr imgReader;
  try{
    imgReader.open(input_opt[0]);
  }
  catch(string errorstring){
    cerr << errorstring << endl;
  }
  ofstream outputFile;
  if(output_opt[0]!="")
    outputFile.open(output_opt[0].c_str(),ios::out);

  ImgReaderOgr inputReader(input_opt[0]);
  if(attribute_opt[0]!="ALL"){
    vector<double> xvector;
    vector<double> yvector;
    if(inputReader.getGeometryType()==wkbPoint)
      inputReader.readXY(xvector,yvector);
    Vector2d<double> theData(attribute_opt.size());
    for(int ifield=0;ifield<attribute_opt.size();++ifield){
      if(verbose_opt[0])
        cout << "field: " << ifield << endl;
      theData[ifield].clear();
      inputReader.readData(theData[ifield],OFTReal,attribute_opt[ifield],0,verbose_opt[0]);
    }
    if(verbose_opt[0]){
      std::cout << "number of fields: " << theData.size() << std::endl;
      std::cout << "number of samples: " << theData[0].size() << std::endl;
    }
    if(transpose_opt[0]){
      if(pos_opt[0]&&(inputReader.getGeometryType()==wkbPoint)){
        if(output_opt[0]!=""){
          outputFile << "X" << " ";
          for(int isample=0;isample<xvector.size();++isample){
            outputFile << xvector[isample];
            if(isample<xvector.size()-1)
              outputFile << " ";
            else
              outputFile << std::endl;
          }
          outputFile << "Y" << " ";
          for(int isample=0;isample<yvector.size();++isample){
            outputFile << yvector[isample];
            if(isample<yvector.size()-1)
              outputFile << " ";
            else
              outputFile << std::endl;
          }
        }
        else{
          std::cout << "X" << " ";
          for(int isample=0;isample<xvector.size();++isample){
            std::cout << xvector[isample];
            if(isample<xvector.size()-1)
              std::cout << " ";
            else
              std::cout << std::endl;
          }
          std::cout << "Y" << " ";
          for(int isample=0;isample<yvector.size();++isample){
            std::cout << yvector[isample];
            if(isample<yvector.size()-1)
              std::cout << " ";
            else
              std::cout << std::endl;
          }
        }
      }
      for(int ifield=0;ifield<theData.size();++ifield){
        if(output_opt[0]!=""){
          outputFile << ifield << " ";
          for(int isample=0;isample<theData[0].size();++isample){
            outputFile << theData[ifield][isample];
            if(isample<theData[0].size()-1)
              outputFile << " ";
            else
              outputFile << std::endl;
          }
        }
        else{
          std::cout << ifield << " ";
          for(int isample=0;isample<theData[0].size();++isample){
            std::cout << theData[ifield][isample];
            if(isample<theData[0].size()-1)
              std::cout << " ";
            else
              std::cout << std::endl;
          }
        }
      }
    }
    else{
      for(int isample=0;isample<theData[0].size();++isample){
        if(output_opt[0]!=""){
          outputFile << isample << " ";
          if(pos_opt[0])
            outputFile << xvector[isample] << " " << yvector[isample] << " ";
          for(int ifield=0;ifield<theData.size();++ifield){
            outputFile << theData[ifield][isample];
            if(ifield<theData.size()-1)
              outputFile << " ";
            else
              outputFile << std::endl;
          }
        }
        else{
          std::cout << isample << " ";
          if(pos_opt[0])
            std::cout  << xvector[isample] << " " << yvector[isample] << " ";
          for(int ifield=0;ifield<theData.size();++ifield){
            std::cout << theData[ifield][isample];
            if(ifield<theData.size()-1)
              std::cout << " ";
            else
              std::cout << std::endl;
          }
        }
      }
    }
    if(output_opt[0]!="")
      outputFile.close();
  }
  else{
    if(output_opt[0]!=""){
      ofstream outputFile(output_opt[0].c_str(),ios::out);
      outputFile << imgReader;
      outputFile.close();
    }
    else
      std::cout << imgReader;
  }
  imgReader.close();
}


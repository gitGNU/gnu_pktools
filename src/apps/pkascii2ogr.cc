/**********************************************************************
pkascii2ogr.cc: program to create vector (ogr) points polygon) from text file
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
#include <string>
#include <fstream>
#include <assert.h>
#include "base/Optionpk.h"
#include "imageclasses/ImgWriterOgr.h"

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
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<string> input_opt("i","input","Input ASCII file");
  Optionpk<string> output_opt("o", "output", "Output file");
  Optionpk<short> colX_opt("x", "x", "column number of x (0)", 0);
  Optionpk<short> colY_opt("y", "y", "column number of y (1)", 1);
  Optionpk<bool> polygon_opt("l", "line", "create OGRPolygon as geometry instead of points.  Fields are taken from first point and polygon is automatically closed (no need to repeat first point at last line). (false: use OGRPoint)", false);
  Optionpk<string> fname_opt("n", "name", "Field names for the columns in the input ascii file");
  Optionpk<string> ftype_opt("ot", "ot", "Field type (Real, Integer, String) for each of the fields as defined by name","Real");
  Optionpk<string> itype_opt("of", "of", "image type string", "ESRI Shapefile");
  Optionpk<string> projection_opt("p", "projection", "use EPSG:<code> or Wkt string", "EPSG:4326");
  Optionpk<char> fs_opt("fs","fs","field separator.",' ');
  Optionpk<int> verbose_opt("v", "verbose", "verbose (0)", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);
  input_opt.retrieveOption(argc,argv);
  output_opt.retrieveOption(argc,argv);
  colX_opt.retrieveOption(argc,argv);
  colY_opt.retrieveOption(argc,argv);
  polygon_opt.retrieveOption(argc,argv);
  fname_opt.retrieveOption(argc,argv);
  ftype_opt.retrieveOption(argc,argv);
  itype_opt.retrieveOption(argc,argv);
  projection_opt.retrieveOption(argc,argv);
  fs_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);

  if(version_opt[0]){
    cout << version_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }

  if(help_opt[0]){
    cout << "usage: pksensormodel -i asciifile [OPTIONS]" << endl;
    exit(0);
  }

  string theProjection;
  theProjection=projection_opt[0];
  int ogr_typecount=11;//hard coded for now!
  while(ftype_opt.size()<fname_opt.size())
    ftype_opt.push_back(ftype_opt[0]);
  // vector<string> fname(fname_opt.size());
  vector<OGRFieldType> ftype(ftype_opt.size());
  if(verbose_opt[0])
    cout << "field types can be: ";
  for(int ifield=0;ifield<fname_opt.size();++ifield){
    // fname[ifield]=fname_opt[ifield];
    for(int iType = 0; iType < ogr_typecount; ++iType){
      if(!ifield&&verbose_opt[0])
        cout << " " << OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType);
      if( OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType) != NULL
          && EQUAL(OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType),
                   ftype_opt[ifield].c_str()))
        ftype[ifield]=(OGRFieldType) iType;
    }
  }
  //todo: what if unknown?
  if(verbose_opt[0]){
    cout << endl << "field types are: ";
    for(int ifield=0;ifield<ftype.size();++ifield)
      cout << OGRFieldDefn::GetFieldTypeName(ftype[ifield]) << " ";
    cout << endl;
  }
  
  ImgWriterOgr imgWriter;
  imgWriter.open(output_opt[0]);
  try{
    if(polygon_opt[0])
      imgWriter.ascii2ogr(input_opt[0], "New Layer", fname_opt, ftype, colX_opt[0], colY_opt[0], theProjection, wkbPolygon, fs_opt[0]);
    else
      imgWriter.ascii2ogr(input_opt[0], "New Layer", fname_opt, ftype, colX_opt[0], colY_opt[0], theProjection, wkbPoint, fs_opt[0]);
  }
  catch(string errorString){
    cout << errorString << endl;
  }
  imgWriter.close();
}

/**********************************************************************
pkegcs.cc: Utility for raster files in European Grid Coordinate System
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
#include "base/Optionpk.h"
#include "imageclasses/ImgReaderGdal.h"
#include "algorithms/Egcs.h"

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
  Optionpk<string> image_opt("i","image","input image to analyse","");
  Optionpk<unsigned short>  band_opt("b", "band", "Band specific information", 0);
  Optionpk<string> cell2bb_opt("c2b","cell2bb","convert cell code to geo coordinates of boundingbox (e.g. 32-AB)","");
  Optionpk<string> cell2mid_opt("c2m","cell2mid","convert cell code to centre in geo coordinates (e.g. 32-AB)","");
  Optionpk<bool> refpixel_opt("\0", "ref", "get reference pixel (lower left corner of centre of gravity pixel)", false);
  Optionpk<double> maskValue_opt("m", "mask", "mask value(s) for no data to calculate reference pixel in image (Default is 0)",0);
  Optionpk<int> dx_opt("dx","dx","resolution (default is 250m)",250);
  Optionpk<bool> geo2cell_opt("g2c", "geo2cell", "get cell code for coordinates in x_opt and y_opt given the resolution in dx_opt", false);
  Optionpk<double> x_opt("x","x","x coordinate in epsg:3035",0);
  Optionpk<double> y_opt("y","y","y coordinate in epsg:3035",0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);

  if(version_opt[0]){
    cout << version_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }

  image_opt.retrieveOption(argc,argv);
  band_opt.retrieveOption(argc,argv);
  cell2bb_opt.retrieveOption(argc,argv);
  cell2mid_opt.retrieveOption(argc,argv);
  geo2cell_opt.retrieveOption(argc,argv);
  refpixel_opt.retrieveOption(argc,argv);
  maskValue_opt.retrieveOption(argc,argv);
  dx_opt.retrieveOption(argc,argv);
  x_opt.retrieveOption(argc,argv);
  y_opt.retrieveOption(argc,argv);
  
  Egcs egcs;
  if(cell2bb_opt[0]!=""){
    int theULX, theULY, theLRX, theLRY;
    egcs.setLevel(egcs.cell2level(cell2bb_opt[0]));
    egcs.cell2bb(cell2bb_opt[0],theULX,theULY,theLRX,theLRY);
    cout << setprecision(12) << "--ulx=" << theULX << " --uly=" << theULY << " --lrx=" << theLRX << " --lry=" << theLRY << endl;
  }
  if(cell2mid_opt[0]!=""){
    double midX, midY;
    egcs.setLevel(egcs.cell2level(cell2mid_opt[0]));
    egcs.cell2mid(cell2mid_opt[0],midX,midY);
    cout << setprecision(12) << "-x=" << midX << " -y=" << midY << endl;
  }
  if(geo2cell_opt[0]){
    egcs.setLevel(egcs.res2level(dx_opt[0]));
    cout << egcs.geo2cell(x_opt[0],y_opt[0]) << endl;
  }
  if(image_opt[0]!=""){
    ImgReaderGdal imgReader;
    imgReader.open(image_opt[0]);
    if(refpixel_opt[0]){
      assert(band_opt[0]<imgReader.nrOfBand());
      for(int inodata=0;inodata<maskValue_opt.size();++inodata)
        imgReader.pushNoDataValue(maskValue_opt[inodata]);
      // if(verbose_opt[0]){
      //   vector<double> noData;
      //   imgReader.getNoDataValues(noData,band_opt[0]);
      //   cout << "number of no data values: " << noData.size() << endl;
      // }
      double refX,refY;
      //get centre of reference (centre of gravity) pixel in image
      imgReader.getRefPix(refX,refY,band_opt[0]);
      cout << setprecision(12) << "--x " << refX << " --y " << refY << endl;
      egcs.setLevel(egcs.res2level(imgReader.getDeltaX()));
      // unsigned short theLevel=egcs.getLevel(imgReader.getDeltaX());
      // egcs.setLevel(theLevel);
      cout << "cell code at level " << egcs.getLevel() << " (resolution is " << egcs.getResolution() << "): " << egcs.geo2cell(refX,refY) << endl;
    }
  }
}

/**********************************************************************
pktestProspect: example program how to use class Prospect
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
#include <iostream>
#include <string>
#include "base/Optionpk.h"
#include "algorithms/Filter.h"
#include "fileclasses/FileReaderAscii.h"
#include "Prospect.h"

int main(int argc, char *argv[])
{
  Optionpk<double> N_opt("N","N","parameter N for Prospect model",1.2);
  Optionpk<double> Cab_opt("Cab","Cab","parameter Cab for Prospect model (g.cm-2)",30);
  Optionpk<double> Car_opt("Car","Car","parameter Car for Prospect model (g.cm-2)",10);
  Optionpk<double> Cbrown_opt("Cbrown","Cbrown","parameter Cbrown for Prospect model (arbitrary units)",0);
  Optionpk<double> Cw_opt("Cw","Cw","parameter Cw for Prospect model (cm)",0.015);
  Optionpk<double> Cm_opt("Cm","Cm","parameter Cm for Prospect model (g.cm-2)",0.009);
  Optionpk<double> step_opt("s","step","step size for interpolation",0.01);
  Optionpk<double> fwhm_opt("fwhm", "fwhm", "list of full width half to apply spectral filtering (-fwhm band1 -fwhm band2 ...)");
  Optionpk<double> wavelengthIn_opt("win", "wavelengthIn", "list of wavelengths in input spectrum (-w band1 -w band2 ...)");
  Optionpk<double> wavelengthOut_opt("wout", "wavelengthOut", "list of wavelengths in output spectrum (-w band1 -w band2 ...)");
  Optionpk<short> verbose_opt("v","verbose","verbose mode",0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=N_opt.retrieveOption(argc,argv);
    doProcess=Cab_opt.retrieveOption(argc,argv);
    doProcess=Car_opt.retrieveOption(argc,argv);
    doProcess=Cbrown_opt.retrieveOption(argc,argv);
    doProcess=Cw_opt.retrieveOption(argc,argv);
    doProcess=Cm_opt.retrieveOption(argc,argv);
    doProcess=fwhm_opt.retrieveOption(argc,argv);
    doProcess=wavelengthIn_opt.retrieveOption(argc,argv);
    doProcess=wavelengthOut_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(std::string predefinedString){//command line option contained license or version
    std::cout << predefinedString << std::endl;//report the predefined string to stdout
    exit(0);//stop processing
  }
  if(!doProcess){//command line option contained help option
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;//provide extra details for help to the user
    exit(0);//stop processing
  }

  prospect::Prospect prospectModel;
  prospectModel.setN(N_opt[0]);
  prospectModel.setCab(Cab_opt[0]);
  prospectModel.setCar(Car_opt[0]);
  prospectModel.setCbrown(Cbrown_opt[0]);
  prospectModel.setCw(Cw_opt[0]);
  prospectModel.setCm(Cm_opt[0]);
  vector<double> leafReflectance;
  vector<double> leafTransmittance;
  prospectModel.getLeafSpectrum(leafReflectance,leafTransmittance);
  assert(leafReflectance.size()==leafTransmittance.size());
  for(int index=0;index<leafReflectance.size();++index)
    std::cout << 400+index << " " << leafReflectance[index] << " " << leafTransmittance[index] << std::endl;
}

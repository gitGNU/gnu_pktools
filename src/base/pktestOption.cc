/**********************************************************************
pktestOption: example program how to use class Optionpk pktestOption.cc 
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

int main(int argc, char *argv[])
{
  Optionpk<std::string> foo_opt("f","foo","command line option **foo** of type string can be invoked with either short (f) or long (foo) option","defaultString");
  Optionpk<int> bar_opt("\0","bar","command line option **bar** of type int has no short option");//bar will only be visible in long help (hide=1)
  Optionpk<bool> easterEgg_opt("egg","egg","this help information is useless");//this option will not be shown in help (hide=2)

  bar_opt.setHide(1);//option only visible with long help (--help)
  easterEgg_opt.setHide(2);//hidden option

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=foo_opt.retrieveOption(argc,argv);
    bar_opt.retrieveOption(argc,argv);
    easterEgg_opt.retrieveOption(argc,argv);
    }
  catch(std::string predefinedString){//command line option contained license or version
    std::cout << predefinedString << std::endl;//report the predefined string to stdout
    exit(0);//stop processing
  }
  if(!doProcess){//command line option contained help option
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;//provide extra details for help to the user
    exit(0);//stop processing
  }
  std::cout << "foo: ";
  for(int ifoo=0;ifoo<foo_opt.size();++ifoo){
    std::cout << foo_opt[ifoo] << " ";
  }
  std::cout << std::endl;
  std::cout << foo_opt << std::endl;//short cut for code above

  if(bar_opt[0]>0)
    std::cout << "long option for bar was used with a positive value" << std::endl;
  
  if(easterEgg_opt[0])
    std::cout << "How did you find this option -egg or --egg? Not through the help info!" << std::endl;
}

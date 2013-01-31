#include <iostream>
#include <string>
#include "base/Optionpk.h"

int main(int argc, char *argv[])
{
  Optionpk<std::string> foo_opt("f","foo","command line option **foo** of type string can be invoked with either short (f) or long (foo) option","defaultString");
  Optionpk<int> bar_opt("\0","bar","command line option **bar** of type int has no short option",false,1);//bar will only be visible in long help (hide=1)
  Optionpk<bool> easterEgg_opt("egg","egg","this help information is useless",false,2);//this option will not be shown in help (hide=2)

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

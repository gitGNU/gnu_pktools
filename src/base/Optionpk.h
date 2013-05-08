/**********************************************************************
Optionpk.h: class to handle command line options (inherits from stl vector class)
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
#ifndef _OPTIONPK_H_
#define _OPTIONPK_H_

#include <vector>
#include <string>
#include <cstdlib>
#include <assert.h>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <typeinfo>
#include "gdal/ogr_feature.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

using namespace std;

///throw this class when syntax error in command line option
class BadConversion : public runtime_error {
 public:
 BadConversion(string const& s)
   : runtime_error(s)
    { }
};

///convert command line option to value of the defined type, throw exception in case of failure
template<typename T> inline T string2type(std::string const& s){
  std::istringstream i(s);
  T x;
  if (!(i >> x) )
     throw BadConversion(s);
  return x;
}

///convert command line option to value of the defined type, throw if entire string could not get converted
template<typename T> inline T string2type(std::string const& s,bool failIfLeftoverChars){
  std::istringstream i(s);
  char c;
  T x;
  if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
     throw BadConversion(s);
  return x;
}

///specialization for string
template<> inline std::string string2type(std::string const& s){
  return s;
}

///specialization for OGRFieldType
template<> inline OGRFieldType string2type(std::string const& s){
  OGRFieldType ftype;
  int ogr_typecount=11;//hard coded for now!
  for(int iType = 0; iType < ogr_typecount; ++iType){
    if( OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType) != NULL
        && EQUAL(OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType),s.c_str()))
      ftype=(OGRFieldType) iType;
  }
  return ftype;
}

///serialization for help or to dump option values to screen in verbose mode
template<typename T> inline std::string type2string(T const& value){
  std::ostringstream oss;
  oss << value;
  return oss.str();
}

///specialization for bool
template<> inline std::string type2string(bool const& value){
  if(value)
    return("true");
  else
    return("false");
}

///specialization for string
template<> inline std::string type2string(std::string const& value){
  // if(value.empty())
  //   return("<empty string>");
  // else
    return(value);
}

///specialization for float
template<> inline std::string type2string(float const& value){
  std::ostringstream oss;
  // oss.precision(1);
  // oss.setf(ios::fixed);
  oss << value;
  return oss.str();
}

///specialization for double
template<> inline std::string type2string(double const& value){
  std::ostringstream oss;
  // oss.precision(1);
  //  oss.setf(ios::fixed);
  oss << value;
  return oss.str();
}

/**
Class to implement command line options. With the constructor you can define an option, in both short `-` and long `--` format, of a specific type, help information and a default value.\n
This class inherits from std::vector, so the option variable is a vector, supporting multiple inputs for the same option (e.g., --input file1 [--input file2 ...].
Several command line option formats are supported: 
- `-shortOption value`
- `-shortOption=value`
- `--longOption value`
- `--longOption=value`
- `-shortOption` (no value for boolean options, which are automatically set by invoking the option)
- `--longOption` (no value for boolean options, which are automatically set by invoking the option)

Option names should have regular characters and no white space in them. Some names are reserved and can not be used either:
- short option `h` or long option `help`: shows help info
- long option `license`: shows license info
- long option `version`: shows current version of pktools
- long option `doxygen`: shows help info in table format, ready to be included in doxygen

A call to member function \ref retrieveOption reads the command line arguments and initializes the object (vector). Make sure to call this member function before using the option object in your main program (or a segmentation error due to an un-initialized vector will occur).

All calls to retrieveOption should reside in a try{} block. If one of the reserved options
- `license`
- `version`
is used, an exception of type std::string is thrown. This can be caught with a catch(string predefinedString) right after the try block, where the message can be sent to stdout and the program can be ended.

Similarly, if help is invoked with the short option `-h` or long option `--help`, the main program is informed by the return value `false` of \ref retrieveOption (for any option). An example how to use Optionpk is shown in \ref pktestOption.cc
**/

template<class T> class Optionpk : public std::vector <T>
{
public:
  ///default constructor
  Optionpk();
  ///constructor for option without default value
  Optionpk(const std::string& shortName, const std::string& longName, const std::string& helpInfo);
  ///constructor for option with default value. Option can be hidden for help info
  Optionpk(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const T& defaultValue, short hide=0);
  ///default destructor
  ~Optionpk();
  ///set help information
  void setHelp(const std::string& helpInfo){m_help=helpInfo;};
  bool retrieveOption(int argc, char ** argv);
  template<class T1> friend ostream& operator<<(ostream & os, const Optionpk<T1>& theOption);

  void setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo);
  void setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const T& defaultValue, short hide);
  void setDefault(const T& defaultValue);
  std::string getDefaultValue() const {return m_defaultValue;};
  void setShortName(const std::string& shortName);
  void setLongName(const std::string& longName);
  std::string getShortName() const {return m_shortName;};
  std::string getLongName() const {return m_longName;};

  std::string getHelp() const {return m_help;};
  static std::string getGPLv3License(){
    return static_cast<std::string>("\n\
    This program is free software: you can redistribute it and/or modify\n\
    it under the terms of the GNU General Public License as published by\n\
    the Free Software Foundation, either version 3 of the License, or\n\
    (at your option) any later version.\n\
    \n\
    This program is distributed in the hope that it will be useful,\n\
    but WITHOUT ANY WARRANTY; without even the implied warranty of\n\
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n\
    GNU General Public License for more details.\n\
                                          \n\
    You should have received a copy of the GNU General Public License\n\
    along with this program.  If not, see <http://www.gnu.org/licenses/>.\n");};
  std::vector<std::string>::const_iterator findSubstring(const std::string& argument) const;
private:
  bool hasArgument() const {return m_hasArgument;};//all options except bools should have arguments
  bool hasShortOption() const {return m_shortName.compare("\0");};
  bool hasLongOption() const {return m_longName.compare("\0");};
  std::string usage() const;
  std::string usageDoxygen() const;

  std::string m_shortName;
  std::string m_longName;
  std::string m_help;
  bool m_hasArgument;
  T m_defaultValue;
  bool m_hasDefault;
  short m_hide;
};

template<class T1> ostream& operator<<(ostream& os, const Optionpk<T1>& theOption)
{
  os << theOption.getLongName() << ": ";
  for(int index=0;index<theOption.size();++index)
    os << type2string<T1>(theOption[index]) << " ";
  os << std::endl;
  return os;
}

template<class T> Optionpk<T>::Optionpk() 
: m_hasDefault(false)
{
}

/**
constructor without default value\n
shortName is option invoked with `-`\n
longName is option invoked with `--`\n
helpInfo is the help message that is shown when option -h or --help is invoked\n
**/
template<class T> Optionpk<T>::Optionpk(const std::string& shortName, const std::string& longName, const std::string& helpInfo)
: m_hasDefault(false)
{
  setAll(shortName,longName,helpInfo);
}

/**
constructor with default value.\n
shortName is option invoked with `-`\n
longName is option invoked with `--`\n
helpInfo is the help message that is shown when option -h or --help is invoked\n
defaultValue is default value of the option (first value of vector: option[0])\n
hide=0 : option is visible for in both short (`-h`) and long (`--help`) help. Typical use: mandatory options\n
hide=1 : option is only visible in long help (`--help`). Typical use: expert options\n
hide=2 : option is hidden for user. Typical use: Easter eggs or options only known to author
**/
template<class T> Optionpk<T>::Optionpk(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const T& defaultValue, short hide)
{
  setAll(shortName,longName,helpInfo,defaultValue, hide);
}

template<class T> std::string Optionpk<T>::usage() const
{
  std::ostringstream helpss;
  std::string shortOption=m_shortName;
  std::string longOption=m_longName;
  shortOption.insert(0,"-");
  longOption.insert(0,"--");
  if(hasShortOption())
    helpss << "   " << setiosflags(ios::left) << setw(4) << shortOption;
  else
    helpss << "   " << setiosflags(ios::left) << setw(4) << " ";
  if(hasLongOption())
    helpss << "   " << setiosflags(ios::left) << setw(20) << longOption;
  else
    helpss << "   " << setiosflags(ios::left) << setw(20) << " ";
  helpss << "   " << m_help;
  if(m_hasDefault)
    helpss << " (default: " << type2string<T>(m_defaultValue) << ")";
  return helpss.str();
}

template<class T> std::string Optionpk<T>::usageDoxygen() const
{
  std::ostringstream helpss;
  std::string shortOption=m_shortName;
  std::string longOption=m_longName;

  if(hasShortOption())
    helpss << " | " << setiosflags(ios::left) << setw(6) << shortOption << " | ";
  else
    helpss << " | " << setiosflags(ios::left) << "       | ";
  if(hasLongOption())
    helpss << setiosflags(ios::left) << setw(20) << longOption << " | ";
  else
    helpss << setiosflags(ios::left) << "                     | ";
  helpss << setiosflags(ios::left) << setw(4) << typeid(T).name() << " | ";
  if(m_hasDefault)
    helpss <<setiosflags(ios::left) << setw(5) << type2string<T>(m_defaultValue) << " |";
  else
    helpss << "      |";
  helpss << m_help << " | ";

  return helpss.str();
}

template<class T> void Optionpk<T>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=true;
  m_help=helpInfo;
  m_hide=0;
}

template<class T> void Optionpk<T>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const T& defaultValue, short hide)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=true;
  m_help=helpInfo;
  m_defaultValue=defaultValue;
  m_hasDefault=true;
  m_hide=hide;
}

///specialization for bool
template<> void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo);

template<> void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=false;
  m_help=helpInfo;
  m_hide=0;
}

///specialization for bool
template<> void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const bool& defaultValue, short hide);

///specialization for bool
template<> void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const bool& defaultValue, short hide)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=false;
  m_help=helpInfo;
  m_defaultValue=defaultValue;
  m_hasDefault=true;
  m_hide=hide;
}

///specialization for bool
template<> Optionpk<bool>::Optionpk(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const bool& defaultValue, short hide)
{
  setAll(shortName,longName,helpInfo,defaultValue, hide);
}

template<class T> Optionpk<T>::~Optionpk() 
{
}

/**
make sure to call this function first before using the option in main program (or segmentation fault will occur...)
**/
template<class T> bool Optionpk<T>::retrieveOption(int argc, char **argv){ 
  bool noHelp=true;//return value, alert main program that hard coded option (help, version, license, doxygen) was invoked
  std::string helpStringShort="-h";//short option for help (hard coded)
  std::string helpStringLong="--help";//long option for help (hard coded)
  std::string helpStringDoxygen="--doxygen";//option to create table of options ready to use for doxygen
  std::string versionString="--version";//option to show current version
  std::string licenseString="--license";//option to show current version
  for(int i = 1; i < argc; ++i ){
    std::string currentArgument;
    std::string currentOption=argv[i];
    std::string shortOption=m_shortName;
    std::string longOption=m_longName;
    shortOption.insert(0,"-");
    longOption.insert(0,"--");
    size_t foundEqual=currentOption.rfind("=");
    if(foundEqual!=std::string::npos){
      currentArgument=currentOption.substr(foundEqual+1);
      currentOption=currentOption.substr(0,foundEqual);
    }
    if(!helpStringShort.compare(currentOption)){
      if(m_hide<1)
        std::cout << usage() << std::endl;
      noHelp=false;
    }
    else if(!helpStringLong.compare(currentOption)){
      if(m_hide<2)
        std::cout << usage() << std::endl;
      noHelp=false;
    }
    else if(!helpStringDoxygen.compare(currentOption)){
      if(m_hide<2)
        std::cout << usageDoxygen() << std::endl;
      noHelp=false;
    }
    else if(!versionString.compare(currentOption)){
      std::string theVersion="version ";
      theVersion+=VERSION;
      theVersion+=", Copyright (C) Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n \
   This is free software, and you are welcome to redistribute it\n      \
   under certain conditions; use option --license for details.";
      throw(theVersion);//no need to continue registering (break prevents from multiplication of version info)
    }
    else if(!licenseString.compare(currentOption)){
      throw(getGPLv3License());
    }
    if(hasShortOption()&&!(shortOption.compare(currentOption))){//for -option
      if(foundEqual!=std::string::npos)
        this->push_back(string2type<T>(currentArgument));
      else if(m_hasArgument && i < argc-1)
        this->push_back(string2type<T>(argv[++i]));
      else
        this->push_back(string2type<T>("1"));
    }
    else if(hasLongOption()&&!(longOption.compare(currentOption))){//for --option
      if(foundEqual!=std::string::npos)
        this->push_back(string2type<T>(currentArgument));
      else if(m_hasArgument && i < argc-1)
        this->push_back(string2type<T>(argv[++i]));
      else
        this->push_back(string2type<T>("1"));
    }
  }
  if(!(this->size())&&m_hasDefault)//only set default value if no options were given
    this->push_back(m_defaultValue);
  return(noHelp);
}

//find a substring in string option (e.g., option is of type -co INTERLEAVE=BAND)
template<> std::vector<std::string>::const_iterator Optionpk<std::string>::findSubstring(const std::string& argument) const{
  std::vector<std::string>::const_iterator opit=this->begin();
  while(opit!=this->end()){
    if(opit->find(argument)!=std::string::npos)
      break;
    ++opit;
  }
  return opit;
}

#endif

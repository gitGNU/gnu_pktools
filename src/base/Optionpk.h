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
#include "ogr_feature.h"

using namespace std;

class BadConversion : public runtime_error {
 public:
 BadConversion(string const& s)
   : runtime_error(s)
    { }
};

template<typename T> inline T string2type(string const& s,bool failIfLeftoverChars){
  istringstream i(s);
  char c;
  T x;
  if (!(i >> x) || (failIfLeftoverChars && i.get(c)))
     throw BadConversion(s);
  return x;
}

template<typename T> inline T string2type(string const& s){
  istringstream i(s);
  T x;
  if (!(i >> x) )
     throw BadConversion(s);
  return x;
}

//specialization for string
template<> inline string string2type(string const& s){
  return s;
}

//specialization for OGRFieldType
template<> inline OGRFieldType string2type(string const& s){
  OGRFieldType ftype;
  int ogr_typecount=11;//hard coded for now!
  for(int iType = 0; iType < ogr_typecount; ++iType){
    if( OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType) != NULL
        && EQUAL(OGRFieldDefn::GetFieldTypeName((OGRFieldType)iType),s.c_str()))
      ftype=(OGRFieldType) iType;
  }
  return ftype;
}

template<typename T> inline string type2string(T const& value){
  ostringstream oss;
  oss << value;
  return oss.str();
}

//specialization for bool
template<> inline string type2string(bool const& value){
  if(value)
    return("true");
  else
    return("false");
}

//specialization for string
template<> inline string type2string(string const& value){
  if(value.empty())
    return("<empty string>");
  else
    return(value);
}

//specialization for float
template<> inline string type2string(float const& value){
  ostringstream oss;
  // oss.precision(1);
  // oss.setf(ios::fixed);
  oss << value;
  return oss.str();
}

//specialization for double
template<> inline string type2string(double const& value){
  ostringstream oss;
  // oss.precision(1);
  //  oss.setf(ios::fixed);
  oss << value;
  return oss.str();
}

template<class T> class Optionpk : public vector <T>
{
public:
  Optionpk();
  Optionpk(const string& shortName, const string& longName, const string& helpInfo);
  Optionpk(const string& shortName, const string& longName, const string& helpInfo,const T& defaultValue);
  ~Optionpk();
  void setHelp(const string& helpInfo){m_help=helpInfo;};
  string usage() const;
  static string getGPLv3License(){
    return static_cast<string>("\n\
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
  string getHelp() const {return m_help;};
  int retrieveOption(int argc, char ** argv);
  std::vector<string>::const_iterator findSubstring(string argument) const;
  template<class T1> friend ostream& operator<<(ostream & os, const Optionpk<T1>& theOption);
private:
  void setAll(const string& shortName, const string& longName, const string& helpInfo);
  void setAll(const string& shortName, const string& longName, const string& helpInfo,const T& defaultValue);
  void setDefault(const T& defaultValue);
  string getDefaultValue() const {return m_defaultValue;};
  void setShortName(const string& shortName);
  void setLongName(const string& longName);
  string getShortName() const {return m_shortName;};
  string getLongName() const {return m_longName;};
  bool hasArgument() const {return m_hasArgument;};//all options except bools should have arguments
  bool hasShortOption() const {return m_shortName.compare("\0");};
  bool hasLongOption() const {return m_longName.compare("\0");};

  string m_shortName;
  string m_longName;
  string m_help;
  bool m_hasArgument;
  T m_defaultValue;
  bool m_hasDefault;
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

template<class T> Optionpk<T>::Optionpk(const string& shortName, const string& longName, const string& helpInfo)
  : m_hasDefault(false)
{
  setAll(shortName,longName,helpInfo);
}

  template<class T> Optionpk<T>::Optionpk(const string& shortName, const string& longName, const string& helpInfo,const T& defaultValue)
{
  setAll(shortName,longName,helpInfo,defaultValue);
}

template<class T> string Optionpk<T>::usage() const
{
  ostringstream helpss;
  string shortOption=m_shortName;
  string longOption=m_longName;
  if(hasShortOption())
    helpss << " | " << setiosflags(ios::left) << setw(6) << shortOption << " | ";
  else 
    helpss << " | " << setiosflags(ios::left) << setw(6) << " | ";
  if(hasLongOption())
    helpss << setiosflags(ios::left) << setw(20) << longOption << " | ";
  else
    helpss << setiosflags(ios::left) << setw(20) << " | ";
  helpss << setiosflags(ios::left) << setw(4) << typeid(T).name() << " | ";
  helpss << m_help << " | ";
  if(m_hasDefault)
    helpss << type2string<T>(m_defaultValue) << "|";
  else
    helpss << "|";
  return helpss.str();
  /* ostringstream helpss; */
  /* string shortOption=m_shortName; */
  /* string longOption=m_longName; */
  /* shortOption.insert(0,"-"); */
  /* longOption.insert(0,"--"); */
  /* if(hasShortOption()) */
  /*   helpss << "   " << setiosflags(ios::left) << setw(4) << shortOption; */
  /* else  */
  /*   helpss << "   " << setiosflags(ios::left) << setw(4) << " "; */
  /* if(hasLongOption()) */
  /*   helpss << "   " << setiosflags(ios::left) << setw(20) << longOption; */
  /* else */
  /*   helpss << "   " << setiosflags(ios::left) << setw(20) << " "; */
  /* helpss << "   " << m_help; */
  /* if(m_hasDefault) */
  /*   helpss << " (default: " << type2string<T>(m_defaultValue) << ")"; */
  /* return helpss.str(); */
}

template<class T> void Optionpk<T>::setAll(const string& shortName, const string& longName, const string& helpInfo,const T& defaultValue)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=true;
  m_help=helpInfo;
  m_defaultValue=defaultValue;
  m_hasDefault=true;
}

template<class T> void Optionpk<T>::setAll(const string& shortName, const string& longName, const string& helpInfo)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=true;
  m_help=helpInfo;
}

template<> void Optionpk<bool>::setAll(const string& shortName, const string& longName, const string& helpInfo,const bool& defaultValue);

template<> void Optionpk<bool>::setAll(const string& shortName, const string& longName, const string& helpInfo);

template<> Optionpk<bool>::Optionpk(const string& shortName, const string& longName, const string& helpInfo,const bool& defaultValue)
{
  setAll(shortName,longName,helpInfo,defaultValue);
}

template<> void Optionpk<bool>::setAll(const string& shortName, const string& longName, const string& helpInfo,const bool& defaultValue)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=false;
  m_help=helpInfo;
  m_defaultValue=defaultValue;
  m_hasDefault=true;
}

template<class T> Optionpk<T>::~Optionpk() 
{
}

template<class T> int Optionpk<T>::retrieveOption(int argc, char **argv){ 
  for(int i = 1; i < argc; ++i ){
    string helpStringShort="-h";
    string helpStringLong="--help";
    string currentArgument;
    string currentOption=argv[i];
    string shortOption=m_shortName;
    string longOption=m_longName;
    shortOption.insert(0,"-");
    longOption.insert(0,"--");
    // size_t foundEqual=currentOption.rfind("=");
    size_t foundEqual=currentOption.rfind("=");
    if(foundEqual!=string::npos){
      currentArgument=currentOption.substr(foundEqual+1);
      currentOption=currentOption.substr(0,foundEqual);
    }
    if(!(helpStringShort.compare(currentOption))||!(helpStringLong.compare(currentOption)))
      cout << usage() << endl;
    if(hasShortOption()&&!(shortOption.compare(currentOption))){//for -option
      if(foundEqual!=string::npos)
        this->push_back(string2type<T>(currentArgument));
      else if(m_hasArgument && i < argc-1)
        this->push_back(string2type<T>(argv[++i]));
      else
        this->push_back(string2type<T>("1"));
    }
    else if(hasLongOption()&&!(longOption.compare(currentOption))){//for --option
      if(foundEqual!=string::npos)
        this->push_back(string2type<T>(currentArgument));
      else if(m_hasArgument && i < argc-1)
        this->push_back(string2type<T>(argv[++i]));
      else
        this->push_back(string2type<T>("1"));
    }
  }
  if(!(this->size())&&m_hasDefault)//only set default value if no options were given
    this->push_back(m_defaultValue);
  return(this->size());
}

template<class T> std::vector<string>::const_iterator Optionpk<T>::findSubstring(string argument) const{
  std::vector<string>::const_iterator opit=this->begin();
  while(opit!=this->end()){
    if(opit->find(argument)!=std::string::npos)
      break;
      ++opit;
  }
  return opit;
}

#endif

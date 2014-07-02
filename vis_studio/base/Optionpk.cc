/**********************************************************************
Optionpk.cc: source file used for specialization of template class 
Optionpk defined in Optionpk.h
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

#include "base/Optionpk.h"

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

///specialization for bool
template<> inline void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo);

template<> inline void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo)
{
  m_shortName=shortName;
  m_longName=longName;
  m_hasArgument=false;
  m_help=helpInfo;
  m_hide=0;
}

///specialization for bool
template<> inline void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const bool& defaultValue, short hide);

///specialization for bool
template<> inline void Optionpk<bool>::setAll(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const bool& defaultValue, short hide)
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
template<> inline Optionpk<bool>::Optionpk(const std::string& shortName, const std::string& longName, const std::string& helpInfo,const bool& defaultValue, short hide)
{
  setAll(shortName,longName,helpInfo,defaultValue, hide);
}

//specialization (only makes sense for T=std::string), generic function throws exception
//find a substring in string option (e.g., option is of type -co INTERLEAVE=BAND)
template<> inline std::vector<std::string>::const_iterator Optionpk<std::string>::findSubstring(const std::string& argument) const{
  std::vector<std::string>::const_iterator opit=this->begin();
  while(opit!=this->end()){
    if(opit->find(argument)!=std::string::npos)
      break;
    ++opit;
  }
  return opit;
}

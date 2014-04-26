/**********************************************************************
Egcs.cc: Conversions from and to european grid coding system
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
#include "Egcs.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <assert.h>

Egcs::Egcs(){
}

// Egcs::Egcs(unsigned short level)
//   : m_level(level){
// }


Egcs::Egcs(unsigned short level)
{
  m_level=level;
}

unsigned short Egcs::cell2level(const std::string& cellCode) const{
  size_t pos=cellCode.find("-");
  std::string TILE=cellCode.substr(0,pos);
  unsigned short level=0;
  int base_level=19-(TILE.size()/2-1)*3;
  int quad_level=0;
  if(pos!=std::string::npos)
    quad_level=cellCode.size()-pos-1;
  if(quad_level>1)
    level=base_level-2;
  else if(quad_level>0)
    level=base_level-1;
  else
    level=base_level;
  return level;
}

Egcs::~Egcs(){
}

unsigned short Egcs::res2level(double resolution) const{
  double base=pow(10,log(resolution*4.0)/log(10.0));
  double diff=base/(2*resolution);
  return 0.5+(log(base)/log(10.0)*3+1-diff);
}

double Egcs::getResolution() const{
  unsigned short exponent=(m_level+1)/3;
  double base=pow(10.0,exponent);
  if((m_level)%3==2)
    return(base/4);
  else if((m_level)%3==0)
    return(base/2);
  else
    return(base);
}

void Egcs::force2grid(double& ulx, double& uly, double& lrx, double &lry) const{
  double dx=getResolution();
  double dy=dx;
  //ulx
  ulx=floor(ulx);
  ulx-=static_cast<int>(ulx)%(static_cast<int>(dx));
  //uly
  uly=ceil(uly);
  if(static_cast<int>(uly)%static_cast<int>(dy))
    uly+=dy;
  uly-=static_cast<int>(uly)%(static_cast<int>(dy));
  //lrx
  lrx=ceil(lrx);
  if(static_cast<int>(lrx)%static_cast<int>(dx))
    lrx+=dx;
  lrx-=static_cast<int>(lrx)%(static_cast<int>(dx));
  //lry
  lry=floor(lry);
  lry-=static_cast<int>(lry)%(static_cast<int>(dy));
}

void Egcs::cell2bb(const std::string& cellCode, int &ulx, int &uly, int &lrx, int &lry) const{
  size_t pos=cellCode.find("-");
  std::string TILE=cellCode.substr(0,pos);
  std::string TILEX=TILE.substr(0,TILE.size()/2);
  std::string TILEY=TILE.substr(TILE.size()/2);
  std::string QUAD;
  char QUAD1,QUAD2;
  std::istringstream stilex(TILEX);
  std::istringstream stiley(TILEY);
  int llx,lly;
  stilex >> llx;
  stiley >> lly;
  llx*=getBaseSize();
  lly*=getBaseSize();
  switch((19-m_level)%3){
  case(0)://there should be no QUAD
    assert(pos==std::string::npos);
    break;
  case(2)://there is a QUAD2
    assert(pos+1!=std::string::npos);
    QUAD=cellCode.substr(pos+1);
    QUAD2=QUAD.substr(1,1).c_str()[0];
    switch(QUAD2){
    case('A'):
      break;
    case('C'):
      lly+=getBaseSize()/4;
      break;
    case('D'):
      lly+=getBaseSize()/4;
    case('B'):
      llx+=getBaseSize()/4;
    break;
    }
  case(1)://QUAD1: deliberate fall through from case(2)!
    if(!QUAD.size()){
      assert(pos+1!=std::string::npos);
      QUAD=cellCode.substr(pos+1);
    }
    QUAD1=QUAD.substr(0,1).c_str()[0];
    switch(QUAD1){
    case('A'):
      break;
    case('C'):
      lly+=getBaseSize()/2;
      break;
    case('D'):
      lly+=getBaseSize()/2;
    case('B'):
      llx+=getBaseSize()/2;
      break;
    }
    break;
  }
  ulx=llx;
  uly=static_cast<int>(lly+getSize());
  lrx=static_cast<int>(llx+getSize());
  lry=lly;
}

void Egcs::cell2mid(const std::string& cellCode, double& midX, double& midY) const{
  int ulx,uly,lrx,lry;
  cell2bb(cellCode,ulx,uly,lrx,lry);
  midX=(ulx+lrx)/2.0;
  midY=(lry+uly)/2.0;
}

std::string Egcs::geo2cell(double geoX, double geoY) const{
  int ndgts=7-(m_level+1)/3;
  double xcel=static_cast<int>(geoX)/getBaseSize();
  double ycel=static_cast<int>(geoY)/getBaseSize();
  std::ostringstream osx;
  std::ostringstream osy, osxy;
  // osx << setprecision(ndgts) << xcel;
  // osy << setprecision(ndgts) << ycel;
  // osx << setprecision(0) << fixed << geoX;
  // osy << setprecision(0) << fixed << geoY;
  osx << std::fixed << geoX;
  osy << std::fixed << geoY;
  std::string quad1="";
  std::string quad2="";
  switch((19-m_level)%3){
  case(2):
    if(static_cast<int>(geoX)%(getBaseSize()/2)>=getBaseSize()/2.0){
      if(static_cast<int>(geoY)%(getBaseSize()/2)>=getBaseSize()/2.0)
        quad2="D";
      else
        quad2="B";
    }
    else if(static_cast<int>(geoY)%getBaseSize()>=getBaseSize()/4.0)
      quad2="C";
    else
      quad2="A";
  case(1)://deliberate fall through!
    if(static_cast<int>(geoX)%getBaseSize()>=getBaseSize()/2.0){
      if(static_cast<int>(geoY)%getBaseSize()>=getBaseSize()/2.0)
        quad1="-D";
      else
        quad1="-B";
    }
    else if(static_cast<int>(geoY)%getBaseSize()>=getBaseSize()/2.0)
      quad1="-C";
    else
      quad1="-A";
    break;
  default:
    break;
  }
  osxy << osx.str().substr(0,ndgts) << osy.str().substr(0,ndgts) << quad1 << quad2;
  return osxy.str();
}


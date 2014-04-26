/**********************************************************************
Egcs.h: Conversions from and to european grid coding system
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

#ifndef _EGCS_H_
#define _EGCS_H_

class Egcs
{
public:
  Egcs();
  Egcs(unsigned short level);
  /* Egcs(unsigned short level); */
  ~Egcs();
  unsigned short cell2level(const std::string& cellCode) const;
  std::string geo2cell(double x, double y) const;
  double getSize() const {return getBaseSize()*pow(2.0,(m_level-19)%3);};
  void setLevel(unsigned short level){m_level=level;};
  unsigned short getLevel() const{return m_level;};
  unsigned short res2level(double resolution) const;
  double getResolution() const;
  void force2grid(double& ulx, double& uly, double& lrx, double &lry) const;
  void cell2bb(const std::string& cellCode, int &ulx, int &uly, int &lrx, int &lry) const;
  void cell2mid(const std::string& cellCode, double& midX, double& midY) const;
private:
  int getBaseSize() const {return pow(10.0,(m_level+1)/3);};
  unsigned short m_level;
// level square scheme         example
// 19    1000km xy             32
// 18     500km xy-q           32-A
// 17     250km xy-qq          32-AB
// 16     100km xxyy           3320
// 5       25m  xxxxxyyyyy-qq  3346720658-DC
// 1        1m  xxxxxxxyyyyyyy 33467652065889
};
#endif // _EGCS_H_


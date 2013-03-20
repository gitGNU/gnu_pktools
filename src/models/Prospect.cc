/**********************************************************************
Prospect.cc: class for radiative transfer leaf PROSPECT
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
#include "Prospect.h"

void prospect::Prospect::getLeafSpectrum(std::vector<double>& leafReflectance, std::vector<double>& leafTransmittance) const{
  double RT[2][2101];//400 nm - 2500 nm
  prospect_5b_(m_N,m_Cab,m_Car,m_Cbrown,m_Cw,m_Cm,RT);
  leafReflectance.resize(2101);
  leafTransmittance.resize(2101);
  for(int index=0;index<2101;++index){
    leafReflectance[index]=RT[0][index];
    leafTransmittance[index]=RT[1][index];
  }
}

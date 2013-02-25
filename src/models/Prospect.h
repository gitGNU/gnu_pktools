/**********************************************************************
Prospect.h: class for radiative transfer leaf PROSPECT
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
#ifndef _PROSPECT_
#define _PROSPECT_
#include <assert.h>
#include <math.h>
#include <vector>

extern "C" {
  void prospect_5b__(const double& N,const double& Cab,const double& Car, const double& Cbrown, const double& Cw, const double& Cm, double RT[][2]);
}

namespace prospect
{
  class Prospect{
  public:
    Prospect(void){};
    ~Prospect(void){};
    void setN(double n){m_N=n;};
    void setCab(double cab){m_Cab=cab;};
    void setCar(double car){m_Car=car;};
    void setCbrown(double cbrown){m_Cbrown=cbrown;};
    void setCw(double cw){m_Cw=cw;};
    void setCm(double cm){m_Cm=cm;};
    double getN() const {return m_N;};
    double getCab() const {return m_Cab;};
    double getCar() const {return m_Car;};
    double getCbrown() const {return m_Cbrown;};
    double getCw() const {return m_Cw;};
    double getCm() const {return m_Cm;};
    void getLeafSpectrum(std::vector<double>& leafReflectance,std::vector<double>& leafTransmittance) const;
  private:
    double m_N;              // structure coefficient
    double m_Cab;            // chlorophyll content (g.cm-2) 
    double m_Car;            // carotenoid content (g.cm-2)
    double m_Cbrown;         // brown pigment content (arbitrary units)
    double m_Cw;             // EWT (cm)
    double m_Cm;             // LMA (g.cm-2)
  };
}
#endif //_PROSPECT_

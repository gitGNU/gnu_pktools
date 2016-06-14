/**********************************************************************
ImgRegression.h: class to calculate regression between two raster datasets
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
#ifndef _IMGREGRESSION_H_
#define _IMGREGRESSION_H_

#include <vector>
#include "imageclasses/ImgRaster.h"
#include "StatFactory.h"

namespace imgregression
{
  class ImgRegression{
  public:
    ImgRegression(void);
    ~ImgRegression(void);
    double getRMSE(ImgRaster& imgReader1, ImgRaster& imgReader2, double &c0, double &c1, unsigned int b1=0, unsigned int b2=0, short verbose=0) const;
    double getRMSE(ImgRaster& imgReader, unsigned int b1, unsigned int b2, double& c0, double& c1, short verbose=0) const;
    double getR2(ImgRaster& imgReader1, ImgRaster& imgReader2, double &c0, double &c1, unsigned int b1=0, unsigned int b2=0, short verbose=0) const;
    double pgetR2(ImgRaster& imgReader1, ImgRaster& imgReader2, double& c0, double& c1, unsigned int band1, unsigned int band2, short verbose=0) const;
    double getR2(ImgRaster& imgReader, unsigned int b1, unsigned int b2, double& c0, double& c1, short verbose=0) const;
    double pgetR2(ImgRaster& imgReader, unsigned int band1, unsigned int band2, double& c0, double& c1, short verbose=0) const;

    void setThreshold(double theThreshold){m_threshold=theThreshold;};
    void setDown(int theDown){m_down=theDown;};
  private:
    int m_down;
    double m_threshold;
  };
}
#endif //_IMGREGRESSION_H_

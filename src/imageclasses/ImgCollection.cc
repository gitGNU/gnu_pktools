/**********************************************************************
ImgCollection.cc: class to read raster files using GDAL API library
Copyright (C) 2008-2016 Pieter Kempeneers

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
#include <vector>
#include <string>
#include <iostream>
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"
#include "algorithms/Egcs.h"
#include "apps/AppFactory.h"
#include "ImgCollection.h"

using namespace std;
using namespace app;

//   ///time period covering the image collection (check http://www.boost.org/doc/libs/1_55_0/doc/html/date_time/examples.html#date_time.examples.time_periods for how to use boost period)
// boost::posix_time::time_period ImgCollection::getTimePeriod(){
//   if(m_time.size()){
//     std::vector<boost::posix_time::time_period>::iterator tit=m_time.begin();
//     boost::posix_time::time_period timePeriod=*(tit++);
//     while(tit!=m_time.end()){
//       timePeriod.span(*(tit++));
//     }
//     return(timePeriod);
//   }
// }    

/**
 * @param ulx upper left coordinate in x
 * @param uly upper left coordinate in y
 * @param lrx lower left coordinate in x
 * @param lry lower left coordinate in y
 **/
void ImgCollection::getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const{
  std::vector<std::shared_ptr<ImgRasterGdal> >::const_iterator it=begin();
  if(it!=end())
    (*(it++))->getBoundingBox(ulx,uly,lrx,lry);
  while(it!=end()){
    double imgulx,imguly,imglrx,imglry;
    (*(it++))->getBoundingBox(imgulx,imguly,imglrx,imglry);
    ulx=(ulx>imgulx)? imgulx : ulx;
    uly=(uly<imguly)? imguly : uly;
    lrx=(lrx<imglrx)? imglrx : lrx;
    lry=(lry>imglry)? imglry : lry;
  }
}    

/**
 * @param x,y georeferenced coordinates in x and y
 * @return true if image covers the georeferenced location
 **/
bool ImgCollection::covers(double x, double  y) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((x > theULX)&&
         (x < theLRX)&&
         (y < theULY)&&
         (y >theLRY));
}

/**
 * @param ulx upper left coordinate in x
 * @param uly upper left coordinate in y
 * @param lrx lower left coordinate in x
 * @param lry lower left coordinate in y
 * @return true if image (partially) covers the bounding box
 **/
bool ImgCollection::covers(double ulx, double  uly, double lrx, double lry) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((ulx < theLRX)&&(lrx > theULX)&&(lry < theULY)&&(uly > theLRY));
}

/**
 * @param noDataValues standard template library (stl) vector containing no data values
 * @return number of no data values in this dataset
 **/
unsigned int ImgCollection::getNoDataValues(std::vector<double>& noDataValues) const
{
  if(m_noDataValues.size()){
    noDataValues=m_noDataValues;
    return m_noDataValues.size();
  }
  else
    return 0;
}

/**
 * @param noDataValue no data value to be pushed for this dataset
 * @return number of no data values in this dataset
 **/
unsigned int ImgCollection::pushNoDataValue(double noDataValue)
{
  if(find(m_noDataValues.begin(),m_noDataValues.end(),noDataValue)==m_noDataValues.end())
    m_noDataValues.push_back(noDataValue);
  return(m_noDataValues.size());
}


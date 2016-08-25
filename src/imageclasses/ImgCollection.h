/**********************************************************************
ImgCollection.h: class to read raster files using GDAL API library
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
#ifndef _IMGCOLLECTION_H_
#define _IMGCOLLECTION_H_

#include <string>
#include <vector>
#include <memory>
// #include "boost/date_time/posix_time/posix_time.hpp"
#include "ImgReaderOgr.h"
#include "ImgRaster.h"
#include "apps/AppFactory.h"

namespace app{
class AppFactory;
}

/**
   This class is used to store a collection of raster images
**/
class ImgCollection : public std::vector<std::shared_ptr<ImgRaster> >
{
public:
  enum CRULE_TYPE {overwrite=0, maxndvi=1, maxband=2, minband=3, validband=4, mean=5, mode=6, median=7,sum=8,minallbands=9,maxallbands=10,stdev=11};
  ///default constructor
  ImgCollection(void) : m_index(0) {};// : std::vector<ImgRaster*>(), m_index(0) {};
  ImgCollection(unsigned int theSize){
    for(unsigned int iimg=0;iimg<theSize;++iimg){
      this->emplace_back(new(ImgRaster));
    }
    m_index=0;
  }
  ///destructor
  ~ImgCollection(void){};

  // ImgCollection(const ImgCollection&) = default;  
  // ImgCollection& operator=(const ImgCollection&) = default;

  ///get bounding box of image collection
  void getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const;
  // ///get begin and last time of image collection
  // boost::posix_time::time_period getTimePeriod();

  // ///filter collection according to period
  // void filterTime(const boost::posix_time::time_period& thePeriod){
  //   unsigned int index=0;
  //   std::vector<std::shared_ptr<ImgRaster>>::iterator it=begin();
  //   std::vector<boost::posix_time::time_period>::iterator tit=m_time.begin();
  //   while(it!=end()&&tit!=m_time.end()){
  //     if(thePeriod.contains(m_time[m_index])){
  //       ++it;
  //       ++tit;
  //     }
  //     else{
  //       it=erase(it);
  //       tit=m_time.erase(tit);
  //     }
  //   }
  //   m_index=0;
  // };
  ///filter collection according to bounding box
  void filterGeo(double ulx, double uly, double lrx, double lry){
    std::vector<std::shared_ptr<ImgRaster>>::iterator it=begin();
    // std::vector<boost::posix_time::time_period>::iterator tit=m_time.begin();
    while(it!=end()){
      if((*it)->covers(ulx, uly, lrx, lry)){
        ++it;
        // if(tit!=m_time.end())
        //   ++tit;
      }
      else{
        it=erase(it);
        // if(tit!=m_time.end())
        //   tit=m_time.erase(tit);
      }
    }
    m_index=0;
  };
  ///filter collection according to position
  void filterGeo(double x, double y){
    std::vector<std::shared_ptr<ImgRaster>>::iterator it=begin();
    // std::vector<boost::posix_time::time_period>::iterator tit=m_time.begin();
    while(it!=end()){
      if((*it)->covers(x,y)){
        ++it;
        // if(tit!=m_time.end())
        //   ++tit;
      }
      else{
        it=erase(it);
        // if(tit!=m_time.end())
        //   tit=m_time.erase(tit);
      }
    }
    m_index=0;
};
  ///push image to collection
  void pushImage(std::shared_ptr<ImgRaster> imgRaster){
    this->emplace_back(imgRaster);
  };
  // ///push image to collection with corresponding period
  // void pushImage(std::shared_ptr<ImgRaster> imgRaster, boost::posix_time::time_period imgPeriod){
  //   this->emplace_back(imgRaster);
  //   // m_time.push_back(imgPeriod);
  // };
  // ///push image period
  // void pushTime(boost::posix_time::time_period imgPeriod){
  //   m_time.push_back(imgPeriod);
  // };
  // ///set image periods for collection
  // void setTime(const std::vector<boost::posix_time::time_period>& timeVector){
  //   m_time=timeVector;
  // };
  // ///get image periods for collection
  // void getTime(std::vector<boost::posix_time::time_period>& timeVector){
  //   timeVector=m_time;
  // };
  ///Check if a geolocation is covered by this dataset. Only the bounding box is checked, irrespective of no data values.
  bool covers(double x, double y) const;
  ///Check if a region of interest is (partially) covered by this dataset. Only the bounding box is checked, irrespective of no data values.
  bool covers(double ulx, double  uly, double lrx, double lry) const;

  // std::shared_ptr<ImgRaster> getNextImage(){
  //   if(m_index<size())
  //     return(this->at(m_index++));
  //   else
  //     return(0);
  // }
  // std::shared_ptr<ImgRaster> getNextImage(boost::posix_time::time_period& imgPeriod){
  //   if(m_index<size()){
  //     if(m_index<m_time.size())
  //       imgPeriod=m_time[m_index];
  //     return(this->at(m_index++));
  //   }
  //   else
  //     return(0);
  // }
  void resetIterator(){m_index=0;};
  void clean(){clear();m_index=0;};
  void close(){
    for(std::vector<std::shared_ptr<ImgRaster>>::iterator it=begin();it!=end();++it)
      (*it)->close();
  }
  ///Get the no data values of this dataset as a standard template library (stl) vector
  unsigned int getNoDataValues(std::vector<double>& noDataValues) const;
  ///Check if value is nodata in this dataset
  bool isNoData(double value) const{if(m_noDataValues.empty()) return false;else return find(m_noDataValues.begin(),m_noDataValues.end(),value)!=m_noDataValues.end();};
  ///Push a no data value for this dataset
  unsigned int pushNoDataValue(double noDataValue);
  ///Set the no data values of this dataset using a standard template library (stl) vector as input
  unsigned int setNoData(const std::vector<double>& nodata){m_noDataValues=nodata; return(m_noDataValues.size());};
  ///composite image
  CPLErr composite(std::shared_ptr<ImgRaster> imgWriter, const app::AppFactory& app);
  ///composite image only for in memory
  std::shared_ptr<ImgRaster> composite(app::AppFactory& app);
  ///crop image
  CPLErr crop(std::shared_ptr<ImgRaster> imgWriter, const app::AppFactory& app);
  ///crop image only for in memory
  std::shared_ptr<ImgRaster> crop(app::AppFactory& app);
private:
  unsigned int m_index;
  std::vector<double> m_noDataValues;
  // std::vector<boost::posix_time::time_period> m_time;
};

#endif // _IMGCOLLECTION_H_

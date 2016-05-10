/**********************************************************************
AppFactory.h: class for application functions
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
#ifndef _APPFACTORY_H_
#define _APPFACTORY_H_

#include <algorithm>
#include <vector>
#include <iostream>
#include <string>
#include "imageclasses/ImgReaderGdal.h"
#include "imageclasses/ImgWriterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"

namespace appfactory
{

namespace crule{
  enum CRULE_TYPE {overwrite=0, maxndvi=1, maxband=2, minband=3, validband=4, mean=5, mode=6, median=7,sum=8,minallbands=9,maxallbands=10,stdev=11};
}
  
class AppFactory{

public:
 AppFactory(unsigned long int memory=-1) : m_memory(memory){};
  virtual ~AppFactory(void){};
  //todo: support arguments as list of arguments, class or struct and xml file?
  bool pkcomposite(std::vector<ImgReaderGdal>& input, ImgWriterGdal& imgWriter);   

private:
    //todo: create member attribute for pointer to memory buffer?
    unsigned long int m_memory;
};


}

#endif /* _APPFACTORY_H_ */

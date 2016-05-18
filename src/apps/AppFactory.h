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
#include "imageclasses/ImgRasterGdal.h"
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"
#include "base/Optionpk.h"
#include "algorithms/StatFactory.h"

namespace app
{

namespace crule{
  enum CRULE_TYPE {overwrite=0, maxndvi=1, maxband=2, minband=3, validband=4, mean=5, mode=6, median=7,sum=8,minallbands=9,maxallbands=10,stdev=11};
}
  
  class AppFactory{

  public:
    AppFactory(void) : m_argc(1), m_argv(std::vector<std::string>(1,"appFactory")){};
    virtual ~AppFactory(void){};
    void setOptions(int argc, char* argv[]);
    ///set bool option (used as flag)
    void setOption(const std::string &key)
    {
      std::ostringstream os;
      os << "--" << key;;
      m_argv.push_back(os.str().c_str());
      ++m_argc;
    };
    ///set key value option
    void setOption(const std::string &key, const std::string &value)
    {
      std::ostringstream os;
      os << "--" << key;
      m_argv.push_back(os.str());
      ++m_argc;
      m_argv.push_back(value);
      ++m_argc;
    };
    void getHelp() {setOption("help");};
    void clearOptions() {m_argc=1;m_argv.clear();m_argv.push_back("appFactory");};
    void showOptions() 
    {
      for(int iarg=1;iarg<m_argv.size();++iarg)
        std::cout << m_argv[iarg] << " ";
      std::cout << std::endl;
    };
    bool pkcrop(std::vector<ImgRasterGdal>& input, ImgRasterGdal& imgWriter);   
    bool pkcomposite(std::vector<ImgRasterGdal>& input, ImgRasterGdal& imgWriter);   

  private:
    //todo: create member attribute for pointer to memory buffer?
    int m_argc;
    std::vector<std::string> m_argv;
  };


}

#endif /* _APPFACTORY_H_ */

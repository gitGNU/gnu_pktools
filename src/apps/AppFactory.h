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

#include <string>
#include <vector>
#include <iostream>
#include "base/Optionpk.h"

namespace app
{

  class AppFactory{

  public:
    AppFactory(void) : m_argc(1), m_argv(std::vector<std::string>(1,"appFactory")){}
    AppFactory(int argc, char* argv[]) : m_argc(1), m_argv(std::vector<std::string>(1,"appFactory")){
      setOptions(argc, argv);
    };
    virtual ~AppFactory(void){};
    bool empty() const {return(m_argv.empty());};
    bool size() const {return(m_argv.size());};
    void setOptions(int argc, const std::vector<std::string> argv){
      m_argc=argc;
      m_argv.clear();
      m_argv=argv;
    }
    void setOptions(int argc, char* argv[]){
      m_argc=argc;
      m_argv.clear();
      for(int iarg=0;iarg<argc;++iarg)
        m_argv.push_back(argv[iarg]);
    }
    ///push bool option (used as flag)
    void pushOption(const std::string &key)
    {
      std::ostringstream os;
      // os << "--" << key;;
      if(key=="help")
        os << "--" << key;
      else
        os << "-" << key;
      m_argv.push_back(os.str().c_str());
      ++m_argc;
    };
    ///set bool option (used as flag)
    void setOption(const std::string &key)
    {
      clearOption(key);
      std::ostringstream os;
      // os << "--" << key;;
      if(key=="help")
        os << "--" << key;
      else
        os << "-" << key;
      m_argv.push_back(os.str().c_str());
      ++m_argc;
    };
    ///template set key value option
    template<typename T> void setOption(const std::string &key, T value){
      setOption(key,type2string<T>(value));
    }
    ///set key value option
    void setOption(const std::string &key, const std::string &value)
    {
      clearOption(key);
      std::ostringstream os;
      os << "-" << key;
      m_argv.push_back(os.str());
      ++m_argc;
      m_argv.push_back(value);
      ++m_argc;
    };
    ///push key value option
    void pushOption(const std::string &key, const std::string &value)
    {
      std::ostringstream os;
      os << "-" << key;
      m_argv.push_back(os.str());
      ++m_argc;
      m_argv.push_back(value);
      ++m_argc;
    };
    void getHelp() {setOption("help");};
    void clearOptions() {m_argc=1;m_argv.clear();m_argv.push_back("appFactory");};
    void clearOption(const std::string &key)
    {
      std::vector<std::string>::iterator opit=m_argv.begin();
      while(opit!=m_argv.end()){
        if(opit->find("-"+key)!=std::string::npos){
          m_argv.erase(opit);
          --m_argc;
          if(opit!=m_argv.end()){
            if(opit->find("-")==std::string::npos){//not a bool option
              m_argv.erase(opit);
              --m_argc;
            }
          }
        }
        else
          ++opit;
      }
    };
    void showOptions() const
    {
      for(int iarg=1;iarg<m_argv.size();++iarg)
        std::cout << m_argv[iarg] << " ";
      std::cout << std::endl;
    };
    int getArgc() const {return m_argc;};
    std::string getArgv(unsigned int i) const {
      if((i>0)&&(i<m_argv.size()))
        return m_argv[i];
      else
        throw(std::string("Error: invalid index"));
    }
    std::vector<std::string> getArgv() const {
      return m_argv;
    }
    std::vector<std::string> getArgv() {
      return m_argv;
    }
  private:
    //todo: create member attribute for pointer to memory buffer?
    int m_argc;
    std::vector<std::string> m_argv;
  };


}

#endif /* _APPFACTORY_H_ */

/**********************************************************************
FileReaderAscii.h: class to read ASCII files using (colum based)
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
#ifndef _IMGREADERASCII_H_
#define _IMGREADERASCII_H_

#include <string>
#include <vector>
#include <fstream>
#include "base/Optionpk.h"
//#include <armadillo>

//--------------------------------------------------------------------------
class FileReaderAscii
{
public:
  FileReaderAscii();
  FileReaderAscii(const std::string& filename);
  FileReaderAscii(const std::string& filename, const char& fieldseparator);
  ~FileReaderAscii();
  void reset(){m_ifstream.clear();m_ifstream.seekg(0,std::ios::beg);};
  void open(const std::string& filename);
  void close(void);
  void setFieldSeparator(const char& fieldseparator){m_fs=fieldseparator;};
  void setMinRow(int minRow){m_minRow=minRow;};
  void setMaxRow(int maxRow){m_maxRow=maxRow;};
  void setComment(char comment){m_comment=comment;};
  unsigned int nrOfCol(bool checkCols=false, bool verbose=false);
  unsigned int nrOfRow(bool checkCols=false, bool verbose=false);
  template<class T> unsigned int readData(std::vector<std::vector<T> > &dataVector, const std::vector<int> &cols, double scale=1.0, double offset=0.0, bool transpose=false, bool verbose=false);
  template<class T> unsigned int readData(std::vector<T> &dataVector, int col, double scale=1.0, double offset=0, bool verbose=false);

  protected:
  std::string m_filename;
  std::ifstream m_ifstream;
  char m_fs;
  char m_comment;
  double m_min;
  double m_max;
  int m_minRow;
  int m_maxRow;
};

template<class T> unsigned int FileReaderAscii::readData(std::vector<T> &dataVector, int col, double scale, double offset, bool verbose){
  reset();
  dataVector.clear();
  int nrow=0;
  bool withinRange=true;
  if(m_fs>' '&&m_fs<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
    if(verbose)
      std::cout << "reading csv file " << m_filename << std::endl;
    std::string csvRecord;
    while(getline(m_ifstream,csvRecord)){//read a line
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        std::istringstream csvstream(csvRecord);
        std::string item;
        int ncol=0;
        bool isComment=false;
        while(getline(csvstream,item,m_fs)){//read a column
          if(verbose)
            std::cout << item << " ";
          size_t pos=item.find(m_comment);
          if(pos!=std::string::npos){
            isComment=true;
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
          }
          if(ncol==col){
            T value=scale*string2type<T>(item)+offset;
            if((value>=m_min&&value<=m_max)||m_max<=m_min)
              dataVector.push_back(value);
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose)
          std::cout << std::endl;
        if(dataVector.size()&&ncol<=col){
          std::ostringstream ess;
          ess << "Error: different number of cols found in line " << nrow << " (" << ncol << ")" << std::endl;
          throw(ess.str());
        }
      }
      ++nrow;
    }
    assert(dataVector.size());
  }
  else{//space or tab delimited fields
    if(verbose)
      std::cout << "space or tab delimited fields" << std::endl;
    std::string spaceRecord;
    while(!getline(m_ifstream, spaceRecord).eof()){
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        if(verbose>1)
          std::cout << spaceRecord << std::endl;
        std::istringstream lineStream(spaceRecord);
        std::string item;
        int ncol=0;
        bool isComment=false;
        while(lineStream >> item){
          if(verbose)
            std::cout << item << " ";
          // std::istringstream itemStream(item);
          size_t pos=item.find(m_comment);
          if(pos!=std::string::npos){
            isComment=true;
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
          }
          T value=scale*string2type<T>(item)+offset;
          // T value=string2type<T>(item);
          if(ncol==col){
            if((value>=m_min&&value<=m_max)||m_max<=m_min)
              dataVector.push_back(value);
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose>1)
          std::cout << std::endl;
        if(verbose)
          std::cout << "number of columns: " << ncol << std::endl;
        if(dataVector.size()&&ncol<=col){
          std::ostringstream ess;
          ess << "Error: different number of cols found in line " << nrow << " (" << ncol << ")" << std::endl;
          throw(ess.str());
        }
      }
      ++nrow;
    }
  }
  return dataVector.size();
}

template<class T> unsigned int FileReaderAscii::readData(std::vector<std::vector<T> > &dataVector, const std::vector<int> &cols, double scale, double offset, bool transpose, bool verbose){
  reset();
  dataVector.clear();
  if(!transpose)
    dataVector.resize(cols.size());
  int nrow=0;
  bool withinRange=true;
  if(m_fs>' '&&m_fs<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
    if(verbose)
      std::cout << "reading csv file " << m_filename << std::endl;
    std::string csvRecord;
    while(getline(m_ifstream,csvRecord)){//read a line
      std::vector<T> sampleVector;
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        std::istringstream csvstream(csvRecord);
        std::string item;
        int ncol=0;
        bool isComment=false;
        while(getline(csvstream,item,m_fs)){//read a column
          if(verbose)
            std::cout << item << " ";
          size_t pos=item.find(m_comment);
          if(pos!=std::string::npos){
            isComment=true;
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
          }
          for(int icol=0;icol<cols.size();++icol){
            if(ncol==cols[icol]){
              T value=scale*string2type<T>(item)+offset;
              // T value=string2type<T>(item);
              if((value>=m_min&&value<=m_max)||m_max<=m_min){
                if(transpose)
                  sampleVector.push_back(value);
                else
                  dataVector[icol].push_back(value);
              }
            }
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose)
          std::cout << std::endl;
        // if(dataVector.back().size())
        //   assert(ncol>=cols[0]);
      }
      if(sampleVector.size()&&transpose)
        dataVector.push_back(sampleVector);
      ++nrow;
    }
    assert(dataVector.size());
  }
  else{//space or tab delimited fields
    if(verbose)
      std::cout << "space or tab delimited fields" << std::endl;
    std::string spaceRecord;
    while(!getline(m_ifstream, spaceRecord).eof()){
      std::vector<T> sampleVector;
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        if(verbose>1)
          std::cout << spaceRecord << std::endl;
        std::istringstream lineStream(spaceRecord);
        std::string item;
        int ncol=0;
        bool isComment=false;
        while(lineStream >> item){
          if(verbose)
            std::cout << item << " ";
          // std::istringstream itemStream(item);
          size_t pos=item.find(m_comment);
          if(pos!=std::string::npos){
            isComment=true;
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
          }
          T value=scale*string2type<T>(item)+offset;
          // T value=string2type<T>(item);
          for(int icol=0;icol<cols.size();++icol){
            if(ncol==cols[icol]){
              if((value>=m_min&&value<=m_max)||m_max<=m_min){
                if(transpose)
                  sampleVector.push_back(value);
                else
                  dataVector[icol].push_back(value);
              }
            }
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose>1)
          std::cout << std::endl;
        if(verbose)
          std::cout << "number of columns: " << ncol << std::endl;
        // if(dataVector.back().size())
        //   assert(ncol>=cols[0]);
      }
      if(sampleVector.size()&&transpose)
        dataVector.push_back(sampleVector);
      ++nrow;
    }
  }
  return dataVector.size();
}

#endif // _IMGREADERASCII_H_

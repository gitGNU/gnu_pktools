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

//--------------------------------------------------------------------------
class FileReaderAscii
{
public:
  FileReaderAscii(void);
  FileReaderAscii(const std::string& filename);
  FileReaderAscii(const std::string& filename, const char& fieldseparator);
  ~FileReaderAscii(void);
  void reset(){m_ifstream.clear();m_ifstream.seekg(0,ios::beg);};
  void open(const std::string& filename);
  void close(void);
  void setFieldSeparator(const char& fieldseparator){m_fs=fieldseparator;};
  void setMinRow(int minRow){m_minRow=minRow;};
  void setMaxRow(int maxRow){m_maxRow=maxRow;};
  void setComment(char comment){m_comment=comment;};
  template<class T> unsigned int readData(vector<vector<T> > &dataVector, const vector<int> &cols);
  template<class T> unsigned int readData(vector<T> &dataVector, int col);
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

FileReaderAscii::FileReaderAscii(void)
  : m_min(0),m_max(0),m_minRow(0),m_maxRow(0),m_fs(' '),m_comment('#'){
}

FileReaderAscii::FileReaderAscii(const std::string& filename)
  : m_min(0),m_max(0),m_minRow(0),m_maxRow(0),m_fs(' '),m_comment('#'){
  open(filename);
}

FileReaderAscii::FileReaderAscii(const std::string& filename, const char& fieldseparator)
  : m_min(0),m_max(0),m_minRow(0),m_maxRow(0),m_fs(' '),m_comment(fieldseparator){
  open(filename);
}

FileReaderAscii::~FileReaderAscii(void)
{
}

void FileReaderAscii::open(const std::string& filename){
  m_filename=filename;
  m_ifstream.open(filename.c_str(),ios_base::in);
  if(!(m_ifstream)){
    string errorString;
    errorString="Error: could not open file ";
    errorString+=filename;
    throw(errorString);
  }
}

void FileReaderAscii::close(){
  m_ifstream.close();
  //  m_ifstream.clear();
}

template<class T> unsigned int FileReaderAscii::readData(vector<T> &dataVector, int col){
  reset();
  bool verbose=false;
  dataVector.clear();
  int nrow=0;
  bool withinRange=true;
  if(m_fs>' '&&m_fs<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
    if(verbose)
      cout << "reading csv file " << m_filename << endl;
    string csvRecord;
    while(getline(m_ifstream,csvRecord)){//read a line
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        istringstream csvstream(csvRecord);
        string item;
        int ncol=0;
        bool isComment=false;
        while(getline(csvstream,item,m_fs)){//read a column
          if(verbose)
            cout << item << " ";
          unsigned pos=item.find(m_comment);
          if(pos!=std::string::npos){
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
            isComment=true;
          }
          if(ncol==col){
            T value=string2type<T>(item);
            if((value>=m_min&&value<=m_max)||m_max<=m_min)
              dataVector.push_back(value);
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose)
          cout << endl;
        if(dataVector.size())
          assert(ncol>=col);
      }
      ++nrow;
    }
    assert(dataVector.size());
  }
  else{//space or tab delimited fields
    if(verbose)
      std::cout << "space or tab delimited fields" << std::endl;
    string spaceRecord;
    while(!getline(m_ifstream, spaceRecord).eof()){
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        if(verbose>1)
          cout << spaceRecord << endl;
        istringstream lineStream(spaceRecord);
        string item;
        int ncol=0;
        bool isComment=false;
        while(lineStream >> item){
          if(verbose)
            cout << item << " ";
          // istringstream itemStream(item);
          unsigned pos=item.find(m_comment);
          if(pos!=std::string::npos){
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
            isComment=true;
          }
          T value=string2type<T>(item);
          if(ncol==col){
            if((value>=m_min&&value<=m_max)||m_max<=m_min)
              dataVector.push_back(value);
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose>1)
          cout << endl;
        if(verbose)
          cout << "number of columns: " << ncol << endl;
        if(dataVector.size())
          assert(ncol>=col);
      }
      ++nrow;
    }
  }
  return dataVector.size();
}

template<class T> unsigned int FileReaderAscii::readData(vector<vector<T> > &dataVector, const vector<int> &cols){
  reset();
  bool verbose=false;
  dataVector.clear();
  dataVector.resize(cols.size());
  int nrow=0;
  bool withinRange=true;
  if(m_fs>' '&&m_fs<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
    if(verbose)
      cout << "reading csv file " << m_filename << endl;
    string csvRecord;
    while(getline(m_ifstream,csvRecord)){//read a line
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        istringstream csvstream(csvRecord);
        string item;
        int ncol=0;
        bool isComment=false;
        while(getline(csvstream,item,m_fs)){//read a column
          if(verbose)
            cout << item << " ";
          unsigned pos=item.find(m_comment);
          if(pos!=std::string::npos){
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
            isComment=true;
          }
          for(int icol=0;icol<cols.size();++icol){
            if(ncol==cols[icol]){
              T value=string2type<T>(item);
              if((value>=m_min&&value<=m_max)||m_max<=m_min)
                dataVector[icol].push_back(value);
            }
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose)
          cout << endl;
        if(dataVector.back().size())
          assert(ncol>=cols[0]);
      }
      ++nrow;
    }
    assert(dataVector.size());
  }
  else{//space or tab delimited fields
    if(verbose)
      std::cout << "space or tab delimited fields" << std::endl;
    string spaceRecord;
    while(!getline(m_ifstream, spaceRecord).eof()){
      withinRange=true;
      if(nrow<m_minRow)
        withinRange=false;
      if(m_maxRow>m_minRow)
        if(nrow>m_maxRow)
          withinRange=false;
      if(withinRange){
        if(verbose>1)
          cout << spaceRecord << endl;
        istringstream lineStream(spaceRecord);
        string item;
        int ncol=0;
        bool isComment=false;
        while(lineStream >> item){
          if(verbose)
            cout << item << " ";
          // istringstream itemStream(item);
          unsigned pos=item.find(m_comment);
          if(pos!=std::string::npos){
            if(pos>0)
              item=item.substr(0,pos-1);
            else
              break;
            if(verbose)
              std::cout << "comment found, string is " << item << std::endl;
            isComment=true;
          }
          T value=string2type<T>(item);
          for(int icol=0;icol<cols.size();++icol){
            if(ncol==cols[icol]){
              if((value>=m_min&&value<=m_max)||m_max<=m_min)
                dataVector[icol].push_back(value);
            }
          }
          ++ncol;
          if(isComment)
            break;
        }
        if(verbose>1)
          cout << endl;
        if(verbose)
          cout << "number of columns: " << ncol << endl;
        if(dataVector.back().size())
          assert(ncol>=cols[0]);
      }
      ++nrow;
    }
  }
  return dataVector.size();
}
#endif // _IMGREADERASCII_H_

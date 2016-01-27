/**********************************************************************
FileReaderLas.h: class to read LAS files using liblas API library
Copyright (C) 2008-2012 Pieter Kempeneers

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
#ifndef _IMGREADERLAS_H_
#define _IMGREADERLAS_H_

#include <string>
#include <vector>
#include "liblas/liblas.hpp"

//--------------------------------------------------------------------------
class LastReturnFilter: public liblas::FilterI
{
public:
  LastReturnFilter();
  bool filter(const liblas::Point& point);

private:
  LastReturnFilter(LastReturnFilter const& other);
  LastReturnFilter& operator=(LastReturnFilter const& rhs);
};

class FileReaderLas
{
public:
  FileReaderLas(void);
  FileReaderLas(const std::string& filename);
  ~FileReaderLas(void);
  void open(const std::string& filename);
  void close(void);
  liblas::Header const& getHeader() const;
  bool isCompressed() const;
  unsigned long int getPointCount() const;
  void las2ascii(const std::string& filename, bool verbose=false) const;
  template<typename T> liblas::Bounds<T> getExtent() const {return getHeader().GetExtent();};
  void getExtent(double& ulx, double& uly, double& lrx, double& lry) const;
  double getMinZ() const;
  double getMaxZ() const;
  liblas::Reader* getReader(){return m_reader;};
  void resetReader(){m_reader->Reset();};
  void setFilter(std::vector<liblas::FilterPtr> const& filters);
  bool const& readNextPoint(){return(m_reader->ReadNextPoint());};
  bool const& readNextPoint(liblas::Point& thePoint);
  liblas::Point const& getPoint(){return m_reader->GetPoint();};
  liblas::Point const& readPointAt(std::size_t n){m_reader->ReadPointAt(n);return m_reader->GetPoint();};
  // void addBoundsFilter(double ulx, double uly, double lrx, double lry);
  void addReturnsFilter(std::vector<unsigned short> const& returns);
  void addClassFilter(std::vector<unsigned short> const& classes);
  void setFilters(const std::vector<liblas::FilterPtr>& filters){m_filters=filters;setFilters();};
  void setFilters(){m_reader->SetFilters(m_filters);};
protected:
  void setCodec(const std::string& filename);
  std::string m_filename;
  std::ifstream *m_ifstream;
  liblas::Reader* m_reader;
  std::vector<liblas::FilterPtr> m_filters;
};

#endif // _IMGREADERLAS_H_

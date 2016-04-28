/**********************************************************************
ImgUpdaterGdal.h: class to read raster files using GDAL API library
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
#ifndef _IMGUPDATERGDAL_H_
#define _IMGUPDATERGDAL_H_

#include "ImgReaderGdal.h"
#include "ImgWriterGdal.h"

//--------------------------------------------------------------------------
class ImgUpdaterGdal : public ImgReaderGdal, public ImgWriterGdal
{
public:
  ///default constructor. Image needs to be opened later with one of the open methods.
  ImgUpdaterGdal(void);
  ///constructor opening an image in update mode (read and write). Caching is supported when memory>0.
  ImgUpdaterGdal(const std::string& filename, const GDALAccess& readMode=GA_Update, unsigned int memory=0);
  ///destructor
  ~ImgUpdaterGdal(void);
  ///Open an image in update mode (read and write). Caching is supported when memory>0.
  void open(const std::string& filename, const GDALAccess& readMode=GA_Update, unsigned int memory=0);
  ///close the image
  void close(void);

protected:
};

#endif // _IMGUPDATERGDAL_H_

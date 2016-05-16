/**********************************************************************
ImgUpdaterGdal.cc: class to read raster files using GDAL API library
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
#include <iostream>
#include "ImgUpdaterGdal.h"

ImgUpdaterGdal::ImgUpdaterGdal(void){};

/**
 * @param filename Open a raster dataset with this filename
 * @param readMode Open dataset in ReadOnly or Update mode
 * @param memory Available memory to cache image raster data (in MB)
 **/
ImgUpdaterGdal::ImgUpdaterGdal(const std::string& filename, const GDALAccess& readMode, unsigned int memory){
  open(filename,readMode,memory);
}

ImgUpdaterGdal::~ImgUpdaterGdal(void)
{
  // delete m_gds;
//   GDALDumpOpenDatasets(stderr);
//   GDALDestroyDriverManager();//could be used by other objects...
}

//--------------------------------------------------------------------------

/**
 * @param filename Open a raster dataset with this filename
 * @param readMode Open dataset in ReadOnly or Update mode
 * @param memory Available memory to cache image raster data (in MB)
 **/
void ImgUpdaterGdal::open(const std::string& filename, const GDALAccess& readMode, unsigned int memory)
{
  ImgReaderGdal::open(filename,readMode,memory);
}

void ImgUpdaterGdal::close(void)
{
  ImgRasterGdal::close();
}

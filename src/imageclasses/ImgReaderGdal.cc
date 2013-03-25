/**********************************************************************
ImgReaderGdal.cc: class to read raster files using GDAL API library
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
#include "ImgReaderGdal.h"
#include <assert.h>
#include <sstream>
#include <iostream>

ImgReaderGdal::ImgReaderGdal(void)
  : m_gds(NULL), m_isGeoRef(false), m_ncol(0), m_nrow(0), m_nband(0)
{}

void ImgReaderGdal::open(const string& filename)//, double magicX, double magicY)
{
  m_filename = filename;
  setCodec();//magicX,magicY);
}

ImgReaderGdal::~ImgReaderGdal(void)
{
  // delete m_gds;
//   GDALDumpOpenDatasets(stderr);
//   GDALDestroyDriverManager();//could be used by other objects...
}

//--------------------------------------------------------------------------
void ImgReaderGdal::close(void)
{
  GDALClose(m_gds);
}

void ImgReaderGdal::setCodec()//double magicX, double magicY)
{
  GDALAllRegister();
  m_gds = (GDALDataset *) GDALOpen(m_filename.c_str(), GA_ReadOnly );
  if(m_gds == NULL){
    string errorString="FileOpenError";
    throw(errorString);
  }
  m_ncol= m_gds->GetRasterXSize();
  m_nrow= m_gds->GetRasterYSize();
  m_nband= m_gds->GetRasterCount();
  m_isGeoRef=( static_cast<string>(m_gds->GetProjectionRef())  != "" );
  // m_magic_x=magicX;
  // m_magic_y=magicY;
  if(m_isGeoRef){
    double adfGeoTransform[6];
    if( m_gds->GetGeoTransform( adfGeoTransform ) == CE_None )
    {
      m_ulx=adfGeoTransform[0];
      m_uly=adfGeoTransform[3];
      m_delta_x=adfGeoTransform[1];
      m_delta_y=-adfGeoTransform[5];
    }
  }
  else{
    m_ulx=0;
    m_uly=nrOfRow();
    m_delta_x=1;
    m_delta_y=1;
  }
}

string ImgReaderGdal::getProjection(void) const 
{
  string theProjection=m_gds->GetProjectionRef();
  // size_t startpos,endpos;
  // while((startpos=theProjection.find(",AUTHORITY"))!=string::npos){
  //   endpos=theProjection.find("]",startpos+1,1)+1;
  //   theProjection.erase(startpos,endpos-startpos);
  // }
  return theProjection;
}

string ImgReaderGdal::getProjectionRef(void) const 
{
  string theProjection;
  if(m_gds->GetProjectionRef())
    return(m_gds->GetProjectionRef());
  else
    return "";
}

GDALDataType ImgReaderGdal::getDataType(int band) const
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1))->GetRasterDataType();
}

GDALRasterBand* ImgReaderGdal::getRasterBand(int band)
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1));
}

GDALColorTable* ImgReaderGdal::getColorTable(int band) const
{
  assert(band<m_nband+1);
  return (m_gds->GetRasterBand(band+1))->GetColorTable();
}

string ImgReaderGdal::getDriverDescription() const
{
  return m_gds->GetDriver()->GetDescription();
}

void ImgReaderGdal::getGeoTransform(double& ulx, double& uly, double& deltaX, double& deltaY, double& rot1, double& rot2) const
{
  double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
  m_gds->GetGeoTransform(adfGeoTransform);
  ulx=adfGeoTransform[0];
  deltaX=adfGeoTransform[1];
  rot1=adfGeoTransform[2];
  uly=adfGeoTransform[3];
  rot2=adfGeoTransform[4];
  deltaY=-adfGeoTransform[5];//convention of GDAL!
}

string ImgReaderGdal::getGeoTransform() const
{
  if(!isGeoRef())
    return("");
  else{
    double adfGeoTransform[6];// { 444720, 30, 0, 3751320, 0, -30 };
    m_gds->GetGeoTransform(adfGeoTransform);
    double ulx=adfGeoTransform[0];
    double deltaX=adfGeoTransform[1];
    double rot1=adfGeoTransform[2];
    double uly=adfGeoTransform[3];
    double rot2=adfGeoTransform[4];
    double deltaY=-adfGeoTransform[5];//convention of GDAL!
    ostringstream s;
    s << "[" << ulx << "," << deltaX << "," << rot1 << "," << uly << "," << rot2 << "," << -deltaY << "]";
    return(s.str());
  }
}

char** ImgReaderGdal::getMetadata()
{
  if(m_gds->GetMetadata()!=NULL)
    return(m_gds->GetMetadata());
  else
    return (char**)"";
}

char** ImgReaderGdal::getMetadata() const
{
  if(m_gds->GetMetadata()!=NULL)
    return(m_gds->GetMetadata());
  else 
    return (char**)"";
}

void ImgReaderGdal::getMetadata(list<string>& metadata) const
{
  char** cmetadata=m_gds->GetMetadata();
  while(*cmetadata!=NULL){
    metadata.push_back(*(cmetadata));
    ++cmetadata;
  }
}

string ImgReaderGdal::getDescription() const
{
  if(m_gds->GetDriver()->GetDescription()!=NULL)
    return m_gds->GetDriver()->GetDescription();
  else
    return "";
}

string ImgReaderGdal::getMetadataItem() const 
{
  if(m_gds->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME )!=NULL)
    return m_gds->GetDriver()->GetMetadataItem( GDAL_DMD_LONGNAME );
  else
    return "";
}
string ImgReaderGdal::getImageDescription() const 
{
  if(m_gds->GetDriver()->GetMetadataItem("TIFFTAG_IMAGEDESCRIPTION")!=NULL)
    return m_gds->GetDriver()->GetMetadataItem("TIFFTAG_IMAGEDESCRIPTION");
  else
    return "";
}

string ImgReaderGdal::getInterleave() const
{
  if(m_gds->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE"))
    return m_gds->GetMetadataItem( "INTERLEAVE", "IMAGE_STRUCTURE");
  else
    return("BAND");
}

string ImgReaderGdal::getCompression() const
{
  if(m_gds->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE"))
    return m_gds->GetMetadataItem( "COMPRESSION", "IMAGE_STRUCTURE");
  else
    return("NONE");
}

bool ImgReaderGdal::getBoundingBox(double& ulx, double& uly, double& lrx, double& lry) const
{
  if(m_isGeoRef){
    //    ulx=m_ulx-(m_magic_x-1.0)*m_delta_x;
    //    uly=m_uly+(m_magic_y-1.0)*m_delta_y;
    ulx=m_ulx;
    uly=m_uly;
    lrx=ulx+nrOfCol()*m_delta_x;
    lry=uly-nrOfRow()*m_delta_y;
    return true;
  }
  else{
    ulx=0;
    uly=nrOfRow()-1;
    lrx=nrOfCol()-1;
    lry=0;
    return false;
  }
}

bool ImgReaderGdal::getCentrePos(double& x, double& y) const
{
  if(m_isGeoRef){
//     x=m_ulx+(nrOfCol()/2.0-(m_magic_x-1.0))*m_delta_x;
//     y=m_uly-(nrOfRow()/2.0-(m_magic_y-1.0))*m_delta_y;
    x=m_ulx+(nrOfCol()/2.0)*m_delta_x;
    y=m_uly-(nrOfRow()/2.0)*m_delta_y;
    return true;
  }
  else{
    x=nrOfCol()/2.0;
    y=nrOfRow()/2.0;
    return false;
  }
}

//i and j represent fraction of pixels, return true if image is georeferenced
bool ImgReaderGdal::geo2image(double x, double y, double& i, double& j) const
{
  //double values are returned, caller is responsible for interpolation step
  if(m_isGeoRef){
//     double ulx=m_ulx-(m_magic_x-1.0)*m_delta_x;
//     double uly=m_uly+(m_magic_y-1.0)*m_delta_y;
    double ulx=m_ulx;
    double uly=m_uly;
    i=(x-ulx)/m_delta_x;
    j=(uly-y)/m_delta_y;
    return true;
  }
  else{
    i=x;
    j=nrOfRow()-y;
    return false;
  }
}

//x and y represent centre of pixel, return true if image is georeferenced
bool ImgReaderGdal::image2geo(double i, double j, double& x, double& y) const
{
  if(m_isGeoRef){
//     x=m_ulx+(1.5-m_magic_x+i)*m_delta_x;
//     y=m_uly-(1.5-m_magic_y+j)*m_delta_y;
    x=m_ulx+(0.5+i)*m_delta_x;
    y=m_uly-(0.5+j)*m_delta_y;
    return true;
  }
  else{
    x=0.5+i;
    y=nrOfRow()-(0.5+j);
    return false;
  }
}

bool ImgReaderGdal::covers(double x, double  y) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((x > theULX)&&
         (x < theLRX)&&
         (y < theULY)&&
         (y >theLRY));
}

bool ImgReaderGdal::covers(double ulx, double  uly, double lrx, double lry) const
{
  double theULX, theULY, theLRX, theLRY;
  getBoundingBox(theULX,theULY,theLRX,theLRY);
  return((ulx < theLRX)&&(lrx > theULX)&&(lry < theULY)&&(uly > theLRY));
}

double ImgReaderGdal::getMin(int& x, int& y, int band) const{
  double minValue=0;
  vector<double> lineBuffer(nrOfCol());
  bool init=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,GDT_Float64,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
      if(valid){
        if(!init){
          y=irow;
          x=icol;
          minValue=lineBuffer[icol];
          init=true;
        }
        else if(minValue>lineBuffer[icol]){
          y=irow;
          x=icol;
          minValue=lineBuffer[icol];
        }
      }
    }
  }
  if(init)
    return minValue;
  else
    throw(static_cast<string>("Warning: not initialized"));
}

double ImgReaderGdal::getMax(int& x, int& y, int band) const{
  double maxValue=0;
  vector<double> lineBuffer(nrOfCol());
  bool init=false;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,GDT_Float64,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
      if(valid){
        if(!init){
          y=irow;
          x=icol;
          maxValue=lineBuffer[icol];
          init=true;
        }
        else if(maxValue<lineBuffer[icol]){
          y=irow;
          x=icol;
          maxValue=lineBuffer[icol];
        }
      }
    }
  }
  if(init)
    return maxValue;
  else
    throw(static_cast<string>("Warning: not initialized"));
}

void ImgReaderGdal::getMinMax(int startCol, int endCol, int startRow, int endRow, int band, double& minValue, double& maxValue) const
{
  GDALRasterBand  *poBand;
  int             bGotMin, bGotMax;
  double          adfMinMax[2];
        
  poBand = m_gds->GetRasterBand(band+1);
  adfMinMax[0] = poBand->GetMinimum( &bGotMin );
  adfMinMax[1] = poBand->GetMaximum( &bGotMax );
  if( ! (bGotMin && bGotMax) )
    GDALComputeRasterMinMax((GDALRasterBandH)poBand, FALSE, adfMinMax);
    // GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
  minValue=adfMinMax[0];
  maxValue=adfMinMax[1];

  vector<double> lineBuffer(endCol-startCol+1);
  bool init=false;
  assert(endRow<nrOfRow());
  for(int irow=startCol;irow<endRow+1;++irow){
    readData(lineBuffer,GDT_Float64,startCol,endCol,irow,band);
    for(int icol=0;icol<lineBuffer.size();++icol){
      bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
      if(valid){
	if(!init){
	  minValue=lineBuffer[icol];
	  maxValue=lineBuffer[icol];
	  init=true;
	}
	else{
	  if(minValue>lineBuffer[icol])
	    minValue=lineBuffer[icol];
	  if(maxValue<lineBuffer[icol])
	    maxValue=lineBuffer[icol];
	}
      }
    }
  }
  if(!init)
    throw(static_cast<string>("Warning: not initialized"));
}

void ImgReaderGdal::getMinMax(double& minValue, double& maxValue, int band, bool exhaustiveSearch) const
{
  GDALRasterBand  *poBand;
  int             bGotMin, bGotMax;
  double          adfMinMax[2];
        
  poBand = m_gds->GetRasterBand(band+1);
  adfMinMax[0] = poBand->GetMinimum( &bGotMin );
  adfMinMax[1] = poBand->GetMaximum( &bGotMax );
  if( ! (bGotMin && bGotMax) )
    GDALComputeRasterMinMax((GDALRasterBandH)poBand, FALSE, adfMinMax);
    // GDALComputeRasterMinMax((GDALRasterBandH)poBand, TRUE, adfMinMax);
  minValue=adfMinMax[0];
  maxValue=adfMinMax[1];
  if(exhaustiveSearch){//force exhaustive search
    vector<double> lineBuffer(nrOfCol());
    bool init=false;
    for(int irow=0;irow<nrOfRow();++irow){
      readData(lineBuffer,GDT_Float64,irow,band);
      for(int icol=0;icol<nrOfCol();++icol){
        bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
        if(valid){
          if(!init){
            minValue=lineBuffer[icol];
            maxValue=lineBuffer[icol];
            init=true;
          }
          else{
            if(minValue>lineBuffer[icol])
              minValue=lineBuffer[icol];
            if(maxValue<lineBuffer[icol])
              maxValue=lineBuffer[icol];
          }
        }
      }
    }
    if(!init)
      throw(static_cast<string>("Warning: not initialized"));
  }
}

void ImgReaderGdal::getRange(vector<short>& range, int band) const
{
  vector<short> lineBuffer(nrOfCol());
  range.clear();
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,GDT_Int16,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      if(find(range.begin(),range.end(),lineBuffer[icol])==range.end())
        range.push_back(lineBuffer[icol]);
    }
  }
  sort(range.begin(),range.end());
}

int ImgReaderGdal::getNoDataValues(vector<double>& noDataValues) const
{
  if(m_noDataValues.size()){
    noDataValues=m_noDataValues;
    return m_noDataValues.size();
  }
  else
    return 0;
}

int ImgReaderGdal::pushNoDataValue(double noDataValue)
{
  if(find(m_noDataValues.begin(),m_noDataValues.end(),noDataValue)==m_noDataValues.end())
    m_noDataValues.push_back(noDataValue);
  return(m_noDataValues.size());
}

// bool ImgReaderGdal::setNoDataValue(double noDataValue,int band)
// {
//   GDALRasterBand  *poBand;
//   poBand = m_gds->GetRasterBand(band+1);
//   if(poBand->SetNoDataValue(noDataValue)!=CE_None)
//     return false;
//   else
//     return true;
// }

void ImgReaderGdal::getRefPix(double& refX, double &refY, int band) const
{
  vector<double> lineBuffer(nrOfCol());
  double validCol=0;
  double validRow=0;
  int nvalidCol=0;
  int nvalidRow=0;
  for(int irow=0;irow<nrOfRow();++irow){
    readData(lineBuffer,GDT_Float64,irow,band);
    for(int icol=0;icol<nrOfCol();++icol){
      bool valid=(find(m_noDataValues.begin(),m_noDataValues.end(),lineBuffer[icol])==m_noDataValues.end());
      if(valid){
        validCol+=icol+1;
        ++nvalidCol;
        validRow+=irow+1;
        ++nvalidRow;
      }
    }
  }
  if(isGeoRef()){
    //reference coordinate is lower left corner of pixel in centre of gravity
    //we need geo coordinates for exactly this location: validCol(Row)/nvalidCol(Row)-0.5
    double cgravi=validCol/nvalidCol-0.5;
    double cgravj=validRow/nvalidRow-0.5;
    double refpixeli=floor(cgravi);
    double refpixelj=ceil(cgravj-1);
    //but image2geo provides location at centre of pixel (shifted half pixel right down)
    image2geo(refpixeli,refpixelj,refX,refY);
    //refX and refY now refer to centre of gravity pixel
    refX-=0.5*getDeltaX();//shift to left corner
    refY-=0.5*getDeltaY();//shift to lower left corner
  }
  else{
    refX=floor(validCol/nvalidCol-0.5);//left corner
    refY=floor(validRow/nvalidRow-0.5);//upper corner
    //shift to lower left corner of pixel
    refY+=1;
      
  }
}

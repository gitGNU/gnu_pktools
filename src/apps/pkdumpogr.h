/**********************************************************************
pkdumpogr.h: dump ogr file to text file or standard output
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
#include <map>
#include "base/Vector2d.h"
#include "imageclasses/ImgReaderOgr.h"

#ifndef _PKDUMPOGR_H_
#define _PKDUMPOGR_H_

using namespace std;

template<typename T> unsigned int readDataImageShape(const string &filename,
                                                     map<int,Vector2d<T> > &mapPixels, //[classNr][pixelNr][bandNr],
                                                     vector<string>& fields,
                                                     double start,
                                                     double end,
                                                     const string& label,
                                                     const string& query="",
                                                     int verbose=false);


template<typename T> unsigned int readDataImageShape(const string &filename,
                                                     map<int,Vector2d<T> > &mapPixels, //[classNr][pixelNr][bandNr],
                                                     vector<string>& fields,
                                                     double start,
                                                     double end,
                                                     const string& label,
                                                     const string& query,
                                                     int verbose)
{
  mapPixels.clear();
  int nsample=0;
  int totalSamples=0;  
  int nband=0;
  if(verbose)
    cout << "reading shape file " << filename  << endl;
  ImgReaderOgr imgReaderShape;
  try{
    imgReaderShape.open(filename);
    bool queryFound=false;
    //only retain bands in fields
    imgReaderShape.getFields(fields);
    vector<string>::iterator fit=fields.begin();
    if(verbose)
      cout << "reading fields: ";
    while(fit!=fields.end()){
      if(verbose)
        cout << *fit << " ";
      size_t pos=(*fit).find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ");
      if(pos==string::npos){
        if(query!=""){
          if((*fit).find(query)!=string::npos)
            queryFound=true;
        }
        fields.erase(fit);
      }
      else{
        string fieldname=(*fit).substr(pos);
          int iband=atoi(fieldname.c_str());
          if((start||end)&&(iband<start||iband>end))
            fields.erase(fit);
          else
            ++fit;
      }
    }
    if(verbose)
      cout << endl;
    if(verbose){
      cout << "fields:";
      for(vector<string>::iterator fit=fields.begin();fit!=fields.end();++fit)
        cout << " " << *fit;
      cout << endl;
    }
    if(!nband){
      if(queryFound){
        ostringstream qs;
        qs << "select * from " << imgReaderShape.getLayerName() << " where " << query << "=1";
        if(verbose)
          cout << "reading with sql: " << qs.str() << endl;
        nband=imgReaderShape.readSql(mapPixels,OFTReal,fields,label,qs.str(),NULL,0,true,false);
      }
      else{
        if(verbose)
          cout << "reading data" << endl;
        nband=imgReaderShape.readData(mapPixels,OFTReal,fields,label,0,true,verbose==2);
      }
    }
    else{
      if(queryFound){
        ostringstream qs;
        qs << "select * from " << imgReaderShape.getLayerName() << " where " << query << "=1";
        if(verbose)
          cout << "reading with sql: " << qs.str() << endl;
        assert(nband==imgReaderShape.readSql(mapPixels,OFTReal,fields,label,qs.str(),NULL,0,true,false));
      }
      else        
        assert(nband==imgReaderShape.readData(mapPixels,OFTReal,fields,label,0,true,false));
    }
  }
  catch(string e){
    ostringstream estr;
    estr << e << " " << filename;
    throw(estr.str());
  }
  nsample=imgReaderShape.getFeatureCount();
  totalSamples+=nsample;
  if(verbose)
    cout << ": " << nsample << " samples read with " << nband << " bands" << endl;
  imgReaderShape.close();
  if(verbose)
    cout << "total number of samples read " << totalSamples << endl;
  return totalSamples;
}
#endif //_PKDUMPOGR_H_

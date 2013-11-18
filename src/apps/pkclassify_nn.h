/**********************************************************************
pkclassify_nn.h: classify raster image using Artificial Neural Network
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
#include "imageclasses/ImgReaderOgr.h"
#include "base/Vector2d.h"

#ifndef _PKCLASSIFY_NN_H_
#define _PKCLASSIFY_NN_H_

using namespace std;

template<typename T> unsigned int readDataImageShape(const string &filename,
                                                     map<string,Vector2d<T> > &mapPixels, //[classNr][pixelNr][bandNr],
                                                     vector<string>& fields,
                                                     const vector<short>& bands,
                                                     const string& label,
                                                     int verbose=false);

template<typename T> unsigned int readDataImageShape(const string &filename,
                                                     map<string,Vector2d<T> > &mapPixels, //[classNr][pixelNr][bandNr],
                                                     vector<string>& fields,
                                                     double start,
                                                     double end,
                                                     const string& label,
                                                     int verbose=false);

template<typename T> unsigned int readDataImageShape(const string &filename,
                                                     map<string,Vector2d<T> > &mapPixels, //[classNr][pixelNr][bandNr],
                                                     vector<string>& fields,
                                                     const vector<short>& bands,
                                                     const string& label,
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
    //only retain bands in fields
    imgReaderShape.getFields(fields);
    vector<string>::iterator fit=fields.begin();
    if(verbose>1)
      cout << "reading fields: ";
    while(fit!=fields.end()){
      if(verbose)
        cout << *fit << " ";
      // size_t pos=(*fit).find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ");
      if((*fit).substr(0,1)=="B"){
	if((*fit).substr(1).find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ")!=string::npos){
	  int theBand=atoi((*fit).substr(1).c_str());
	  if(bands.size()){
	    bool validBand=false;
	    for(int iband=0;iband<bands.size();++iband){
	      if(theBand==bands[iband])
		validBand=true;
	    }
	    if(validBand)
	      ++fit;
	    else
	      fields.erase(fit);
	  }
	  else
	    ++fit;
	}
	else if((*fit)=="B" || (*fit)=="Band")//B is only band
	  ++fit;
      }
      else
        fields.erase(fit);
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
      if(verbose)
        cout << "reading data" << endl;
      nband=imgReaderShape.readData(mapPixels,OFTReal,fields,label,0,true,verbose==2);

    }
    else
      assert(nband==imgReaderShape.readData(mapPixels,OFTReal,fields,label,0,true,false));
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

template<typename T> unsigned int readDataImageShape(const string &filename,
                                                     map<string,Vector2d<T> > &mapPixels, //[classNr][pixelNr][bandNr],
                                                     vector<string>& fields,
                                                     double start,
                                                     double end,
                                                     const string& label,
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
    //only retain bands in fields
    imgReaderShape.getFields(fields);
    vector<string>::iterator fit=fields.begin();
    if(verbose)
      cout << "reading fields: ";
    while(fit!=fields.end()){
      if(verbose)
        cout << *fit << " ";
      // size_t pos=(*fit).find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ");
      if((*fit).substr(0,1)=="B"){
	if((*fit).substr(1).find_first_not_of("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ_ ")!=string::npos){
	  int iband=atoi((*fit).substr(1).c_str());
	  if((start||end)&&(iband<start||iband>end))
	    fields.erase(fit);
	  else
	    ++fit;
	}
	else if(*fit=="B" || *fit=="Band")
	  ++fit;
      }
      else
        fields.erase(fit);
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
      if(verbose)
        cout << "reading data" << endl;
      nband=imgReaderShape.readData(mapPixels,OFTReal,fields,label,0,true,verbose==2);

    }
    else
      assert(nband==imgReaderShape.readData(mapPixels,OFTReal,fields,label,0,true,false));
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
#endif //_PKCLASSIFY_NN_H_

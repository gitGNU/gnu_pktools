/**********************************************************************
Filter.cc: class for filtering
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
#include "Filter.h"
#include <assert.h>
#include <iostream>

filter::Filter::Filter(void)
{
}


filter::Filter::Filter(const vector<double> &taps)
  : m_taps(taps)
{
  assert(m_taps.size()%2);
}

void filter::Filter::setTaps(const vector<double> &taps)
{
  m_taps=taps;
  assert(m_taps.size()%2);
}

void filter::Filter::morphology(const ImgReaderGdal& input, ImgWriterGdal& output, int method, int dim, short down, int offset)
{
  Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
  Vector2d<double> lineOutput(input.nrOfBand(),input.nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int y=0;y<input.nrOfRow();++y){
    for(int iband=0;iband<input.nrOfBand();++iband)
      input.readData(lineInput[iband],GDT_Float64,y,iband);
    vector<double> pixelInput(input.nrOfBand());
    vector<double> pixelOutput(input.nrOfBand());
    for(int x=0;x<input.nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      morphology(pixelInput,pixelOutput,method,dim,down,offset);
      for(int iband=0;iband<input.nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(int iband=0;iband<input.nrOfBand();++iband){
      try{
        output.writeData(lineOutput[iband],GDT_Float64,y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::doit(const ImgReaderGdal& input, ImgWriterGdal& output, short down, int offset)
{
  Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
  Vector2d<double> lineOutput(input.nrOfBand(),input.nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int y=0;y<input.nrOfRow();++y){
    for(int iband=0;iband<input.nrOfBand();++iband)
      input.readData(lineInput[iband],GDT_Float64,y,iband);
    vector<double> pixelInput(input.nrOfBand());
    vector<double> pixelOutput(input.nrOfBand());
    for(int x=0;x<input.nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      doit(pixelInput,pixelOutput,down,offset);
      for(int iband=0;iband<input.nrOfBand();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(int iband=0;iband<input.nrOfBand();++iband){
      try{
        output.writeData(lineOutput[iband],GDT_Float64,y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}

void filter::Filter::applyFwhm(const ImgReaderGdal& input, ImgWriterGdal& output, const vector<double> &wavelengthIn, const vector<double> &wavelengthOut, const vector<double> &fwhm, bool verbose){
  Vector2d<double> lineInput(input.nrOfBand(),input.nrOfCol());
  Vector2d<double> lineOutput(wavelengthOut.size(),input.nrOfCol());
  const char* pszMessage;
  void* pProgressArg=NULL;
  GDALProgressFunc pfnProgress=GDALTermProgress;
  double progress=0;
  pfnProgress(progress,pszMessage,pProgressArg);
  for(int y=0;y<input.nrOfRow();++y){
    for(int iband=0;iband<input.nrOfBand();++iband)
      input.readData(lineInput[iband],GDT_Float64,y,iband);
    vector<double> pixelInput(input.nrOfBand());
    vector<double> pixelOutput;
    for(int x=0;x<input.nrOfCol();++x){
      pixelInput=lineInput.selectCol(x);
      applyFwhm<double>(pixelInput,wavelengthIn,pixelOutput,wavelengthOut, fwhm, verbose);
      for(int iband=0;iband<pixelOutput.size();++iband)
        lineOutput[iband][x]=pixelOutput[iband];
    }
    for(int iband=0;iband<output.nrOfBand();++iband){
      try{
        output.writeData(lineOutput[iband],GDT_Float64,y,iband);
      }
      catch(string errorstring){
        cerr << errorstring << "in band " << iband << ", line " << y << endl;
      }
    }
    progress=(1.0+y)/output.nrOfRow();
    pfnProgress(progress,pszMessage,pProgressArg);
  }
}



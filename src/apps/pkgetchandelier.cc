/**********************************************************************
pkgetchandelier.cc: program to optimize model parameters for brdf correction
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
#include <assert.h>
#include <string>
#include <iostream>
#include <nlopt.hpp>
#include "base/PointData.h"
#include "algorithms/StatFactory.h"
#include "imageclasses/ImgReaderOgr.h"
#include "Optionpk.h"
#include "pkgetchandelier.h"

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
  assert(grad.empty());
  double rmse=0;
  vector<vector<PointData> > *blockData=reinterpret_cast<vector< vector <PointData> > *> (my_func_data);
  for(unsigned int ipoint=0;ipoint<blockData->size();++ipoint){
    for(unsigned int ifile1=0;ifile1<blockData->at(ipoint).size()-1;++ifile1){
      double reflectance1=0;
      if(blockData->at(ipoint)[ifile1].getImage()<0)
        reflectance1=blockData->at(ipoint)[ifile1].getReflectance();
      else{
        State state1;
        state1.k=x[blockData->at(ipoint)[ifile1].getImage()*4];
        state1.e=x[blockData->at(ipoint)[ifile1].getImage()*4+1];
        state1.a=x[blockData->at(ipoint)[ifile1].getImage()*4+2];
        state1.haze=x[blockData->at(ipoint)[ifile1].getImage()*4+3];
        reflectance1=blockData->at(ipoint)[ifile1].correctReflectance(state1);
      }
      for(unsigned int ifile2=ifile1;ifile2<blockData->at(ipoint).size();++ifile2){
        double reflectance2=0;
        if(blockData->at(ipoint)[ifile2].getImage()<0)
          reflectance2=blockData->at(ipoint)[ifile2].getReflectance();
        else{
          State state2;
          state2.k=x[blockData->at(ipoint)[ifile2].getImage()*4];
          state2.e=x[blockData->at(ipoint)[ifile2].getImage()*4+1];
          state2.a=x[blockData->at(ipoint)[ifile2].getImage()*4+2];
          state2.haze=x[blockData->at(ipoint)[ifile2].getImage()*4+3];
          reflectance2=blockData->at(ipoint)[ifile2].correctReflectance(state2);
        }
        // if(reflectance2>2*reflectance1||reflectance2<reflectance1/2)
        //   continue;
        double diff=(reflectance1-reflectance2)*(reflectance1-reflectance2);
        if(PointData::m_residual>0){
          if(sqrt(diff/reflectance1/reflectance2)<PointData::m_residual)
            rmse+=diff;
        }
        else
          rmse+=diff;
      }
    }
  }
  rmse=sqrt(rmse/blockData->size());
  // std::cout << "difference: " << diff << std::endl;
  return rmse;
}

int main(int argc, char *argv[])
{
  std::string versionString="version ";
  versionString+=VERSION;
  versionString+=", Copyright (C) 2008-2012 Pieter Kempeneers.\n\
   This program comes with ABSOLUTELY NO WARRANTY; for details type use option -h.\n\
   This is free software, and you are welcome to redistribute it\n\
   under certain conditions; use option --license for details.";
  Optionpk<bool> version_opt("\0","version",versionString,false);
  Optionpk<bool> license_opt("lic","license","show license information",false);
  Optionpk<bool> help_opt("h","help","shows this help info",false);
  Optionpk<bool> todo_opt("\0","todo","",false);
  Optionpk<string> input_opt("i", "input", "Reflectance imageinput (ogr vector) files","");
  // Optionpk<string> sample_opt("s", "sample", "sample (ogr vector) file","");
  Optionpk<double> grid_opt("grid", "grid", "grid for tiepoints (in m). Smaller values provide more tiepoints",100);
  Optionpk<string> geom_opt("g", "geom", "Geometry image files","");
  Optionpk<unsigned short> band_opt("b", "band", "field name of Reflectance band in reflectance file",0);
  Optionpk<unsigned short> sza_opt("sza", "sza", "band number (starting from 0) for Sun Zenith Angle in geometry image file",3);
  // Optionpk<unsigned short> vza_opt("vza", "vza", "band number (starting from 0) for Sun Zenith Angle in geometry image file",5);
  Optionpk<unsigned short> saa_opt("saa", "saa", "band number (starting from 0) for Sun Azimuth Angle in geometry image file",6);
  Optionpk<unsigned short> vaa_opt("vaa", "vaa", "band number (starting from 0) for View Azimuth Angle in geometry image file",7);
  Optionpk<double> tau_opt("tau", "tau", "Optical depth",0.2);
  Optionpk<double> deltaT_opt("deltaT", "deltaT", "Exposure time",1.0);
  Optionpk<double> residual_opt("res", "residual", "residual error for filtering outlier tiepoints: maximum relative error to take tie point into account (use 0 if not filtering is required)",0.0);
  Optionpk<unsigned int> maxit_opt("maxit","maxit","maximum number of iterations",500);
  Optionpk<double> tolerance_opt("tol","tolerance","relative tolerance for stopping criterion",0.0001);
  Optionpk<double> init_opt("is","init","initial state vector: -is k -is e -is a -s haze",0);
  Optionpk<double> lb_opt("lb","lb","lower bounds for k,e and a: use -lb k -lb e -lb a -lb haze",0);
  Optionpk<double> ub_opt("ub","ub","upper bounds for k,e and a: use -ub k -ub e -ub a -ub haze",0);
  Optionpk<short> dim_opt("dim","dim","window size to calculate mean reflectance (use odd number, e.g., 3, 5, 7)",3);
  Optionpk<double> var_opt("var","var","maximum variance in window to take tiepoint into account (use 0 if don't care)",0);
  Optionpk<string> mask_opt("m", "mask", "Mask image(s) (single mask for all input images or one mask for each input image", "");
  Optionpk<double> invalid_opt("t", "invalid", "Mask value where image is invalid.", 0);
  // Optionpk<bool> single_opt("\0","single","create virtual tiepoints for single images with no overlap",false);
  Optionpk<unsigned long int> mintp_opt("min", "min", "Minimum number of tiepoints", 1000);
  Optionpk<short> verbose_opt("v", "verbose", "verbose", 0);

  version_opt.retrieveOption(argc,argv);
  license_opt.retrieveOption(argc,argv);
  help_opt.retrieveOption(argc,argv);
  todo_opt.retrieveOption(argc,argv);

  input_opt.retrieveOption(argc,argv);
  // sample_opt.retrieveOption(argc,argv);
  grid_opt.retrieveOption(argc,argv);
  geom_opt.retrieveOption(argc,argv);
  band_opt.retrieveOption(argc,argv);
  sza_opt.retrieveOption(argc,argv);
  saa_opt.retrieveOption(argc,argv);
  vaa_opt.retrieveOption(argc,argv);
  tau_opt.retrieveOption(argc,argv);
  deltaT_opt.retrieveOption(argc,argv);
  residual_opt.retrieveOption(argc,argv);
  maxit_opt.retrieveOption(argc,argv);
  tolerance_opt.retrieveOption(argc,argv);
  init_opt.retrieveOption(argc,argv);
  lb_opt.retrieveOption(argc,argv);
  ub_opt.retrieveOption(argc,argv);
  dim_opt.retrieveOption(argc,argv);
  var_opt.retrieveOption(argc,argv);
  mask_opt.retrieveOption(argc,argv);
  invalid_opt.retrieveOption(argc,argv);
  // single_opt.retrieveOption(argc,argv);
  verbose_opt.retrieveOption(argc,argv);
  mintp_opt.retrieveOption(argc,argv);

  if(version_opt[0]||todo_opt[0]){
    cout << version_opt.getHelp() << endl;
    cout << "todo: " << todo_opt.getHelp() << endl;
    exit(0);
  }
  if(license_opt[0]){
    cout << Optionpk<bool>::getGPLv3License() << endl;
    exit(0);
  }
  if(help_opt[0]){
    cout << "usage: pktrirad -s tiepoints -i input [-i input] -g geom [-g geom] [OPTIONS]" << endl;
    exit(0);
  }

  assert(geom_opt.size()==input_opt.size());

  assert(init_opt.size()==4*input_opt.size());//k,e,a,haze
  assert(lb_opt.size()==4*input_opt.size());
  assert(ub_opt.size()==4*input_opt.size());

  PointData::m_tau=tau_opt[0];
  PointData::m_deltaT=deltaT_opt[0];
  PointData::m_residual=residual_opt[0];
  //construct list of points defined in sample file
  vector<vector< PointData> > blockData;//[point][image]
  vector<PointData> pdVector;//list of points (contains all image data for this tiepoint)
  //get union bounding box
  double maxLRX=0;
  double maxULY=0;
  double minULX=0;
  double minLRY=0;
  for(int ifile=0;ifile<input_opt.size();++ifile){
    double theULX, theULY, theLRX, theLRY;
    if(verbose_opt[0]>1)
      std::cout << "opening input file " << input_opt[ifile] << std::endl;
    ImgReaderGdal inputReader(input_opt[ifile]);
    inputReader.getBoundingBox(theULX,theULY,theLRX,theLRY);
    if(verbose_opt[0])
      std::cout << setprecision(12) << "--ulx=" << theULX << " --uly=" << theULY << " --lrx=" << theLRX << " --lry=" << theLRY << " ";
    if(!ifile){
      maxLRX=theLRX;
      maxULY=theULY;
      minULX=theULX;
      minLRY=theLRY;
    }
    else{
      maxLRX=(theLRX>maxLRX)?theLRX:maxLRX;
      maxULY=(theULY>maxULY)?theULY:maxULY;
      minULX=(theULX<minULX)?theULX:minULX;
      minLRY=(theLRY<minLRY)?theLRY:minLRY;
    }
    inputReader.close();
  }
  if(verbose_opt[0])
    std::cout << "union bounding box: " << setprecision(12) << "--ulx=" << minULX << " --uly=" << maxULY << " --lrx=" << maxLRX << " --lry=" << minLRY << std::endl;

  if(verbose_opt[0])
    cout << "number of mask images: " << mask_opt.size() << endl;
  vector<double> oldRowMask(mask_opt.size());
  int imask=0;
  ImgReaderGdal maskReader;
  if(mask_opt[0]!=""){
    if(mask_opt.size()!=input_opt.size())//single mask for all input images: open now
      maskReader.open(mask_opt[imask]);
  }
  for(double y=maxULY-grid_opt[0]/2.0;y>minLRY;y-=grid_opt[0]){
    for(double x=minULX+grid_opt[0]/2.0;x<maxLRX;x+=grid_opt[0]){
      for(int ifile=0;ifile<input_opt.size();++ifile){
        PointData pd;
        pd.setImage(ifile);
        if(verbose_opt[0]>1)
          std::cout << "opening input file " << input_opt[ifile] << std::endl;
        ImgReaderGdal inputReader(input_opt[ifile]);
        if(verbose_opt[0]>1)
          std::cout << "opening geom input file " << geom_opt[ifile] << std::endl;
        ImgReaderGdal geoReader(geom_opt[ifile]);
        double reflectance;
        double col,row;
        inputReader.geo2image(x,y,col,row);
        if(col<0||row<0||col>inputReader.nrOfCol()||row>inputReader.nrOfRow()){
          inputReader.close();
          geoReader.close();
          continue;
        }
        //begin mask
        if(mask_opt[0]!=""){
          bool masked=false;
          if(mask_opt.size()==input_opt.size())//one mask for each input image: open now
            maskReader.open(mask_opt[ifile]);
          double colMask,rowMask;
          maskReader.geo2image(x,y,colMask,rowMask);
          colMask=static_cast<int>(colMask);
          rowMask=static_cast<int>(rowMask);
          double maskValue=0;
          if(rowMask>=0&&rowMask<maskReader.nrOfRow()&&colMask>=0&&colMask<maskReader.nrOfCol()){
            try{
              maskReader.readData(maskValue,GDT_Float64,static_cast<int>(colMask),static_cast<int>(rowMask));
            }
            catch(string errorstring){
              cerr << errorstring << endl;
              exit(1);
            }
            if(maskValue==invalid_opt[0])
              masked=true;
          }
          if(mask_opt.size()==input_opt.size())//one mask for each input image: close now
            maskReader.close();
          if(masked){
            inputReader.close();
            geoReader.close();
            continue;
          }
        }
        //end mask
        statfactory::StatFactory stat;
        vector<double> windowBuffer;
        for(int windowJ=-dim_opt[0]/2;windowJ<(dim_opt[0]+1)/2;++windowJ){
          double j=row+windowJ;
          if(j<0||j>=inputReader.nrOfRow())
            continue;
          for(int windowI=-dim_opt[0]/2;windowI<(dim_opt[0]+1)/2;++windowI){
            double i=col+windowI;
            if(i<0||i>=inputReader.nrOfCol())
              continue;
            inputReader.readData(reflectance,GDT_Float64,i,j,band_opt[0]);
            windowBuffer.push_back(reflectance);
          }
        }
        double variance=0;
        if(windowBuffer.size()>1)
          stat.meanVar(windowBuffer,reflectance,variance);
        if(var_opt[0]>0){
          if(verbose_opt[0]>1)
            std::cout << "reflectance (" << windowBuffer.size() << " at " << col << "," << row << ") mean, var: " << reflectance << ", " << variance << std::endl;
          if(variance>var_opt[0]){
            inputReader.close();
            geoReader.close();
            continue;
          }
        }
        pd.setReflectance(reflectance);
        double sza;
        double saa;
        double vaa;
        geoReader.geo2image(x,y,col,row);
        geoReader.readData(sza,GDT_Float64,col,row,sza_opt[0]);
        geoReader.readData(saa,GDT_Float64,col,row,saa_opt[0]);
        geoReader.readData(vaa,GDT_Float64,col,row,vaa_opt[0]);
        double theRelativeAzimuth=saa-vaa;
        pd.setCosFi(theRelativeAzimuth);
        pd.setCosSZA(sza);
        double centreX,centreY;
        inputReader.getCentrePos(centreX,centreY);
        double ulx=inputReader.getUlx();
        double uly=inputReader.getUly();
        double r0=((ulx-centreX)*(ulx-centreX)+(uly-centreY)*(uly-centreY));
        double r=((x-centreX)*(x-centreX)+(y-centreY)*(y-centreY));
        pd.setR(r/r0);
        pdVector.push_back(pd);
        inputReader.close();
        geoReader.close();
      }
      if(pdVector.size()>1)//tiepoints must cover at least two images
        blockData.push_back(pdVector);
      // else if(single_opt[0]){
      //   PointData pd=pdVector.back();
      //   pd.setImage(-1);
      //   pdVector.push_back(pd);
      //   blockData.push_back(pdVector);
      // }
      pdVector.clear();
      if(mask_opt.size()!=input_opt.size())//single mask for all input images: close now
        maskReader.close();
    }
  }

  if(verbose_opt[0])
    std::cout << "Number of tiepoints: " << blockData.size() << std::endl;
  if(blockData.size()<mintp_opt[0]){
    std::cerr << "Error: Number of tiepoints=" << blockData.size() << " is smaller than " << mintp_opt[0] << std::endl;
    exit;
  }

  //todo: make nlopt::LN_COBYLA as a command line option
  //nlopt::opt opt(nlopt::LN_COBYLA,4*input_opt.size());//k,e,a,haze
  nlopt::opt optimizer(nlopt::LN_SBPLX,4*input_opt.size());//k,e,a,haze

  optimizer.set_min_objective(objFunction, &blockData);

  if(verbose_opt[0])
    std::cout << "set lower and upper bounds" << std::endl;
  optimizer.set_lower_bounds(lb_opt);
  optimizer.set_upper_bounds(ub_opt);
  
  if(verbose_opt[0])
    std::cout << "set stopping criteria" << std::endl;
  //set stopping criteria
  if(maxit_opt[0])
    optimizer.set_maxeval(maxit_opt[0]);
  else
    optimizer.set_xtol_rel(tolerance_opt[0]);
  double minf=0;
  std::vector<double> x=init_opt;

  if(verbose_opt[0])
    std::cout << "optimizing with " << optimizer.get_algorithm_name() << "..." << std::endl;
  optimizer.optimize(x, minf);
  for(int index=0;index<x.size();++index){
    std::cout << " -s " << x[index];
    if(index%4==3)
      std::cout << std::endl;
  }
  if(verbose_opt[0])
    std::cout << "minf: " << minf << std::endl;
}


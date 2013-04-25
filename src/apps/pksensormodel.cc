/**********************************************************************
pksensormodel.cc: program to calculate geometric position based on row (sensor), col (sensor), roll, pitch, yaw and lens coordinates
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
#include <sstream>
#include <fstream>
#include <vector>
#include <math.h>
#include <nlopt.hpp>
#include "base/Optionpk.h"
#include "algorithms/OptFactory.h"
#include "algorithms/StatFactory.h"
#include "pksensormodel.h"

double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data){
  assert(grad.empty());
  double error=0;
  DataModel *dm=reinterpret_cast<DataModel*> (my_func_data);
  arma::vec bc_att(3);
  bc_att(0)=x[0];
  bc_att(1)=x[1];
  bc_att(2)=x[2];
  dm->setBoresightAtt(bc_att);
  for(unsigned int ipoint=0;ipoint<dm->getSize();++ipoint){
    double e=dm->getDistGeo(ipoint);
    error+=e;
  }
  error/=dm->getSize();
  return error;
}

int main(int argc, char *argv[])
{
  Optionpk<string> input_opt("i","input","name of the input text file");
  Optionpk<string> datum_opt("datum","datum","GPS datum of the input points","WGS84");
  Optionpk<int> s_srs_opt("s_srs","s_srs","source EPSG (integer) code of GCP input coordinates",4326);
  Optionpk<int> t_srs_opt("t_srs","t_srs","target EPSG (integer) code to project output coordinates",4326);
  Optionpk<double> focal_opt("f","focal","focal length of the camera (in m)",0.041517);
  Optionpk<double> dx_opt("dx","dx","pixel size in CCD (in micro m)",20);
  Optionpk<double> dy_opt("dy","dy","pixel size in CCD (in micro m)",20);
  Optionpk<double> x_opt("x","x","gcp x coordinate)");
  Optionpk<double> y_opt("y","y","gcp y coordinate)");
  Optionpk<double> z_opt("z","z","gcp z coordinate in m)");
  Optionpk<double> errorZ_opt("ez","ez","offset (error) in z coordinate of the GCP");
  Optionpk<int> col_opt("c","col","column of the pixel on CCD");
  Optionpk<int> row_opt("r","row","row of the pixel on CCD");
  Optionpk<double> roll_opt("roll","roll","platform attitude roll");
  Optionpk<double> pitch_opt("pitch","pitch","platform attitude pitch");
  Optionpk<double> yaw_opt("yaw","yaw","platform attitude yaw");
  Optionpk<double> xl_opt("xl","xl","platform x position in epsg:4326");
  Optionpk<double> yl_opt("yl","yl","platform y position in epsg:4326");
  Optionpk<double> zl_opt("zl","zl","platform z position (in m)");
  // Optionpk<double> bcpos_opt("bcp","bcp","bore sight offset position",0);
  Optionpk<double> bcatt_opt("bca","bca","bore sight attitude calibration offset for roll, pitch, yaw of the platform",0);
  Optionpk<double> nx_opt("nx","nx","number of columns in CCD",1441);
  Optionpk<double> ny_opt("ny","ny","number of rows in CCD",1);
  Optionpk<double> fov_opt("fov","fov","field of view (in degrees)",39.74584);
  Optionpk<double> ppx_opt("ppx","ppx","scan line principal point in x",711);
  Optionpk<double> ppy_opt("ppy","ppy","scan line principal point in y",0);
  Optionpk<char> fs_opt("fs","fs","field separator.",' ');
  Optionpk<string> output_opt("o", "output", "Output ascii file (empty: use stdout");
  Optionpk<short> sensor_opt("s", "sensor", "Sensor type (0: whiskbroom, 1: pushbroom, 2: frame, 3: apex (predefined settings))",0);
  Optionpk<double> polynome_opt("pol","polynome","coefficients for polynome to correct accross track");
  Optionpk<bool> mean_opt("m","mean","calculate mean error",false);
  Optionpk<bool> optimize_opt("opt","opt","optimize bore sight angles using GCP points in input file",false);
  Optionpk<double> lb_opt("lb","lb","lower bounds for offset bore sight angles: use -lb offset_roll -lb offset_pitch -lb offset_yaw",0);
  Optionpk<double> ub_opt("ub","ub","upper bounds for offset bore sight angles: use -ub offset_roll -ub offset_pitch -ub offset_yaw",0);
  Optionpk<unsigned int> maxit_opt("maxit","maxit","maximum number of iterations",500);
  Optionpk<double> tolerance_opt("tol","tolerance","relative tolerance for stopping criterion",0.0001);
  Optionpk<double> init_opt("is","init","initial state vector for bore sight angles: -is init_roll -is init_pitch -is init_yaw",0);
  Optionpk<double> threshold_opt("t","t","threshold for GCP error (in m)",0);
  Optionpk<bool> gcprad_opt("gcprad","gcprad","gcp coordinates",false);
  Optionpk<bool> pplrad_opt("prad","prad","platform pos coordinates",false);
  Optionpk<bool> aplrad_opt("arad","arad","platform attitude angles",false);
  Optionpk<bool> bcrad_opt("brad","brad","boresight attitude angles",false);
  Optionpk<bool> getzenith_opt("gz","getzenith","get zenith angle from platform",false);
  Optionpk<string> algorithm_opt("a", "algorithm", "optimization algorithm (see http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms)","LN_COBYLA"); 
  Optionpk<short> verbose_opt("v", "verbose", "verbose mode when > 0", 0);

  bool doProcess;//stop process when program was invoked with help option (-h --help)
  try{
    doProcess=input_opt.retrieveOption(argc,argv);
    datum_opt.retrieveOption(argc,argv);
    s_srs_opt.retrieveOption(argc,argv);
    t_srs_opt.retrieveOption(argc,argv);
    focal_opt.retrieveOption(argc,argv);
    dx_opt.retrieveOption(argc,argv);
    dy_opt.retrieveOption(argc,argv);
    nx_opt.retrieveOption(argc,argv);
    ny_opt.retrieveOption(argc,argv);
    x_opt.retrieveOption(argc,argv);
    y_opt.retrieveOption(argc,argv);
    z_opt.retrieveOption(argc,argv);
    errorZ_opt.retrieveOption(argc,argv);
    xl_opt.retrieveOption(argc,argv);
    yl_opt.retrieveOption(argc,argv);
    zl_opt.retrieveOption(argc,argv);
    roll_opt.retrieveOption(argc,argv);
    pitch_opt.retrieveOption(argc,argv);
    yaw_opt.retrieveOption(argc,argv);
    col_opt.retrieveOption(argc,argv);
    row_opt.retrieveOption(argc,argv);
    // bcpos_opt.retrieveOption(argc,argv);
    bcatt_opt.retrieveOption(argc,argv);
    fov_opt.retrieveOption(argc,argv);
    ppx_opt.retrieveOption(argc,argv);
    ppy_opt.retrieveOption(argc,argv);
    fs_opt.retrieveOption(argc,argv);
    output_opt.retrieveOption(argc,argv);
    sensor_opt.retrieveOption(argc,argv);
    polynome_opt.retrieveOption(argc,argv);
    mean_opt.retrieveOption(argc,argv);
    optimize_opt.retrieveOption(argc,argv);
    lb_opt.retrieveOption(argc,argv);
    ub_opt.retrieveOption(argc,argv);
    maxit_opt.retrieveOption(argc,argv);
    tolerance_opt.retrieveOption(argc,argv);
    init_opt.retrieveOption(argc,argv);
    threshold_opt.retrieveOption(argc,argv);
    gcprad_opt.retrieveOption(argc,argv);
    pplrad_opt.retrieveOption(argc,argv);
    aplrad_opt.retrieveOption(argc,argv);
    bcrad_opt.retrieveOption(argc,argv);
    getzenith_opt.retrieveOption(argc,argv);
    algorithm_opt.retrieveOption(argc,argv);
    verbose_opt.retrieveOption(argc,argv);
  }
  catch(string predefinedString){
    std::cout << predefinedString << std::endl;
    exit(0);
  }
  if(!doProcess){
    std::cout << "short option -h shows basic options only, use long option --help to show all options" << std::endl;
    exit(0);//help was invoked, stop processing
  }

  //ID, row of GCP, column of GCP, x_GCP, y_GCP, Z_GCP, x_platform, y_platform, z_platform, platform roll, platform pitch, platform yaw of point.
  vector<int> vID;
  SensorModel::SensorModel theModel;
  
  switch(sensor_opt[0]){
  case(0):
    theModel.setModel(SensorModel::WHISKBROOM);
    break;
  case(1):
    theModel.setModel(SensorModel::PUSHBROOM);
    break;
  case(2):
    theModel.setModel(SensorModel::FRAME);
    break;
  case(3)://APEX
    theModel.setModel(SensorModel::WHISKBROOM);
    fov_opt[0]=27.984271234;
    nx_opt[0]=1000;
    ny_opt[0]=1;
    dx_opt[0]=40;
    dy_opt[0]=0;
    focal_opt[0]=0.08103;
    ppx_opt[0]=498;
    ppy_opt[0]=0;
    break;
  default:
    std::cerr << "Error: sensor type " << sensor_opt[0] << " not supported" << std::endl;
    exit(0);
    break;
  }
  theModel.setFOV(fov_opt[0]);
  theModel.setNcol(nx_opt[0]);
  theModel.setNrow(ny_opt[0]);
  theModel.setDx(dx_opt[0]);
  theModel.setDy(dy_opt[0]);
  theModel.setF(focal_opt[0]);
  theModel.setPPx(ppx_opt[0]);
  theModel.setPPy(ppy_opt[0]);
  theModel.setPolynome(polynome_opt);
  theModel.setDatum(datum_opt[0]);
  // arma::vec bc_pos(3);
  arma::vec bc_att(3);
  // while(bcpos_opt.size()<3)
  //   bcpos_opt.push_back(bcpos_opt[0]);
  while(bcatt_opt.size()<3)
    bcatt_opt.push_back(bcatt_opt[0]);
  if(bcrad_opt[0]){
    // bc_pos[0]=theModel.rad2deg(bcpos_opt[0]);
    // bc_pos[1]=theModel.rad2deg(bcpos_opt[1]);
    bc_att(0)=theModel.rad2deg(bcatt_opt[0]);
    bc_att(1)=theModel.rad2deg(bcatt_opt[1]);
    bc_att(2)=theModel.rad2deg(bcatt_opt[2]);
  }
  else{
    // bc_pos(0)=bcpos_opt(0);
    // bc_pos(1)=bcpos_opt(1);
    bc_att(0)=bcatt_opt[0];
    bc_att(1)=bcatt_opt[1];
    bc_att(2)=bcatt_opt[2];
  }
  // bc_pos[2]=bcpos_opt[2];
  // theModel.setBoresightPos(bc_pos);
  theModel.setBoresightAtt(bc_att);
  if(verbose_opt[0]>1){
    // std::cout << "bore sight position offset: " << theModel.getBoresightPos() << std::endl;
    std::cout << "bore sight attitude offset: " << theModel.getBoresightAtt() << std::endl;
  }
  DataModel theDataModel(theModel);
  theDataModel.setThreshold(threshold_opt[0]);

  ifstream dataFile;

  int nrow=0;
  int ivalue=0;
  double dvalue=0;

  OGRSpatialReference oSourceSRS, oTargetSRS;
  OGRCoordinateTransformation *poCT;
  //input for sensor model should be universal lat lon (epsg:4326)
  oSourceSRS.importFromEPSG(s_srs_opt[0]);
  oTargetSRS.importFromEPSG(4326);
  poCT = OGRCreateCoordinateTransformation( &oSourceSRS,
                                            &oTargetSRS );

  //frame ID, row, column, GCP X, GCP Y, GCP Z, platform X, platformY, platform Z, platform roll, platform pitch, platform yaw.
  for(int ifile=0;ifile<input_opt.size();++ifile){
    if(verbose_opt[0]>1)
      std::cout << "opening file " << input_opt[ifile] << std::endl;
    dataFile.open(input_opt[ifile].c_str());
    assert(dataFile);
    if(fs_opt[0]>' '&&fs_opt[0]<='~'){//field separator is a regular character (minimum ASCII code is space, maximum ASCII code is tilde)
      string csvRecord;
      while(getline(dataFile,csvRecord)){//read a line
        istringstream csvstream(csvRecord);
        string item;
        int icol=0;//first column is id
        arma::vec thePosGCP(3);
        arma::vec thePosPlatform(3);
        arma::vec theAttPlatform(3);
        while(getline(csvstream,item,fs_opt[0])){//read a column
          if(verbose_opt[0]>1)
            std::cout << item << " ";
          switch(icol){
          case(0)://id
            ivalue=atoi(item.c_str());
            vID.push_back(ivalue);
            break;
          case(1)://row
            ivalue=atoi(item.c_str());
            theDataModel.pushRow(ivalue);
            break;
          case(2)://col
            ivalue=atoi(item.c_str());
            theDataModel.pushCol(ivalue);
            break;
          case(3)://X
            dvalue=atof(item.c_str());
            if(gcprad_opt[0])
              thePosGCP(0)=theModel.rad2deg(dvalue);
            else
              thePosGCP(0)=dvalue;
            break;
          case(4)://Y
            dvalue=atof(item.c_str());
            if(gcprad_opt[0])
              thePosGCP(1)=theModel.rad2deg(dvalue);
            else
              thePosGCP(1)=dvalue;
            break;
          case(5)://Z
            dvalue=atof(item.c_str());
            thePosGCP(2)=dvalue;
            break;
          case(6)://XL
            dvalue=atof(item.c_str());
            if(pplrad_opt[0])
              thePosPlatform(0)=theModel.rad2deg(dvalue);
            else
              thePosPlatform(0)=dvalue;
            break;
          case(7)://YL
            dvalue=atof(item.c_str());
            if(pplrad_opt[0])
              thePosPlatform(1)=theModel.rad2deg(dvalue);
            else
              thePosPlatform(1)=dvalue;
            break;
          case(8)://ZL
            dvalue=atof(item.c_str());
            thePosPlatform(2)=dvalue;
            break;
          case(9)://Roll
            dvalue=atof(item.c_str());
            if(aplrad_opt[0])
              theAttPlatform(0)=theModel.rad2deg(dvalue);
            else
              theAttPlatform(0)=dvalue;
            break;
          case(10)://Pitch
            dvalue=atof(item.c_str());
            if(aplrad_opt[0])
              theAttPlatform(1)=theModel.rad2deg(dvalue);
            else
              theAttPlatform(1)=dvalue;
            break;
          case(11)://Yaw
            dvalue=atof(item.c_str());
            if(aplrad_opt[0])
              theAttPlatform(2)=theModel.rad2deg(dvalue);
            else
              theAttPlatform(2)=dvalue;
            break;
          case(12)://track
            break;
          default:
            std::cerr << "Error: too many columns in input file " << input_opt[ifile] << std::endl;
            break;
          }
          ++icol;
        }
        if(s_srs_opt[0]!=4326){
          if(verbose_opt[0]>1)
            std::cout << "transforming: " << thePosGCP(0) << ", " << thePosGCP(1) << ", " << thePosGCP(2) << std::endl;
          poCT->Transform(1,&(thePosGCP(0)),&(thePosGCP(1)),&(thePosGCP(2)));
          if(verbose_opt[0]>1)
            std::cout << "into: " << thePosGCP(0) << ", " << thePosGCP(1) << ", " << thePosGCP(2) << std::endl;
          // poCT->Transform(1,&(thePosPlatform[0]),&(thePosPlatform[1]),&(thePosPlatform[2]));
        }
        theDataModel.pushPosGCP(thePosGCP);
        theDataModel.pushPosPlatform(thePosPlatform);
        theDataModel.pushAttPlatform(theAttPlatform);
        if(verbose_opt[0]>1)
          std::cout << endl;
        ++nrow;
      }
    }//comma separated value (or ASCII character other than space or TAB)
    else{//space or tab delimited fields
      string spaceRecord;
      while(!getline(dataFile, spaceRecord).eof()){
        if(verbose_opt[0]>1)
          std::cout << spaceRecord << std::endl;
        istringstream lineStream(spaceRecord);
        string item;
        int icol=0;
        arma::vec thePosGCP(3);
        arma::vec thePosPlatform(3);
        arma::vec theAttPlatform(3);
        while(lineStream >> item){
          if(verbose_opt[0]>1)
            std::cout << item << " ";
          istringstream itemStream(item);
          switch(icol){
          case(0)://id
            itemStream >> ivalue;
            vID.push_back(ivalue);
            break;
          case(1)://row
            itemStream >> ivalue;
            theDataModel.pushRow(ivalue);
            break;
          case(2)://col
            itemStream >> ivalue;
            theDataModel.pushCol(ivalue);
            break;
          case(3)://X
            itemStream >> dvalue;
            if(gcprad_opt[0])
              thePosGCP(0)=theModel.rad2deg(dvalue);
            else
              thePosGCP(0)=dvalue;
            break;
          case(4)://Y
            itemStream >> dvalue;
            if(gcprad_opt[0])
              thePosGCP(1)=theModel.rad2deg(dvalue);
            else
              thePosGCP(1)=dvalue;
            break;
          case(5)://Z
            itemStream >> dvalue;
            thePosGCP(2)=dvalue;
            break;
          case(6)://XL
            itemStream >> dvalue;
            if(pplrad_opt[0])
              thePosPlatform(0)=theModel.rad2deg(dvalue);
            else
              thePosPlatform(0)=dvalue;
            break;
          case(7)://YL
            itemStream >> dvalue;
            if(pplrad_opt[0])
              thePosPlatform(1)=theModel.rad2deg(dvalue);
            else
              thePosPlatform(1)=dvalue;
            break;
          case(8)://ZL
            itemStream >> dvalue;
            thePosPlatform(2)=dvalue;
            break;
          case(9)://Roll
            itemStream >> dvalue;
            if(aplrad_opt[0])
              theAttPlatform(0)=theModel.rad2deg(dvalue);
            else
              theAttPlatform(0)=dvalue;
            break;
          case(10)://Pitch
            itemStream >> dvalue;
            if(aplrad_opt[0])
              theAttPlatform(1)=theModel.rad2deg(dvalue);
            else
              theAttPlatform(1)=dvalue;
            break;
          case(11)://Yaw
            itemStream >> dvalue;
            if(aplrad_opt[0])
              theAttPlatform(2)=theModel.rad2deg(dvalue);
            else
              theAttPlatform(2)=dvalue;
            break;
          case(12)://track
            break;
          default:
            std::cerr << "Error: too many columns in input file " << input_opt[ifile] << std::endl;
            break;
          }
          ++icol;
        }
        if(s_srs_opt[0]!=4326){
          if(verbose_opt[0]>1)
            std::cout << "transforming: " << thePosGCP(0) << ", " << thePosGCP(1) << ", " << thePosGCP(2) << std::endl;
          poCT->Transform(1,&(thePosGCP(0)),&(thePosGCP(1)),&(thePosGCP(2)));
          if(verbose_opt[0]>1)
            std::cout << "into: " << thePosGCP(0) << ", " << thePosGCP(1) << ", " << thePosGCP(2) << std::endl;
          // poCT->Transform(1,&(thePosPlatform[0]),&(thePosPlatform[1]),&(thePosPlatform[2]));
        }
        theDataModel.pushPosGCP(thePosGCP);
        theDataModel.pushPosPlatform(thePosPlatform);
        theDataModel.pushAttPlatform(theAttPlatform);
        if(verbose_opt[0]>1)
          std::cout << std::endl;
        if(verbose_opt[0]>1)
          std::cout << "number of columns: " << icol << std::endl;
        ++nrow;
      }
    }
    //todo: assert sizes are all equal...
    dataFile.close();
  }//for ifile

  for(int ipoint=0;ipoint<xl_opt.size();++ipoint){
    if(verbose_opt[0]>1)
      std::cout << "ipoint: " << ipoint << std::endl;
    vID.push_back(ipoint);
    arma::vec thePosPlatform(3);
    arma::vec theAttPlatform(3);
    assert(row_opt.size()>ipoint);
    assert(col_opt.size()>ipoint);
    theDataModel.pushRow(row_opt[ipoint]);
    theDataModel.pushCol(col_opt[ipoint]);
    if(z_opt.size()>ipoint){
      // assert(x_opt.size()>ipoint);
      // assert(y_opt.size()>ipoint);
      arma::vec thePosGCP(3);
      if(x_opt.size()>ipoint&&y_opt.size()>ipoint){
        if(gcprad_opt[0]){
          thePosGCP(0)=theModel.rad2deg(x_opt[ipoint]);
          thePosGCP(1)=theModel.rad2deg(y_opt[ipoint]);
        }
        else{
          thePosGCP(0)=x_opt[ipoint];
          thePosGCP(1)=y_opt[ipoint];
        }
      }
      thePosGCP(2)=z_opt[ipoint];
      if(s_srs_opt[0]!=4326)
        poCT->Transform(1,&(thePosGCP(0)),&(thePosGCP(1)),&(thePosGCP(2)));
      theDataModel.pushPosGCP(thePosGCP);
    }
    assert(yl_opt.size()>ipoint);
    assert(zl_opt.size()>ipoint);
    if(pplrad_opt[0]){
      thePosPlatform(0)=theModel.rad2deg(xl_opt[ipoint]);
      thePosPlatform(1)=theModel.rad2deg(yl_opt[ipoint]);
    }
    else{
      thePosPlatform(0)=xl_opt[ipoint];
      thePosPlatform(1)=yl_opt[ipoint];
    }
    thePosPlatform(2)=zl_opt[ipoint];
    // if(s_srs_opt[0]!=4326)
    //   poCT->Transform(1,&(thePosPlatform[0]),&(thePosPlatform[1]),&(thePosPlatform[2]));
    theDataModel.pushPosPlatform(thePosPlatform);

    assert(roll_opt.size()>ipoint);
    assert(pitch_opt.size()>ipoint);
    assert(yaw_opt.size()>ipoint);
    if(aplrad_opt[0]){
      theAttPlatform(0)=theModel.rad2deg(roll_opt[ipoint]);
      theAttPlatform(1)=theModel.rad2deg(pitch_opt[ipoint]);
      theAttPlatform(2)=theModel.rad2deg(yaw_opt[ipoint]);
    }
    else{
      theAttPlatform(0)=roll_opt[ipoint];
      theAttPlatform(1)=pitch_opt[ipoint];
      theAttPlatform(2)=yaw_opt[ipoint];
    }
    theDataModel.pushAttPlatform(theAttPlatform);
  }

  //todo: remove GCP points with error above threshold
  unsigned int nremoved=0;
  if(threshold_opt[0]>0){
    for(int ipoint=0;ipoint<vID.size();++ipoint){
      if(verbose_opt[0]>1)
        std::cout << "point: " << ipoint << std::endl;
      arma::vec pos_platform=theDataModel.getPosPlatform(ipoint);
      arma::vec att_platform=theDataModel.getAttPlatform(ipoint);
      arma::vec pos_calc=theDataModel.getPos(ipoint);
      arma::vec pos_gcp=theDataModel.getPosGCP(ipoint);
      ostringstream gcpss;
      double e=theDataModel.getModel().getDistGeo(pos_gcp,pos_calc);
      if(e>=theDataModel.getThreshold()){
        vID.erase(vID.begin()+ipoint);
        theDataModel.erase(ipoint);
        ++nremoved;
      }
    }
    if(!theDataModel.getSize())
      std::cerr << "Error: no data after filtering, check if input is correcly set (rad or degrees) or try to lower threshold" << std::endl;
  }

  if(optimize_opt[0]){
    while(lb_opt.size()<3)
      lb_opt.push_back(lb_opt[0]);
    while(ub_opt.size()<3)
      ub_opt.push_back(ub_opt[0]);
    while(init_opt.size()<3)
      init_opt.push_back(init_opt[0]);
    //todo: make nlopt::LN_COBYLA as a command line option
    //nlopt::opt opt(nlopt::LN_COBYLA,4*input_opt.size());//k,e,a,haze
    nlopt::opt optimizer=OptFactory::getOptimizer(algorithm_opt[0],3);//bore sight angles

    optimizer.set_min_objective(objFunction, &theDataModel);

    if(verbose_opt[0]>1)
      std::cout << "set lower and upper bounds" << std::endl;
    assert(lb_opt.size()==ub_opt.size());
    for(int index=0;index<lb_opt.size();++index){
      if(bcrad_opt[0]){
        lb_opt[index]=theModel.rad2deg(lb_opt[index]);
        ub_opt[index]=theModel.rad2deg(ub_opt[index]);
      }
    }
    optimizer.set_lower_bounds(lb_opt);
    optimizer.set_upper_bounds(ub_opt);
  
    if(verbose_opt[0]>1)
      std::cout << "set stopping criteria" << std::endl;
    //set stopping criteria
    if(maxit_opt[0])
      optimizer.set_maxeval(maxit_opt[0]);
    else
      optimizer.set_xtol_rel(tolerance_opt[0]);
    double minf=0;
    std::vector<double> x=init_opt;

    optimizer.optimize(x, minf);
    if(verbose_opt[0]){
      std::cout << "optimized with " << optimizer.get_algorithm_name() << "..." << std::endl;
      for(int index=0;index<x.size();++index){
        if(bcrad_opt[0])
          std::cout << setprecision(12) << " -bca " << theModel.deg2rad(x[index]);
        else
          std::cout << setprecision(12) << " -bca " << x[index];
      }
      std::cout << std::endl;
    }
  }

  double posX=0;
  double posY=0;
  ofstream outputStream;
  if(output_opt.size()){
    if(verbose_opt[0]>1)
      std::cout << "opening output file: " << output_opt[0] << std::endl;
    outputStream.open(output_opt[0].c_str());
  }
  vector<double> errorv;

  oSourceSRS.importFromEPSG(4326);
  oTargetSRS.importFromEPSG(t_srs_opt[0]);
  poCT = OGRCreateCoordinateTransformation( &oSourceSRS,
                                            &oTargetSRS );
  for(int ipoint=0;ipoint<vID.size();++ipoint){
    if(verbose_opt[0]>1)
      std::cout << "point: " << ipoint << std::endl;
    arma::vec pos_platform=theDataModel.getPosPlatform(ipoint);
    arma::vec att_platform=theDataModel.getAttPlatform(ipoint);

    arma::vec pos_calc;
    arma::vec pos_gcp;
    if(input_opt.size()||x_opt.size())
      pos_gcp=theDataModel.getPosGCP(ipoint);
    else
      pos_gcp=theDataModel.getPos(ipoint);
    if(errorZ_opt.size())
      pos_calc=theDataModel.getModel().getPos(pos_platform,att_platform,theDataModel.getRow(ipoint),theDataModel.getCol(ipoint),theDataModel.getHeight(ipoint)+errorZ_opt[0]);
    else
      pos_calc=theDataModel.getPos(ipoint);
    double e=theDataModel.getModel().getDistGeo(pos_gcp,pos_calc);
    if(theDataModel.getThreshold()&&e>theDataModel.getThreshold())
      continue;
    errorv.push_back(e);

    ostringstream gcpss;
    gcpss.precision(12);
    poCT->Transform(1,&(pos_gcp(0)),&(pos_gcp(1)),&(pos_gcp(2)));
    gcpss << pos_gcp(0) << " " << pos_gcp(1) << " " << pos_gcp(2) << " ";
    gcpss << errorv.back() << " ";

    poCT->Transform(1,&(pos_calc(0)),&(pos_calc(1)),&(pos_calc(2)));

    if(output_opt.size()){
      outputStream << setprecision(12) << vID[ipoint] << " " << theDataModel.getRow(ipoint) << " " << theDataModel.getCol(ipoint) << " " << pos_calc(0) << " " << pos_calc(1) << " " << pos_calc(2) << " " << gcpss.str();
      if(getzenith_opt[0])
        outputStream << " " << theDataModel.getModel().getZenith(att_platform,theDataModel.getRow(ipoint),theDataModel.getCol(ipoint));
      outputStream << std::endl;
    }
    else{
      std::cout << setprecision(12) << vID[ipoint] << " " << theDataModel.getRow(ipoint) << " " << theDataModel.getCol(ipoint) << " " << pos_calc(0) << " " << pos_calc(1) << " " << pos_calc(2) << " " << gcpss.str();
      if(getzenith_opt[0])
        std::cout << " " << theDataModel.getModel().getZenith(att_platform,theDataModel.getRow(ipoint),theDataModel.getCol(ipoint));
      std::cout << std::endl;
    }
  }
  if(verbose_opt[0]){
    if(output_opt.size())
      outputStream << "Number of GCP above threshold removed: " << nremoved << std::endl;
    else
      std::cout << "Number of GCP above threshold removed: " << nremoved << std::endl;
  }
  if(mean_opt[0]){
    double theMean=0;
    double theVar=0;
    statfactory::StatFactory stat;
    stat.meanVar(errorv,theMean,theVar);
    if(output_opt.size())
      outputStream << setprecision(12) << "mean stddev: " << theMean << " " << sqrt(theVar) << std::endl;
    else
      std::cout << setprecision(12) << "mean stdev: " << theMean << " " << sqrt(theVar) << std::endl;
  }
  if(output_opt.size())
    outputStream.close();
}      

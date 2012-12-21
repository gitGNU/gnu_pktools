/**********************************************************************
SensorModel.h: class for sensor model (calculate position on earth from pos and attidude information)
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
#ifndef _SENSORMODEL_
#define _SENSORMODEL_
#include <assert.h>
#include <math.h>
#include <vector>
#include <iostream>
// #include <gsl/gsl_matrix.h>
// #include <gsl/gsl_vector.h>
// #include <gsl/gsl_permutation.h>
// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_blas.h>
// #include <gsl/gsl_rng.h>
// #include <gsl/gsl_randist.h>
#include <gslwrap/vector_double.h>
#include <gslwrap/matrix_double.h>
#include <gslwrap/matrix_vector_operators.h>
#include "ogr_spatialref.h"

#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif

// #ifndef DEG2RAD
// #define DEG2RAD(DEG) (DEG/180.0*PI)
// #endif

namespace SensorModel
{
  enum Type { PUSHBROOM=0, WHISKBROOM=1 , FRAME=2};

  class SensorModel{
  public:
    SensorModel(void){
      m_bcpos.resize(3);
      m_bcatt.resize(3);
      m_bcpos[0]=0;
      m_bcpos[1]=0;
      m_bcpos[2]=0;
      m_bcatt[0]=0;
      m_bcatt[1]=0;
      m_bcatt[2]=0;
    };
    SensorModel(int theModel) : m_model(theModel){
      m_bcpos.resize(3);
      m_bcatt.resize(3);
      m_bcpos[0]=0;
      m_bcpos[1]=0;
      m_bcpos[2]=0;
      m_bcatt[0]=0;
      m_bcatt[1]=0;
      m_bcatt[2]=0;
    };
    ~SensorModel(){};
    // double deg2rad(double angle) const{return acos(-1)*angle/180.0;};
    static double deg2rad(double angle) {return PI*angle/180.0;};
    static double rad2deg(double angle) {return 180*angle/PI;};
    void setModel(int theModel){m_model=theModel;};
    int getModel() const {return m_model;};
    void setFOV(double fovDEG){m_fov=deg2rad(fovDEG);};
    void setNcol(int ncol){m_ncol=ncol;};
    void setNrow(int nrow){m_nrow=nrow;};
    void setDx(double dx){m_dx=dx;};
    void setDy(double dy){m_dy=dy;};
    void setF(double fc){m_fc=fc;};
    double getF() const {return m_fc;};
    void setPPx(double ppx){m_ppx=ppx;};
    void setPPy(double ppy){m_ppy=ppy;};
    void setBoresightPos(const gsl::vector& bcpos){m_bcpos=bcpos;};
    void setBoresightAtt(const gsl::vector& bcatt){m_bcatt=bcatt;};
    gsl::vector getBoresightPos() const{return m_bcpos;};
    gsl::vector getBoresightAtt() const{return m_bcatt;};
    void setPolynome(const std::vector<double>& polynome) {m_polynome=polynome;};
    void setDatum(const std::string& theDatum="WGS84"){m_spatialref.SetWellKnownGeogCS(theDatum.c_str());};
    double getZenith(const gsl::vector& att_platform, int row, int column) const{
      gsl::vector normallevel(3);
      gsl::vector normalplatform(3);

      gsl::vector apl_deg=att_platform;
      apl_deg+=m_bcatt;
      gsl::vector apl_rad(3);
      apl_rad[0]=deg2rad(apl_deg[0]);//roll
      apl_rad[1]=deg2rad(apl_deg[1]);//pitch
      apl_rad[2]=deg2rad(apl_deg[2]);//yaw

      if(getModel()==PUSHBROOM){
        gsl::vector scanAngle(2);
        scanAngle=scanAngle_PB(row,column);
        apl_rad[0]+=scanAngle[1];
        apl_rad[1]+=scanAngle[0];
      }
      else if(getModel()==WHISKBROOM){
        apl_rad[0]+=0;
        apl_rad[1]=tan(scanAngle_WB(column));
      }
      else//not implemented?
        assert(0);
      normallevel[0]=0;
      normallevel[1]=0;
      normallevel[2]=-1;
      gsl::matrix rotM(3,3);
      rotM=getRz(apl_rad[2])*getRy(apl_rad[1])*getRx(apl_rad[0]);
      normalplatform=rotM*normallevel;
      return rad2deg(acos(-normalplatform[2]));
    };
    gsl::vector getPos(const gsl::vector& pos_platform, const gsl::vector& att_platform, int row, int column, double elevation) const{
      gsl::vector thePosition(3);
      gsl::vector pos_ellips(3);
      gsl::vector ppl_deg=pos_platform;
      gsl::vector apl_deg=att_platform;
      ppl_deg+=m_bcpos;
      apl_deg+=m_bcatt;
      gsl::vector ppl_rad(3);
      gsl::vector apl_rad(3);
      ppl_rad[0]=deg2rad(ppl_deg[0]);
      ppl_rad[1]=deg2rad(ppl_deg[1]);
      ppl_rad[2]=ppl_deg[2];//add geoid elevation if necessary...
      apl_rad[0]=deg2rad(apl_deg[0]);//roll
      apl_rad[1]=deg2rad(apl_deg[1]);//pitch
      apl_rad[2]=deg2rad(apl_deg[2]);//yaw
      gsl::vector pos_ecef=pECEF(ppl_rad,apl_rad,row,column);
      pos_ellips=ecef2geo(pos_ecef);
      pos_ellips[2]=0;
      thePosition=getXYatZ(elevation,ppl_deg,pos_ellips);
      return thePosition;
    };
    double getDistGeo(const gsl::vector& pos1, const gsl::vector& pos2) const{
      double lon1=pos1[0];
      double lat1=pos1[1];
      double lon2=pos2[0];
      double lat2=pos2[1];
      double result;
      //simplified formula (spherical approximation)
      // result=2*asin(sqrt(pow(sin((lat1-lat2)/2),2) + cos(lat1)*cos(lat2)*pow(sin((lon1-lon2)/2),2)));
      //using loxodromes
      result=(pos1[1]==pos2[1]) ? lox2(deg2rad(pos1[1]),deg2rad(pos1[0]),deg2rad(pos2[0])) : lox1(deg2rad(pos1[1]),deg2rad(pos1[0]),deg2rad(pos2[1]),deg2rad(pos2[0]));
      return result;
    };
    gsl::vector geo2ecef(const gsl::vector& pos_geo) const{
      gsl::vector pos_ecef(3);
      double nu=getNu(pos_geo[1]);
      double f1=(nu+pos_geo[2])*cos(pos_geo[1]);
      double e1=getE1();
      pos_ecef[0]=f1*cos(pos_geo[0]);
      pos_ecef[1]=f1*sin(pos_geo[0]);
      pos_ecef[2]=(nu*(1-e1*e1)+pos_geo[2])*sin(pos_geo[1]);
      return pos_ecef;
    };
    gsl::vector ecef2geo(const gsl::vector& pos_ecef) const{
      gsl::vector pos_geo(3);
      pos_geo[0]=atan2(pos_ecef[1],pos_ecef[0]);
      while(pos_geo[0]>2*PI)
        pos_geo[0]-=2*PI;
      double r=sqrt(pos_ecef[0]*pos_ecef[0]+pos_ecef[1]*pos_ecef[1]);
      double f_earth=1.0/m_spatialref.GetInvFlattening();
      double e2=f_earth*(2-f_earth);
      double Ae=m_spatialref.GetSemiMajor();
      pos_geo[1]=atan(pos_ecef[2]/r);
      double c=1;
      double lat=1E+30;
      int iterations=0;
      while(sqrt((pos_geo[1]-lat)*(pos_geo[1]-lat))>1E-15){
        ++iterations;
        lat=pos_geo[1];
        c=1.0/sqrt(1-e2*sin(lat)*sin(lat));
        pos_geo[1]=atan((pos_ecef[2]+Ae*c*e2*sin(lat))/r);
      }
      pos_geo[2]=r/cos(pos_geo[1])-Ae*c;
      pos_geo[0]=rad2deg(pos_geo[0]);
      pos_geo[1]=rad2deg(pos_geo[1]);
      return pos_geo;
    };
  private:
    //line function; In this function, point 0 is the platform position, point 1 is the ellipsoid position
    gsl::vector getXYatZ (double inZ, gsl::vector p0, gsl::vector p1) const{
      gsl::vector posatz(3);
      posatz[0]=p0[0]+(p0[0]-p1[0])/(p0[2]-p1[2])*(inZ-p0[2]);
      posatz[1]=p0[1]+(p0[1]-p1[1])/(p0[2]-p1[2])*(inZ-p0[2]);
      posatz[2]=inZ;
      return posatz;
    };
    //get First eccentricity of the Earth's surface 
    double getE1() const{double f_earth=1.0/m_spatialref.GetInvFlattening();return sqrt(f_earth*(2-f_earth));};
    //get Second eccentricity of the Earth's surface 
    double getEearth() const{double Ae=m_spatialref.GetSemiMajor(); double Be=m_spatialref.GetSemiMinor();return sqrt(Ae*Ae-Be*Be)/Ae;};
    //get Radius of curvature normal to the meridian
    double getNu(double lat) const{double sinlat=sin(lat);return m_spatialref.GetSemiMajor()/sqrt(1-getE1()*getE1()*sinlat*sinlat);};
    //get Ellipsoid radius in the plane of the meridian
    double getRo(double lat) const{double e1=getE1();double sinlat=sin(lat);return (m_spatialref.GetSemiMajor()*(1-e1*e1))/sqrt(1-getE1()*getE1()*sinlat*sinlat);};
    double getEccentricity() const{return sqrt(1-m_spatialref.GetSemiMinor()*m_spatialref.GetSemiMinor()/m_spatialref.GetSemiMajor()/m_spatialref.GetSemiMajor());};
    //Loxodromes are closely related with Mercator projection as the main feature of this projection is that loxodromes are projected as straight lines. To find the angle of the loxodrome between two points one has simply to project the two points using a Mercator projection and calculate the angle of the line joining the two resulting points. Formulas for Mercator projection on the ellipsoid can be found in many places (for instance "Map Projections - A Working Manual" by J.P. Snyder, USGS Professional Paper 1395 pg. 38-47)
    //project lat lon to Mercator projection
    void getMxy(double& lon, double& lat) const{
      // lon=m_spatialref.GetSemiMajor()*lon;
      // double sinlat=sin(lat);
      // double f1=(1+sinlat)/(1-sinlat);
      // double eps=getEccentricity();
      // double f2=pow((1-eps*sinlat)/(1+eps*sinlat),eps);
      // lat=m_spatialref.GetSemiMajor()/2.0*log(f1*f2);
      // //test
      // std::cout << "Jan Mercator lat, lon: " << lat << ", " << lon << std::endl;
      OGRSpatialReference oSourceSRS, oTargetSRS;
      OGRCoordinateTransformation *poCT;
      oSourceSRS.importFromEPSG(4326);
      oTargetSRS.importFromEPSG(3395);
      poCT = OGRCreateCoordinateTransformation( &oSourceSRS,
                                                &oTargetSRS );
      lon=rad2deg(lon);
      lat=rad2deg(lat);
      poCT->Transform(1,&lon,&lat);
    };

    gsl::vector pECEF(const gsl::vector& pos, const gsl::vector& attitude, int row, int column) const{
      gsl::vector A(3);
      gsl::matrix B(3,3);
      B.set_element(0,0,-sin(pos[1])*cos(pos[0]));
      B.set_element(0,1,-sin(pos[0]));
      B.set_element(0,2,-cos(pos[1])*cos(pos[0]));
      B.set_element(1,0,-sin(pos[1])*sin(pos[0]));
      B.set_element(1,1,cos(pos[0]));
      B.set_element(1,2,-cos(pos[1])*sin(pos[0]));
      B.set_element(2,0,cos(pos[1]));
      B.set_element(2,1,0);
      B.set_element(2,2,-sin(pos[1]));
      A=geo2ecef(pos);
      gsl::vector C(3);
      A+=B*SV(row,column,attitude,pos[2]);
      return A;
    };
    gsl::vector getIV(int row, int column) const{
      gsl::vector iv(3);
      iv[0]=0;
      gsl::vector scanAngle(2);
      if(getModel()==PUSHBROOM){
        scanAngle=scanAngle_PB(row,column);
        iv[0]=scanAngle[1];
        iv[1]=tan(scanAngle[0]);
      }
      else if(getModel()==WHISKBROOM){
        iv[0]=0;
        iv[1]=tan(scanAngle_WB(column));
      }
      else//not implemented?
        assert(0);
      iv[2]=-1.0;
      return iv;
    };
    double scanAngle_WB(int column) const{
      return (-m_fov/2.0+column*(m_fov/(m_ncol-1)));
    };
    gsl::vector SV(int row, int column, const gsl::vector& attitude, double z) const{
      gsl::vector rv(3);
      gsl::matrix rotM(3,3);
      rotM=getRz(attitude[2])*getRy(attitude[1])*getRx(attitude[0]);
      rv=rotM*getIV(row,column);
      double height=rv[2];
      // gsl_blas_dscal(z/height,&rv);
      // return rv;
      rv*=z/height;
      return rv;
    };
    gsl::vector scanAngle_PB(int row, int column) const{
      gsl::vector alpha(2);
      double r=corr_along(column);
      double theAx=(column<m_ppx) ? m_dx*0.000001*(m_ppx-(column+corr_across(column))) : m_dx*0.000001*(m_ppx-(column-corr_across(column)));
      double theAy=m_dy*0.000001*(m_ppy-r);
      alpha[0]=(getModel()==FRAME) ? atan(theAx/getF()) : atan(-theAx/getF());
      alpha[1]=atan(-theAy/getF());
      return alpha;
    };
    //geocentric coordinate system is right-handed, orthogogal Cartesian system with its origin at the centre of the earth. The direct Helmert transformation from lat/lon and ellipsoid height to geocentric x,y,z is returned in posX, posY, posZ
    //Euclidian distance between points on sphere
    //Euclidian distance between two vectors
    double getDist(const gsl::vector& v1, const gsl::vector& v2) const{gsl::vector dv=v1;dv-=v2;return dv.norm2();};
    double lox1(double lat1, double lon1, double lat2, double lon2) const{
      double result=(M(lat2)-M(lat1))/cos(Az(lat1,lon1,lat2,lon2));
      return result;
    };
    double lox2(double lat, double lon1, double lon2) const{
      double eps2=getEccentricity();
      eps2*=eps2;
      double sin2=sin(lat);
      sin2*=sin2;
      return m_spatialref.GetSemiMajor()*(lon2-lon1)*cos(lat)/sqrt(1-eps2*sin2);
    };
    //length of the meridian from Equator to point
    double M(double lat) const{
      double eps2=getEccentricity();
      eps2*=eps2;
      double eps4=eps2*eps2;
      double eps6=eps4*eps2;
      double t1=(1-eps2/4.0-3*eps4/64.0-5*eps6/256.0)*lat;
      double t2=(3*eps2/8.0+3*eps4/32.0+45*eps6/1024.0)*sin(2*lat);
      double t3=(15*eps4/256.0+45*eps6/1024.0)*sin(4*lat);
      double t4=35*eps6/3072.0*sin(6*lat);
      double semiMajor=m_spatialref.GetSemiMajor();
      return semiMajor*(t1-t2+t3-t4);
    };
    double atan33(double x, double y) const{
      double angle=atan2(y,x);
      while(angle<0)
        angle+=2*PI;
      angle=PI/2.0+2*PI-angle;
      while(angle>2*PI)
        angle-=2*PI;
      return angle;
    };
    //The azimuth from North (clockwise) of the loxodrome from point 1 to point 2 will than be given by:
    double Az(double lat1, double lon1, double lat2, double lon2) const{
      double x1=lon1;
      double x2=lon2;
      double y1=lat1;
      double y2=lat2;
      getMxy(x1,y1);
      getMxy(x2,y2);
      return atan33(x2-x1,y2-y1);
    };
    double corr_along(double c) const{
      return 0;
    };
    double corr_across(double c) const{
      double result=0;
      if(m_polynome.size()){
        double absc=1;
        result+=absc*m_polynome[0];
        for(int degree=1;degree<m_polynome.size();++degree){
          absc*=sqrt((c-m_ppx)*(c-m_ppx));
          result+=absc*m_polynome[degree];
        }
      }
      return result;
    };
    gsl::matrix getRx(double theta) const{
      gsl::matrix Rx(3,3);
      Rx.set_element(0,0,1);
      Rx.set_element(0,1,0);
      Rx.set_element(0,2,0);
      Rx.set_element(1,0,0);
      Rx.set_element(1,1,cos(theta));
      Rx.set_element(1,2,-sin(theta));
      Rx.set_element(2,0,0);
      Rx.set_element(2,1,sin(theta));
      Rx.set_element(2,2,cos(theta));
      return Rx;
    };    
    gsl::matrix getRy(double theta) const{
      gsl::matrix Ry(3,3);
      Ry.set_element(0,0,cos(theta));
      Ry.set_element(0,1,0);
      Ry.set_element(0,2,sin(theta));
      Ry.set_element(1,0,0);
      Ry.set_element(1,1,1);
      Ry.set_element(1,2,0);
      Ry.set_element(2,0,-sin(theta));
      Ry.set_element(2,1,0);
      Ry.set_element(2,2,cos(theta));
      return Ry;
    };    
    gsl::matrix getRz(double theta) const{
      gsl::matrix Rz(3,3);
      Rz.set_element(0,0,cos(theta));
      Rz.set_element(0,1,-sin(theta));
      Rz.set_element(0,2,0);
      Rz.set_element(1,0,sin(theta));
      Rz.set_element(1,1,cos(theta));
      Rz.set_element(1,2,0);
      Rz.set_element(2,0,0);
      Rz.set_element(2,1,0);
      Rz.set_element(2,2,1);
      return Rz;
    };    
    int m_model;
    OGRSpatialReference m_spatialref;
    //Ae=m_spatialref.GetSemiMajor()
    //Be=m_spatialref.GetSemiMinor()
    //f_earth=1.0/m_spatialref.GetInvFlattening()
    std::string m_datum;
    double m_ppx;
    double m_ppy;
    double m_fc;
    double m_fov;
    int m_ncol;
    int m_nrow;
    double m_dx;
    double m_dy;
    gsl::vector m_bcpos;
    gsl::vector m_bcatt;
    std::vector<double> m_polynome;
  };
}
#endif //_SENSORMODEL_

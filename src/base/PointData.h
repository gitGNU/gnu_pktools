#ifndef _POINTDATA_
#define _POINTDATA_
#include <assert.h>
#include <math.h>
#include <vector>

struct State{
  double k;
  double e;
  double a;
  double haze;
};

class PointData{
 public: 
  PointData();//constructor
  ~PointData();//destructor
  void setCosFi(double fi){
    if(fi>180)
      fi-=360;
    else if(fi<-180)
      fi+=360;
    m_cosfi=cos(deg2rad(fi));
  };
  void setCosSZA(double sza){
    assert(sza>=0);
    m_cosza=cos(deg2rad(sza));
  };
  double getReflectance() const{return m_reflectance;};
  void setReflectance(double reflectance) {m_reflectance=reflectance;};
  double getImage() const{return m_image;}; 
  void setImage(int image) {m_image=image;}; 
  void setR(double r){m_r=r;};
  double getR() const{return m_r;};
  double correctReflectance(const State &state);
  double modelReflectance(const State &state);
  static double m_deltaT;
  static double m_tau;
  static double m_epsilon;
  static double m_flag;
  static double m_residual;
 private:
  double deg2rad(double angle) const{return acos(-1)*angle/180.0;};
  double m_reflectance;
  double m_cosza;
  double m_cosfi;
  double m_r;
  int m_image;
};

#endif //_POINTDATA_

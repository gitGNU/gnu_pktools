#include "PointData.h"

double PointData::m_tau=0;
double PointData::m_deltaT=0;
double PointData::m_epsilon=0.00001;
double PointData::m_flag=1000;
double PointData::m_residual=0;

PointData::PointData(){
}

PointData::~PointData(){
}

double PointData::correctReflectance(const State &state){
  double ecosfi=state.e*m_cosfi;
  double optical_depth=(m_cosza*m_cosza>m_epsilon*m_epsilon) ? exp(-m_tau/m_cosza) : 0;
  double hotspot=(ecosfi<1) ? (state.k+(1-state.k)*(1-state.e*m_cosza)/(1-ecosfi)) : 1;
  double radial=(1+state.a*m_r);
  double correctedReflectance;
  correctedReflectance=(m_reflectance-state.haze)/m_cosza/optical_depth/m_deltaT/hotspot/radial;
  return correctedReflectance;
  // return (m_reflectance-state.haze)/((m_cosza*exp(-m_tau/m_cosza)*m_deltaT*(state.k+(1-state.k)*(1-state.e*m_cosza)/(1-ecosfi))*(1+state.a*m_r)));
}

double PointData::modelReflectance(const State &state){
  double ecosfi=state.e*m_cosfi;
  double optical_depth=exp(-m_tau/m_cosza);
  double hotspot=state.k+(1-state.k)*(1-state.e*m_cosza)/(1-ecosfi);
  double radial=(1+state.a*m_r);
  return m_reflectance*m_cosza*optical_depth*m_deltaT*hotspot*radial+state.haze;
  // return (m_reflectance-state.haze)/((m_cosza*exp(-m_tau/m_cosza)*m_deltaT*(state.k+(1-state.k)*(1-state.e*m_cosza)/(1-ecosfi))*(1+state.a*m_r)));
}


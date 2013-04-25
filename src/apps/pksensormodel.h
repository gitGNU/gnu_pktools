/**********************************************************************
pksensormodel.h: program to calculate geometric position based on row (sensor), col (sensor), roll, pitch, yaw and lens coordinates
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
#ifndef _PKSENSORMODEL_H_
#define _PKSENSORMODEL_H_
#include <vector>
#include <gslwrap/matrix_double.h>
#include "models/SensorModel.h"

double objFunction(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data);

class DataModel{
   public:
  DataModel() : m_threshold(0){};
  DataModel(const SensorModel::SensorModel& theModel) : m_model(theModel), m_threshold(0){};
  ~DataModel(){};
  void setModel(const SensorModel::SensorModel& theModel){m_model=theModel;};
  const SensorModel::SensorModel& getModel() const {return m_model;};
  int getSize() const{return m_posGCP.size();};
  void setThreshold(double theThreshold){m_threshold=theThreshold;};
  double getThreshold(){return m_threshold;};
  int erase(int index){
    m_attPlatform.erase(m_attPlatform.begin()+index);
    m_posPlatform.erase(m_posPlatform.begin()+index);
    m_posGCP.erase(m_posGCP.begin()+index);
    m_row.erase(m_row.begin()+index);
    m_col.erase(m_col.begin()+index);
  };
  int pushAttPlatform(const gsl::vector& atp){m_attPlatform.push_back(atp); return m_attPlatform.size();};
  int pushPosPlatform(const gsl::vector& ppl){m_posPlatform.push_back(ppl); return m_posPlatform.size();};
  int pushPosGCP(const gsl::vector& pgcp){m_posGCP.push_back(pgcp); return m_posGCP.size();};
  int pushRow(int r){m_row.push_back(r); return m_row.size();};
  int pushCol(int c){m_col.push_back(c); return m_col.size();};
  gsl::vector getPosPlatform(int index) const{assert(index>=0);assert(index<m_posPlatform.size());return(m_posPlatform[index]);};
  gsl::vector getAttPlatform(int index) const{assert(index>=0);assert(index<m_attPlatform.size());return(m_attPlatform[index]);};
  gsl::vector getPosGCP(int index) const{assert(index>=0);assert(index<m_posGCP.size());return(m_posGCP[index]);};
  gsl::vector getPos(int index) const{
    assert(index>=0);
    assert(index<m_posPlatform.size());
    assert(index<m_attPlatform.size());
    assert(index<m_row.size());
    assert(index<m_col.size());
    assert(index<m_posGCP.size());
    return(m_model.getPos(m_posPlatform[index],m_attPlatform[index],m_row[index],m_col[index],m_posGCP[index][2]));
  };
  double getDistGeo(int index) const{assert(index>=0);assert(index<m_posGCP.size());return(m_model.getDistGeo(m_posGCP[index],getPos(index)));};
  int getRow(int index) const{assert(index>=0);assert(index<m_row.size());return(m_row[index]);};
  int getCol(int index) const{assert(index>=0);assert(index<m_col.size());return(m_col[index]);};
  double getHeight(int index) const{assert(index>=0);assert(index<m_posGCP.size());return(m_posGCP[index][2]);};
  void setBoresightAtt(const gsl::vector& bc_att){
    m_model.setBoresightAtt(bc_att);
    // for(int index=0;index<m_attPlatform.size();++index)
    //   m_attPlatform[index]+=bc_att;
  };
   private:
  SensorModel::SensorModel m_model;
  vector<gsl::vector> m_posPlatform;
  vector<gsl::vector> m_posGCP;
  vector<gsl::vector> m_attPlatform;
  vector<int> m_row;
  vector<int> m_col;
  double m_threshold;
};

#endif //_PKSENSORMODEL_H_

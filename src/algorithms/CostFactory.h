/**********************************************************************
CostFactory.h: select features, typical use: feature selection for classification
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
#ifndef _COSTFACTORY_H_
#define _COSTFACTORY_H_

#include <math.h>
#include <vector>
#include <map>
#include "ConfusionMatrix.h"
#include "base/Vector2d.h"


class CostFactory{
public:
  CostFactory(void){};
  CostFactory(unsigned short cv, short verbose) : m_cv(cv), m_verbose(verbose){};

  virtual ~CostFactory(void){};
  void setCv(unsigned short cv){m_cv=cv;};
  void setClassValueMap(const std::string& classname, short classvalue){ m_classValueMap[classname]=classvalue;};
  std::map<std::string,short> getClassValueMap(){return m_classValueMap;};
  std::vector<std::string> getNameVector(){return m_nameVector;};
  void setNameVector(std::vector<std::string>& nameVector){m_nameVector=nameVector;};
  int getClassIndex(std::string classname) const {return m_cm.getClassIndex(classname);};
  //pushBackClassName is for confusion matrix
  void pushBackClassName(std::string classname){m_cm.pushBackClassName(classname,true);};//doSort=true
  //pushBackName is for nameVector in CostFactory
  void pushBackName(std::string classname){m_nameVector.push_back(classname);};
  void setNcTraining(const std::vector<unsigned int> nctraining){m_nctraining=nctraining;};
  void setNcTest(const std::vector<unsigned int> nctest){m_nctest=nctest;};
  //getCost needs to be implemented case by case (e.g., SVM, ANN)
  virtual double getCost(const std::vector<Vector2d<float> > &trainingFeatures)=0;
  
protected:
  confusionmatrix::ConfusionMatrix m_cm;
  std::map<std::string,short> m_classValueMap;
  std::vector<std::string> m_nameVector;
  std::vector<unsigned int> m_nctraining;
  std::vector<unsigned int> m_nctest;
  unsigned short m_cv;
  /* std::string m_classname; */
  short m_classvalue;
  short m_verbose;
};
#endif

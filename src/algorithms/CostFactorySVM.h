/**********************************************************************
CostFactorySVM.h: select features, typical use: feature selection for classification
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
#ifndef _COSTFACTORYSVM_H_
#define _COSTFACTORYSVM_H_

#include <math.h>
#include <vector>
#include <map>
#include "base/Vector2d.h"
#include "CostFactory.h"

namespace svm{
  enum SVM_TYPE {C_SVC=0, nu_SVC=1,one_class=2, epsilon_SVR=3, nu_SVR=4};
  enum KERNEL_TYPE {linear=0,polynomial=1,radial=2,sigmoid=3};
}

class CostFactorySVM : public CostFactory
{
public:
CostFactorySVM();
CostFactorySVM(std::string svm_type, std::string kernel_type, unsigned short kernel_degree, float gamma, float coef0, float ccost, float nu,  float epsilon_loss, int cache, float epsilon_tol, bool shrinking, bool prob_est, unsigned short cv, short verbose);
~CostFactorySVM();
double getCost(const std::vector<Vector2d<float> > &trainingFeatures);
  
private:
std::string m_svm_type;
std::string m_kernel_type;
unsigned short m_kernel_degree;
float m_gamma;
float m_coef0;
float m_ccost;
float m_nu;
float m_epsilon_loss;
int m_cache;
float m_epsilon_tol;
bool m_shrinking;
bool m_prob_est;
};
#endif

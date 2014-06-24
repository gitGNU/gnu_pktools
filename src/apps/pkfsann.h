/**********************************************************************
pkfsann.h: feature selection for ann classifier
Copyright (C) 2008-2014 Pieter Kempeneers

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
#include <string>
#include <vector>
#include "base/Vector2d.h"

#ifndef _PKFSANNH_H_
#define _PKFSANNH_H_

enum SelectorValue  { NA=0, SFFS=1, SFS=2, SBS=3, BFS=4};

class CostFactoryANN : public CostFactory
{
 public:
  CostFactoryANN();
  CostFactoryANN(const std::vector<unsigned int>& nneuron, float connection, const std::vector<float> weights, float learning, unsigned int maxit, unsigned short cv, bool verbose);
  ~CostFactoryANN();
  double getCost(const std::vector<Vector2d<float> > &trainingFeatures);
  
 private:
  std::vector<unsigned int> m_nneuron;
  float m_connection;
  const std::vector<float> m_weights;
  float m_learning;
  unsigned int m_maxit;
};


#endif

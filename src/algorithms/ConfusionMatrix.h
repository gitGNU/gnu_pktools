/**********************************************************************
ConfusionMatrix.h: class for (classification accuracy) confusion matrix
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
#ifndef _CONFUSIONMATRIX_H_
#define _CONFUSIONMATRIX_H_

#include <sstream>
#include <vector>
#include "base/Vector2d.h"
#include "base/Optionpk.h"

class ConfusionMatrix{
public:
  ConfusionMatrix();
  ConfusionMatrix(short nclass);
  ConfusionMatrix(const std::vector<std::string>& classNames);
  ConfusionMatrix(const ConfusionMatrix& cm);
  ConfusionMatrix& operator=(const ConfusionMatrix& cm);
  short size() const {return m_results.size();};
  void resize(short nclass);
  void setClassNames(const std::vector<std::string>& classNames, bool doSort=false);
  void pushBackClassName(const std::string& className, bool doSort=false);
  void setResults(const Vector2d<double>& theResults);
  void setResult(const std::string& theRef, const std::string& theClass, double theResult);
  void incrementResult(const std::string& theRef, const std::string& theClass, double theIncrement);
  void clearResults();
  double nReference(const std::string& theRef) const;
  double nReference() const;
  double nClassified(const std::string& theRef) const;
  int nClasses() const {return m_classes.size();};
  std::string getClass(int iclass) const {assert(iclass>=0);assert(iclass<m_classes.size());return m_classes[iclass];};
  int getClassIndex(std::string className) const {
    int index=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),className));
    assert(index>=0);
    if(index<m_results.size())
      return(index);
    else
      return(-1);
  }
  std::vector<std::string> getClassNames() const {return m_classes;};
  ~ConfusionMatrix();
  double pa(const std::string& theClass, double* se95=NULL) const;
  double ua(const std::string& theClass, double* se95=NULL) const;
  double oa(double* se95=NULL) const;
  int pa_pct(const std::string& theClass, double* se95=NULL) const;
  int ua_pct(const std::string& theClass, double* se95=NULL) const;
  int oa_pct(double* se95=NULL) const;
  double kappa() const;
  ConfusionMatrix& operator*=(double weight);
  ConfusionMatrix operator*(double weight);
  ConfusionMatrix& operator+=(const ConfusionMatrix &cm);
  ConfusionMatrix operator+(const ConfusionMatrix &cm){
    return ConfusionMatrix(*this)+=cm;
  }
  void sortClassNames();
  friend ostream& operator<<(ostream& os, const ConfusionMatrix &cm){
    for(int iclass=0;iclass<cm.nClasses();++iclass)
      os << "\t" << cm.m_classes[iclass];
    os << endl;
    assert(cm.m_classes.size()==cm.m_results.size());
    for(int irow=0;irow<cm.m_results.size();++irow){
      os << cm.m_classes[irow];
      for(int icol=0;icol<cm.m_results[irow].size();++icol)
        os << "\t" << cm.m_results[irow][icol];
      os << endl;
    }
    return os;
  };
private:
  std::vector<std::string> m_classes;
  Vector2d<double> m_results;
};

#endif /* _CONFUSIONMATRIX_H_ */

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

namespace confusionmatrix
{
  enum CM_FORMAT { ASCII = 0, LATEX = 1, HTML = 2 };

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
    unsigned int index=0;
    for(index=0;index<m_classes.size();++index){
      if(m_classes[index]==className)
	break;
    }
    if(index>=m_classes.size())
      index=-1;
    return index;
    //    int index=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),className));
    //    assert(index>=0);
    //    if(index<m_results.size())
    //      return(index);
    //    else
    //      return(-1);
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

  void reportSE95(bool doReport) {m_se95=doReport;};
  void setFormat(const CM_FORMAT& theFormat) {m_format=theFormat;};
  void setFormat(const std::string theFormat) {m_format=getFormat(theFormat);};
  CM_FORMAT getFormat() const {return m_format;};
 
  static const CM_FORMAT getFormat(const std::string theFormat){
    if(theFormat=="ascii") return(ASCII);
    else if(theFormat=="latex") return(LATEX);
    else{
      std::string errorString="Format not supported: ";
      errorString+=theFormat;
      errorString+=" use ascii or latex";
      throw(errorString);
    }
  };

  friend std::ostream& operator<<(std::ostream& os, const ConfusionMatrix &cm){
    std::ostringstream streamLine;
    /* streamosclass << iclass; */
    /* m_classes[iclass]=osclass.str(); */

    std::string fieldSeparator=" ";
    std::string lineSeparator="";
    std::string mathMode="";
    switch(cm.getFormat()){
    case(LATEX):
      fieldSeparator=" & ";
      lineSeparator="\\\\";
      mathMode="$";
      break;
    case(ASCII):
    default:
      fieldSeparator="\t";
      lineSeparator="";
      mathMode="";
      break;
    }

    double se95_ua=0;
    double se95_pa=0;
    double se95_oa=0;
    double dua=0;
    double dpa=0;
    // double doa=0;
    // doa = cm.oa(&se95_oa);

    if(cm.getFormat()==LATEX){
      os << "\\documentclass{article}" << std::endl;
      os << "\\begin{document}" << std::endl;
    }
    os << "Kappa = " << mathMode << cm.kappa() << mathMode ;
    os << ", Overall Acc. = " << mathMode << 100.0*cm.oa() << mathMode ;
    if(cm.m_se95)
      os << " (" << mathMode << se95_oa << mathMode << ")";
    os << std::endl;
    os << std::endl;
    if(cm.getFormat()==LATEX){
      os << "\\begin{tabular}{@{}l";
      for(int iclass=0;iclass<cm.nClasses();++iclass)
	os << "l";
      os << "}" << std::endl;
      os << "\\hline" << std::endl;
    }
    
    os << "Class";
    for(int iclass=0;iclass<cm.nClasses();++iclass)
      os << fieldSeparator << cm.m_classes[iclass];
    os << lineSeparator << std::endl;
    if(cm.getFormat()==LATEX)
      os << "\\hline" << std::endl;
    assert(cm.m_classes.size()==cm.m_results.size());
    for(unsigned int irow=0;irow<cm.m_results.size();++irow){
      os << cm.m_classes[irow];
      for(unsigned int icol=0;icol<cm.m_results[irow].size();++icol)
        os << fieldSeparator << cm.m_results[irow][icol];
      os << lineSeparator<< std::endl;
    }
    if(cm.getFormat()==LATEX){
      os << "\\hline" << std::endl;
    }
    else
      os << std::endl;

    os << "User' Acc.";
    for(int iclass=0;iclass<cm.nClasses();++iclass){
      dua=cm.ua_pct(cm.m_classes[iclass],&se95_ua);
      os << fieldSeparator << dua;
      if(cm.m_se95)
	os << " (" << se95_ua << ")";
    }
    os << lineSeparator<< std::endl;
    os << "Prod. Acc.";
    for(int iclass=0;iclass<cm.nClasses();++iclass){
      dpa=cm.pa_pct(cm.m_classes[iclass],&se95_ua);
      os << fieldSeparator << dpa;
      if(cm.m_se95)
	os << " (" << se95_pa << ")";
    }
    os << lineSeparator<< std::endl;
    if(cm.getFormat()==LATEX){
      os << "\\end{tabular}" << std::endl;
      os << "\\end{document}" << std::endl;
    }
    return os;
  };
private:
  std::vector<std::string> m_classes;
  Vector2d<double> m_results;
  CM_FORMAT m_format;
  bool m_se95;
};
}
#endif /* _CONFUSIONMATRIX_H_ */

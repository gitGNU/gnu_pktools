/**********************************************************************
ConfusionMatrix.cc: class for (classification accuracy) confusion matrix
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
#include "ConfusionMatrix.h"
#include <iostream>
#include <numeric>

using namespace confusionmatrix;

bool compareClass(const std::string& string1, const std::string& string2){
  int int1=string2type<int>(string1);
  int int2=string2type<int>(string2);
  return(int1<int2);
};

ConfusionMatrix::ConfusionMatrix()
  : m_classes(),m_results(),m_se95(true),m_format(ASCII)
{
}

ConfusionMatrix::~ConfusionMatrix()
{
}

//constructor where class names are 0,1,...,nclass-1
ConfusionMatrix::ConfusionMatrix(short nclass){
  resize(nclass);
}

ConfusionMatrix::ConfusionMatrix(const std::vector<std::string>& classNames){
  setClassNames(classNames);
}

//copy constructor
ConfusionMatrix::ConfusionMatrix(const ConfusionMatrix& cm){
  setClassNames(cm.m_classes);
  setResults(cm.m_results);
}

//assignment operator
ConfusionMatrix& ConfusionMatrix::operator=(const ConfusionMatrix& cm){
  //check for self-assignment by comparing the address of the implicit object and parameter
  if(this==&cm)
    return *this;
  else{
    setClassNames(cm.m_classes);
    setResults(cm.m_results);
  }
  return *this;
}

ConfusionMatrix& ConfusionMatrix::operator+=(const ConfusionMatrix &cm)
{
  if(cm.m_classes.size()!=this->m_classes.size()){
    std::cerr << "error0: "<< cm.m_classes.size() << "!=" << this->m_classes.size() << std::endl;
    exit(0);
  }
  if(cm.m_results.size()!=this->m_results.size()){
    std::cerr << "error1: "<< cm.m_results.size() << "!=" << this->m_results.size() << std::endl;
    exit(1);
  }
  for(int irow=0;irow<m_results.size();++irow){
    if(cm.m_results[irow].size()!=this->m_results[irow].size()){
      std::cerr << "error2: " << cm.m_results[irow].size() << "!=" << this->m_results[irow].size() << std::endl;
      exit(2);
    }
    for(int icol=0;icol<m_results[irow].size();++icol)
      this->m_results[irow][icol]+=cm.m_results[irow][icol];
  }
  return *this;
}

ConfusionMatrix& ConfusionMatrix::operator*=(double weight)
{
  for(int irow=0;irow<m_results.size();++irow){
    for(int icol=0;icol<m_results[irow].size();++icol)
      m_results[irow][icol]*=weight;
  }
  return *this;
}

void ConfusionMatrix::sortClassNames(){
  sort(m_classes.begin(),m_classes.end(),compareClass);
}

ConfusionMatrix ConfusionMatrix::operator*(double weight)
{
  ConfusionMatrix result = *this;//make a copy of myself
  result*=weight;
  return result;
}

void ConfusionMatrix::resize(short nclass){
  m_classes.resize(nclass);
  for(short iclass=0;iclass<nclass;++iclass){
    std::ostringstream osclass;
    osclass << iclass;
    m_classes[iclass]=osclass.str();
  }
  m_results.resize(nclass,nclass);
}

void ConfusionMatrix::setClassNames(const std::vector<std::string>& classNames, bool doSort){
  m_classes=classNames;
  if(doSort)
    sortClassNames();
  if(m_results.size()!=m_classes.size())
    m_results.resize(m_classes.size(),m_classes.size());
}

void ConfusionMatrix::pushBackClassName(const std::string& className, bool doSort){
  m_classes.push_back(className);
  if(doSort)
    sortClassNames();
  if(m_results.size()!=m_classes.size())
    m_results.resize(m_classes.size(),m_classes.size());
}  


void ConfusionMatrix::setResults(const Vector2d<double>& theResults){
  m_results=theResults;
}

void ConfusionMatrix::clearResults(){
  m_results.clear();
  m_results.resize(m_classes.size(),m_classes.size());
}

void ConfusionMatrix::setResult(const std::string& theRef, const std::string& theClass, double theResult){
  // int ir=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theRef));
  // int ic=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theClass));
  // assert(ir>=0);
  // assert(ir<m_results.size());
  // assert(ic>=0);
  // assert(ic<m_results[ir].size());
  int ir=getClassIndex(theRef);
  int ic=getClassIndex(theClass);
  m_results[ir][ic]=theResult;
}

void ConfusionMatrix::incrementResult(const std::string& theRef, const std::string& theClass, double theIncrement){
  // int ir=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theRef));
  // int ic=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theClass));
  int ir=getClassIndex(theRef);
  int ic=getClassIndex(theClass);
  assert(ir>=0);
  if(ir>=m_results.size())
    std::cerr << "Error: " << theRef << " not found in class ConfusionMatrix when incrementing for class " << theClass << std::endl;
  assert(ir<m_results.size());
  assert(ic>=0);
  assert(ic<m_results[ir].size());
  m_results[ir][ic]+=theIncrement;
}

double ConfusionMatrix::nReference(const std::string& theRef) const{
  // int ir=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theRef));
  int ir=getClassIndex(theRef);
  return accumulate(m_results[ir].begin(),m_results[ir].end(),0);
}

double ConfusionMatrix::nReference() const{
  double nref=0;
  for(int ir=0;ir<m_classes.size();++ir)
    nref+=accumulate(m_results[ir].begin(),m_results[ir].end(),0);
  return nref;
}

double ConfusionMatrix::nClassified(const std::string& theClass) const{
  // int ic=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theClass));
  int ic=getClassIndex(theClass);
  double nclassified=0;
  for(int iref=0;iref<m_results.size();++iref){
    assert(ic<m_results[iref].size());
    nclassified+=m_results[iref][ic];
  }
  return(nclassified);
}

double ConfusionMatrix::pa(const std::string& theClass, double* se95) const{
  assert(m_results.size());
  assert(m_results.size()==m_classes.size());
  double producer=0;
  // int ir=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theClass));
  int ir=getClassIndex(theClass);
  assert(ir>=0);
  assert(ir<m_results.size());
  assert(!theClass.compare(m_classes[ir]));
  for(int iclass=0;iclass<m_results.size();++iclass){
    assert(iclass<m_results[ir].size());
    producer+=m_results[ir][iclass];
  }
  double dpa=(producer>0)? static_cast<double>(m_results[ir][ir])/producer : 0;
  double dqa=1.0-dpa;
  if(se95!=NULL)
    *se95=(dpa<1&&dpa>0)? sqrt(dpa*dqa/(producer-1)) : 0;
  return dpa;
}

int ConfusionMatrix::pa_pct(const std::string& theClass, double* se95) const{
  double dpa=pa(theClass,se95);
  if(se95!=NULL)
    *se95=static_cast<double>(static_cast<int>(0.5+1000*(*se95)))/10.0;
  return static_cast<int>(0.5+100.0*dpa);
}
    

double ConfusionMatrix::ua(const std::string& theClass, double* se95) const{
  assert(m_results.size());
  assert(m_results.size()==m_classes.size());
  double user=0;
  // int ic=distance(m_classes.begin(),find(m_classes.begin(),m_classes.end(),theClass));
  int ic=getClassIndex(theClass);
  assert(ic>=0);
  assert(ic<m_results.size());
  assert(!theClass.compare(m_classes[ic]));
  for(int iref=0;iref<m_results.size();++iref){
    assert(ic<m_results[iref].size());
    user+=m_results[iref][ic];
  }
  double dua=(user>0)? static_cast<double>(m_results[ic][ic])/user : 0;
  double dva=1.0-dva;
  if(se95!=NULL)
    *se95=(dua<1&&dua>0)? sqrt(dua*dva/(user-1)) : 0;
  return dua;
}

int ConfusionMatrix::ua_pct(const std::string& theClass,double* se95) const{
  double dua=ua(theClass,se95);
  if(se95!=NULL)
    *se95=static_cast<double>(static_cast<int>(0.5+1000*(*se95)))/10.0;
  return static_cast<int>(0.5+100.0*dua);
}

double ConfusionMatrix::oa(double* se95) const{
  double ntotal=m_results.sum();
  double pCorrect=0;
  for(int iclass=0;iclass<m_classes.size();++iclass)
    pCorrect+=static_cast<double>(m_results[iclass][iclass])/ntotal;
  double qCorrect=1-pCorrect;
  if(se95!=NULL)
    *se95=(pCorrect<1&&pCorrect>0)? sqrt(pCorrect*qCorrect/(ntotal-1)) : 0;
  if(ntotal>0)
    return(pCorrect);
  else
    return(0);
}

int ConfusionMatrix::oa_pct(double* se95) const{
  double doa=oa(se95);
  if(se95!=NULL)
    *se95=static_cast<double>(static_cast<int>(0.5+1000*(*se95)))/10.0;
  return static_cast<int>(0.5+100.0*doa);
}

double ConfusionMatrix::kappa() const{
  double ntotal=m_results.sum();
  double pChance=0;
  double pCorrect=0;
  for(int iclass=0;iclass<m_classes.size();++iclass){
    pChance+=nClassified(m_classes[iclass])*nReference(m_classes[iclass])/ntotal/ntotal;
    pCorrect+=static_cast<double>(m_results[iclass][iclass])/ntotal;
  }
  if(pChance<1)
    return((pCorrect-pChance)/(1-pChance));
  else
    return(0);
}

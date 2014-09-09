/**********************************************************************
CostFactorySVM.cc: select features, typical use: feature selection for classification
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
#include "CostFactorySVM.h"
#include "svm.h"

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))

CostFactorySVM::CostFactorySVM()
    : CostFactory(2,0), m_svm_type("C_SVC"), m_kernel_type("radial"), m_kernel_degree(3), m_gamma(1.0), m_coef0(0), m_ccost(1000), m_nu(0.5),  m_epsilon_loss(100), m_cache(100), m_epsilon_tol(0.001), m_shrinking(false), m_prob_est(true){
}

CostFactorySVM::~CostFactorySVM(){
}

CostFactorySVM::CostFactorySVM(std::string svm_type, std::string kernel_type, unsigned short kernel_degree, float gamma, float coef0, float ccost, float nu,  float epsilon_loss, int cache, float epsilon_tol, bool shrinking, bool prob_est, unsigned short cv, short verbose)
    : CostFactory(cv,verbose), m_svm_type(svm_type), m_kernel_type(kernel_type), m_kernel_degree(kernel_degree), m_gamma(gamma), m_coef0(coef0), m_ccost(ccost), m_nu(nu),  m_epsilon_loss(epsilon_loss), m_cache(cache), m_epsilon_tol(epsilon_tol), m_shrinking(shrinking), m_prob_est(prob_est){};

double CostFactorySVM::getCost(const std::vector<Vector2d<float> > &trainingFeatures){
  std::map<std::string, svm::SVM_TYPE> svmMap;

  svmMap["C_SVC"]=svm::C_SVC;
  svmMap["nu_SVC"]=svm::nu_SVC;
  svmMap["one_class"]=svm::one_class;
  svmMap["epsilon_SVR"]=svm::epsilon_SVR;
  svmMap["nu_SVR"]=svm::nu_SVR;

  std::map<std::string, svm::KERNEL_TYPE> kernelMap;

  kernelMap["linear"]=svm::linear;
  kernelMap["polynomial"]=svm::polynomial;
  kernelMap["radial"]=svm::radial;
  kernelMap["sigmoid;"]=svm::sigmoid;

  unsigned short nclass=trainingFeatures.size();
  unsigned int ntraining=0;
  unsigned int ntest=0;
  for(int iclass=0;iclass<nclass;++iclass){
    ntraining+=m_nctraining[iclass];
    ntest+=m_nctest[iclass];
  }
  if(ntest)
    assert(!m_cv);
  if(!m_cv)
    assert(ntest);
  unsigned short nFeatures=trainingFeatures[0][0].size();

  struct svm_parameter param;
  param.svm_type = svmMap[m_svm_type];
  param.kernel_type = kernelMap[m_kernel_type];
  param.degree = m_kernel_degree;
  param.gamma = (m_gamma>0)? m_gamma : 1.0/nFeatures;
  param.coef0 = m_coef0;
  param.nu = m_nu;
  param.cache_size = m_cache;
  param.C = m_ccost;
  param.eps = m_epsilon_tol;
  param.p = m_epsilon_loss;
  param.shrinking = (m_shrinking)? 1 : 0;
  param.probability = (m_prob_est)? 1 : 0;
  param.nr_weight = 0;//not used: I use priors and balancing
  param.weight_label = NULL;
  param.weight = NULL;
  param.verbose=(m_verbose>1)? true:false;
  struct svm_model* svm;
  struct svm_problem prob;
  struct svm_node* x_space;

  prob.l=ntraining;
  prob.y = Malloc(double,prob.l);
  prob.x = Malloc(struct svm_node *,prob.l);
  x_space = Malloc(struct svm_node,(nFeatures+1)*ntraining);
  unsigned long int spaceIndex=0;
  int lIndex=0;
  for(int iclass=0;iclass<nclass;++iclass){
    // for(int isample=0;isample<trainingFeatures[iclass].size();++isample){
    for(int isample=0;isample<m_nctraining[iclass];++isample){
      prob.x[lIndex]=&(x_space[spaceIndex]);
      for(int ifeature=0;ifeature<nFeatures;++ifeature){
        x_space[spaceIndex].index=ifeature+1;
        x_space[spaceIndex].value=trainingFeatures[iclass][isample][ifeature];
        ++spaceIndex;
      }
      x_space[spaceIndex++].index=-1;
      prob.y[lIndex]=iclass;
      ++lIndex;
    }
  }

  assert(lIndex==prob.l);
  if(m_verbose>2)
    std::cout << "checking parameters" << std::endl;
  svm_check_parameter(&prob,&param);
  if(m_verbose>2)
    std::cout << "parameters ok, training" << std::endl;
  svm=svm_train(&prob,&param);
  if(m_verbose>2)
    std::cout << "SVM is now trained" << std::endl;

  m_cm.clearResults();
  if(m_cv>1){
    double *target = Malloc(double,prob.l);
    svm_cross_validation(&prob,&param,m_cv,target);
    assert(param.svm_type != EPSILON_SVR&&param.svm_type != NU_SVR);//only for regression
    for(int i=0;i<prob.l;i++){
      std::string refClassName=m_nameVector[prob.y[i]];
      std::string className=m_nameVector[target[i]];
      if(m_classValueMap.size())
	m_cm.incrementResult(type2string<short>(m_classValueMap[refClassName]),type2string<short>(m_classValueMap[className]),1.0);
      else
	m_cm.incrementResult(m_cm.getClass(prob.y[i]),m_cm.getClass(target[i]),1.0);
    }
    free(target);
  }
  else{
    struct svm_node *x_test;
    std::vector<double> result(nclass);
    x_test = Malloc(struct svm_node,(nFeatures+1));
    for(int iclass=0;iclass<nclass;++iclass){
      for(int isample=0;isample<m_nctest[iclass];++isample){
	for(int ifeature=0;ifeature<nFeatures;++ifeature){
	  x_test[ifeature].index=ifeature+1;
	  x_test[ifeature].value=trainingFeatures[iclass][m_nctraining[iclass]+isample][ifeature];
	}
	x_test[nFeatures].index=-1;
	double predict_label=0;
	assert(svm_check_probability_model(svm));
	predict_label = svm_predict_probability(svm,x_test,&(result[0]));
	// predict_label = svm_predict(svm,x_test);
	std::string refClassName=m_nameVector[iclass];
	std::string className=m_nameVector[static_cast<short>(predict_label)];
	if(m_classValueMap.size())
	  m_cm.incrementResult(type2string<short>(m_classValueMap[refClassName]),type2string<short>(m_classValueMap[className]),1.0);
	else
	  m_cm.incrementResult(refClassName,className,1.0);
      }
    }
    free(x_test);
  }
  if(m_verbose>1)
    std::cout << m_cm << std::endl;
  assert(m_cm.nReference());
  // if(m_verbose)

  // std::cout << m_cm << std::endl;
  // std::cout << "Kappa: " << m_cm.kappa() << std::endl;
  // double se95_oa=0;
  // double doa=0;
  // doa=m_cm.oa_pct(&se95_oa);
  // std::cout << "Overall Accuracy: " << doa << " (" << se95_oa << ")"  << std::endl;

  // *NOTE* Because svm_model contains pointers to svm_problem, you can
  // not free the memory used by svm_problem if you are still using the
  // svm_model produced by svm_train(). 
  // however, we will re-train the svm later on after the feature selection
  free(prob.y);
  free(prob.x);
  free(x_space);
  svm_free_and_destroy_model(&(svm));

  return(m_cm.kappa());
}

/**********************************************************************
FeatureSelector.h: select features, typical use: feature selection for classification
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
#ifndef _FEATURESELECTOR_H_
#define _FEATURESELECTOR_H_

#include <math.h>
#include <vector>
#include <list>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "ConfusionMatrix.h"
#include "base/IndexValue.h"
#include "base/Vector2d.h"
#include "gsl/gsl_combination.h"
#include "CostFactory.h"

class FeatureSelector
{
 public:
  FeatureSelector(){};
  ~FeatureSelector(){};
  template<class T> double forward(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures, short verbose=0);
  template<class T> double backward(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int minFeatures, short verbose=0);
  template<class T> double floating(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures=0, double epsilon=0.001, short verbose=0);
  template<class T> double bruteForce(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures=0, short verbose=0);
  
  private:
    template<class T> double addFeature(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, short verbose=0);
  template<class T> double removeFeature(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int& r, short verbose=0);  
  template<class T> double forwardUnivariate(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures, short verbose=0);
  /* template<class T> double forward(std::vector< Vector2d<T> >& v, double (*getCost)(const std::vector< Vector2d<T> >&), std::list<int>& subset, int maxFeatures, short verbose=0); */
  /* template<class T> double backward(std::vector< Vector2d<T> >& v, double (*getCost)(const std::vector< Vector2d<T> >&), std::list<int>& subset, int minFeatures, short verbose=0); */
  /* template<class T> double floating(std::vector< Vector2d<T> >& v, double (*getCost)(const std::vector< Vector2d<T> >&), std::list<int>& subset, int maxFeatures=0, double epsilon=0.001, short verbose=0); */
  /* template<class T> double bruteForce(std::vector< Vector2d<T> >& v, double (*getCost)(const std::vector< Vector2d<T> >&), std::list<int>& subset, int maxFeatures=0, short verbose=0); */
  
  /* private: */
  /*   template<class T> double addFeature(std::vector< Vector2d<T> >& v, double (*getCost)(const std::vector< Vector2d<T> >&), std::list<int>& subset, short verbose=0); */
  /* template<class T> double removeFeature(std::vector< Vector2d<T> >& v, double (*getCost)(const std::vector< Vector2d<T> >&), std::list<int>& subset, int& r, short verbose=0);   */
};

//sequential forward selection Univariate (N single best features)
template<class T> double FeatureSelector::forwardUnivariate(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures=0, short verbose){
  int maxLevels=v[0][0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
  std::vector<IndexValue> cost(maxLevels);
  std::list<int> tmpset=subset;//temporary set of selected features (levels)
  std::vector< Vector2d<T> > tmp(v.size());
  for(int ilevel=0;ilevel<maxLevels;++ilevel){
    if(find(tmpset.begin(),tmpset.end(),ilevel)==tmpset.end()){
      tmpset.push_back(ilevel);
      for(int iclass=0;iclass<v.size();++iclass){
// 	tmp[iclass].resize(v[iclass].size());
	v[iclass].selectCols(tmpset,tmp[iclass]);
      }
      try{
	IndexValue pv;
	pv.position=ilevel;
	pv.value=theCostFactory.getCost(tmp);
	cost[ilevel]=pv;
      }
      catch(...){
	IndexValue pv;
	pv.position=ilevel;
	pv.value=-1;
	cost[ilevel]=pv;
      }
      tmpset.pop_back();
    }
  }
  sort(cost.begin(),cost.end(),Compare_IndexValue());//decreasing order
  int ilevel=0;
  while((subset.size()<maxFeatures)&&(ilevel<maxLevels)){
    if(cost[ilevel].value>0)
      subset.push_back(cost[ilevel].position);
    if(verbose)
      std::cout << "feature " << subset.back() << " has cost: " << cost[ilevel].value << std::endl;
    ++ilevel;
  }
  double maxCost=-1;
  while(subset.size()){
    for(int iclass=0;iclass<v.size();++iclass){
//       tmp[iclass].resize(v[iclass].size());
      v[iclass].selectCols(subset,tmp[iclass]);
    }
    try{
      maxCost=theCostFactory.getCost(tmp);
    }
    catch(...){
      subset.pop_back();
      continue;
    }
    return maxCost;
  }
}

//sequential forward selection Multivariate (Combination of N best features)
template<class T> double FeatureSelector::forward(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures=0, short verbose){
  //Select feature with the best value (get maximal cost for 1 feature)
  double maxCost=0;
  int maxLevels=v[0][0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  while(subset.size()<maxFeatures){
    maxCost=addFeature(v,theCostFactory,subset,verbose);
    if(verbose){
      for(std::list<int>::const_iterator lit=subset.begin();lit!=subset.end();++lit)
        std::cout << *lit << " ";
      std::cout << std::endl;
      // std::cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << maxCost << std::endl;
    }
  }//while
  return maxCost;
}

//sequential backward selection
template<class T> double FeatureSelector::backward(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int minFeatures, short verbose){
  //Select features with least effect on cost when removed (obtain minFeatures eventually)
  double maxCost=0;
  int removedFeature;
  if(subset.empty()){
    for(int iFeature=0;iFeature<v[0][0].size();++iFeature)
      subset.push_back(iFeature);
  }
  if(subset.size()==minFeatures)
    maxCost=theCostFactory.getCost(v);
  while(subset.size()>minFeatures){
    maxCost=removeFeature(v,theCostFactory,subset,removedFeature,verbose);
    if(verbose){
      for(std::list<int>::const_iterator lit=subset.begin();lit!=subset.end();++lit)
        std::cout << *lit << " ";
      std::cout << std::endl;
      // std::cout << "removed " << removedFeature << ", " << subset.size() << "/" << minFeatures << " features remain with cost: " << maxCost << std::endl;
    }
  }//while
  return maxCost;
}

//floating search
template<class T> double FeatureSelector::floating(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures, double epsilon, short verbose){
  std::vector<T> cost;
  int maxLevels=v[0][0].size();
  if(maxFeatures<1)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
  while(cost.size()<subset.size())
    cost.push_back(1);//init original features as cost 1
  cost.push_back(addFeature(v,theCostFactory,subset,verbose));
  ++k;
  if(verbose>1)
    std::cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << cost.back() << std::endl;
  else if(verbose){
    for(std::list<int>::const_iterator lit=subset.begin();lit!=subset.end();++lit)
      std::cout << *lit << " ";
    std::cout << std::endl;
  }
  while(k<maxFeatures){
    cost.push_back(addFeature(v,theCostFactory,subset,verbose));
    ++k;
    if(verbose>1)
    std::cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << cost.back() << std::endl;
    else if(verbose){
      for(std::list<int>::const_iterator lit=subset.begin();lit!=subset.end();++lit)
        std::cout << *lit << " ";
      std::cout << " (cost: " << cost.back() << ")" << std::endl;
    }

    while(k>1){
      int x_r;
      double cost_r=removeFeature(v,theCostFactory,subset,x_r,verbose);
      if(cost_r>cost[k-1]+epsilon){
	--k;
	cost[k]=cost_r;
	cost.pop_back();
	if(verbose>1)
	  std::cout << "removed " << x_r << ", " << subset.size() << "/" << maxFeatures << " features remain with cost: " << cost_r << std::endl;
        else if(verbose){
          for(std::list<int>::const_iterator lit=subset.begin();lit!=subset.end();++lit)
            std::cout << *lit << " ";
          std::cout << " (cost: " << cost.back() << ")" << std::endl;
        }
	continue;
      }
      else if(cost_r>=0){
	subset.push_back(x_r);
	break;
      }
      else if(verbose)
        std::cout << "could not remove any feature" << std::endl;
      cost.pop_back();
    }
  }
  return cost.back();
}

//brute force search (search for all possible combinations and select best)
template<class T> double FeatureSelector::bruteForce(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int maxFeatures, short verbose){
  int maxLevels=v[0][0].size();
  if(maxFeatures<1)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
//   gslmm::combination c(v1[0].size(),maxFeatures,false);
  gsl_combination *c;
  c=gsl_combination_calloc(v[0][0].size(),maxFeatures);
  
  std::list<int> tmpset;//temporary set of selected features (levels)
  std::vector< Vector2d<T> > tmpv(v.size());
  std::list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
  double cost;
  if(subset.size()>=maxLevels)
    return maxCost;
  gsl_combination_next(c);
  do{
    for(int ifeature=0;ifeature<maxFeatures;++ifeature)
      tmpset.push_back(c->data[ifeature]);
    for(int iclass=0;iclass<v.size();++iclass)
      v[iclass].selectCols(tmpset,tmpv[iclass]);
    try{
      cost=theCostFactory.getCost(tmpv);
    }
    catch(...){ //singular matrix encountered
      catchset=tmpset;//this tmpset resulted in failure of getCost
      if(verbose){
        std::cout << "Could not get cost from set: " << std::endl;
	gsl_combination_fprintf(stdout,c," %u");
	printf("\n");
      }      
      tmpset.clear();
      continue;
    }
    if(maxCost<cost){ //set with better cost is found
      maxCost=cost;
      subset=tmpset;
    }
    tmpset.clear();
  }while(gsl_combination_next(c)==GSL_SUCCESS);
  gsl_combination_free(c);
//   }while(c.next());  
  if(maxCost<0) //no level added to better maxCost than current subset (catchset)
    subset=catchset;
  //test: assert list contains no duplicate elements
  for(std::list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    std::list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::addFeature(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, short verbose){
  //select feature with the best value (get maximal cost for 1 feature)
  std::list<int> tmpset=subset;//temporary set of selected features (levels)
  std::vector< Vector2d<T> > tmp(v.size());
  std::list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
  double cost;
  int maxLevels=v[0][0].size();
  if(subset.size()>=maxLevels)
    return maxCost;
  for(int ilevel=0;ilevel<maxLevels;++ilevel){
    if(find(tmpset.begin(),tmpset.end(),ilevel)!=tmpset.end())
      continue;
    tmpset.push_back(ilevel);
    for(int iclass=0;iclass<v.size();++iclass){
//       tmp[iclass].resize(v[iclass].size());
      v[iclass].selectCols(tmpset,tmp[iclass]);
    }
    try{
      cost=theCostFactory.getCost(tmp);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
        std::cout << "Could not add feature " << tmpset.back() << std::endl;
      tmpset.pop_back();
      continue;
    }
    if(maxCost<cost){ //level with better cost is found
      maxCost=cost;
      subset=tmpset;
    }
    tmpset.pop_back();
  }
  if(maxCost<0) //no level added to better maxCost than current subset (catchset)
    subset=catchset;
  //test: assert list contains no duplicate elements
  for(std::list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    std::list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::removeFeature(std::vector< Vector2d<T> >& v, CostFactory& theCostFactory, std::list<int>& subset, int& r, short verbose){
  //find the feature that has the least effect on the cost when it is removed from subset
  std::list<int> tmpset=subset;//temporary set of selected features (levels)
  std::vector< Vector2d<T> > tmp(v.size());
  int nFeatures=subset.size();
  std::list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
  int last;
  double cost;
  int maxLevels=v[0][0].size();
  if(subset.size()>maxLevels||subset.empty()){
    return maxCost;
  }
  std::list<int>::const_iterator it;
  for(int i=0;i<nFeatures;++i){
    last=tmpset.back();
    tmpset.pop_back();
    for(int iclass=0;iclass<v.size();++iclass){
//       tmp[iclass].resize(v[iclass].size());
      v[iclass].selectCols(tmpset,tmp[iclass]);
    }
    try{
      cost=theCostFactory.getCost(tmp);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
        std::cout << "Could not remove feature " << last << std::endl;
      tmpset.push_front(last);
      continue;
    }
    if(maxCost<cost){ //level with better cost is found
      maxCost=cost;
      subset=tmpset;
      r=last;
    }
    tmpset.push_front(last);
  }
  if(maxCost<0){//all levels removed were catched
    subset=catchset;
  }
  //test: assert list contains no duplicate elements
  for(std::list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    std::list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

#endif /* _FEATURESELECTOR_H_ */

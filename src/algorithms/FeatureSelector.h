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
#include <map>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include "PosValue.h"
#include "Vector2d.h"
#include "gsl/gsl_combination.h"

using namespace std;

class FeatureSelector
{
 public:
  FeatureSelector(){};
  ~FeatureSelector(){};
  template<class T> double forwardUnivariate(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures, double prior1=0.5, double prior2=0.5, bool verbose=false);
  template<class T> double forwardUnivariate(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int maxFeatures, bool verbose=false);
  template<class T> double forwardUnivariate(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int maxFeatures, bool verbose=false);
    template<class T> double forward(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures, double prior1=0.5, double prior2=0.5, bool verbose=false);
    template<class T> double forward(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int maxFeatures, bool verbose=false);
  template<class T> double forward(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int maxFeatures, bool verbose=false);
    template<class T> double backward(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int minFeatures, bool verbose=false);
    template<class T> double backward(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int minFeatures, double prior1=0.5, double prior2=0.5,bool verbose=false);
    template<class T> double backward(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int minFeatures, bool verbose=false);
    template<class T> double floating(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int maxFeatures=0, bool verbose=false);
    template<class T> double floating(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures=0, double prior1=0.5, double prior2=0.5,bool verbose=false);
    template<class T> double floating(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int maxFeatures, bool verbose=false);
  template<class T> double bruteForce(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures=0, double prior1=0.5, double prior2=0.5, bool verbose=false);
  
  private:
    template<class T> double addFeature(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, bool verbose=false);
    template<class T> double addFeature(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, double prior1=0.5, double prior2=0.5, bool verbose=false);
    template<class T> double addFeature(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, bool verbose=false);
  template<class T> double removeFeature(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int& r, bool verbose=false);  
    template<class T> double removeFeature(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset,int& r, double prior1=0.5, double prior2=0.5, bool verbose=false);
    template<class T> double removeFeature(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset,int& r, bool verbose=false);
};

//sequential forward selection Univariate (N single best features)
template<class T> double FeatureSelector::forwardUnivariate(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures, double prior1, double prior2, bool verbose){
  int maxLevels=v1[0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
  vector<PosValue> cost(maxLevels);
  list<int> tmpset=subset;//temporary set of selected features (levels)
  Vector2d<T> tmp1(v1.size());
  Vector2d<T> tmp2(v2.size());
  for(int ilevel=0;ilevel<maxLevels;++ilevel){
    if(find(tmpset.begin(),tmpset.end(),ilevel)==tmpset.end()){
      tmpset.push_back(ilevel);
      v1.selectCols(tmpset,tmp1);
      v2.selectCols(tmpset,tmp2);
      try{
	PosValue pv;
	pv.position=ilevel;
	pv.value=getCost(tmp1,tmp2,prior1,prior2);
	cost[ilevel]=pv;
      }
      catch(...){
	PosValue pv;
	pv.position=ilevel;
	pv.value=-1;
	cost[ilevel]=pv;
      }
      tmpset.pop_back();
    }
  }
  sort(cost.begin(),cost.end(),Compare_PosValue());//decreasing order
  int ilevel=0;
  while((subset.size()<maxFeatures)&&(ilevel<maxLevels)){
    if(cost[ilevel].value>0)
      subset.push_back(cost[ilevel].position);
    if(verbose)
      cout << "feature " << subset.back() << " has cost: " << cost[ilevel].value << endl;
    ++ilevel;
  }
  double maxCost=-1;
  while(subset.size()){
    v1.selectCols(subset,tmp1);
    v2.selectCols(subset,tmp2);
    try{
      maxCost=getCost(tmp1,tmp2,prior1,prior2);
    }
    catch(...){
      subset.pop_back();
      continue;
    }
    return maxCost;
  }
}

template<class T> double FeatureSelector::forwardUnivariate(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int maxFeatures=0, bool verbose){
  int maxLevels=v[0][0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
  vector<PosValue> cost(maxLevels);
  list<int> tmpset=subset;//temporary set of selected features (levels)
  vector< Vector2d<T> > tmp(v.size());
  for(int ilevel=0;ilevel<maxLevels;++ilevel){
    if(find(tmpset.begin(),tmpset.end(),ilevel)==tmpset.end()){
      tmpset.push_back(ilevel);
      for(int iclass=0;iclass<v.size();++iclass){
// 	tmp[iclass].resize(v[iclass].size());
	v[iclass].selectCols(tmpset,tmp[iclass]);
      }
      try{
	PosValue pv;
	pv.position=ilevel;
	pv.value=getCost(tmp);
	cost[ilevel]=pv;
      }
      catch(...){
	PosValue pv;
	pv.position=ilevel;
	pv.value=-1;
	cost[ilevel]=pv;
      }
      tmpset.pop_back();
    }
  }
  sort(cost.begin(),cost.end(),Compare_PosValue());//decreasing order
  int ilevel=0;
  while((subset.size()<maxFeatures)&&(ilevel<maxLevels)){
    if(cost[ilevel].value>0)
      subset.push_back(cost[ilevel].position);
    if(verbose)
      cout << "feature " << subset.back() << " has cost: " << cost[ilevel].value << endl;
    ++ilevel;
  }
  double maxCost=-1;
  while(subset.size()){
    for(int iclass=0;iclass<v.size();++iclass){
//       tmp[iclass].resize(v[iclass].size());
      v[iclass].selectCols(subset,tmp[iclass]);
    }
    try{
      maxCost=getCost(tmp);
    }
    catch(...){
      subset.pop_back();
      continue;
    }
    return maxCost;
  }
}

//sequential forward selection Multivariate (Combination of N best features)
template<class T> double FeatureSelector::forward(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures, double prior1, double prior2, bool verbose){
  //Select feature with the best value (get maximal cost for 1 feature)
  double maxCost=0;
  int maxLevels=v1[0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  while(subset.size()<maxFeatures){
    maxCost=addFeature(v1,v2,*getCost,subset,verbose);
    if(verbose)
      cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << maxCost << endl;
  }//while
  return maxCost;
}


//sequential forward selection Multivariate (Combination of N best features)
template<class T> double FeatureSelector::forward(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int maxFeatures=0, bool verbose){
  //Select feature with the best value (get maximal cost for 1 feature)
  double maxCost=0;
  int maxLevels=v[0][0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  while(subset.size()<maxFeatures){
    maxCost=addFeature(v,*getCost,subset,verbose);
    if(verbose)
      cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << maxCost << endl;
  }//while
  return maxCost;
}

template<class T> double FeatureSelector::forward(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int maxFeatures, bool verbose)
{
  //Select feature with the best value (get maximal cost for 1 feature)
  double maxCost=0;
  int maxLevels=(v.begin()->second)[0].size();
  if(!maxFeatures)
    maxFeatures=maxLevels;
  while(subset.size()<maxFeatures){
    maxCost=addFeature(v,y,*getCost,subset,verbose);
    if(verbose)
      cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << maxCost << endl;
  }//while
  return maxCost;
}

//sequential backward selection
template<class T> double FeatureSelector::backward(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int minFeatures, double prior1, double prior2, bool verbose){
  //Select features with least effect on cost when removed (obtain minFeatures eventually)
  double maxCost=0;
  int removedFeature;
  if(subset.empty()){
    for(int iFeature=0;iFeature<v1[0].size();++iFeature)
      subset.push_back(iFeature);
  }//if
  if(subset.size()==minFeatures)
    maxCost=getCost(v1,v2,prior1,prior2);
  while(subset.size()>minFeatures){
    maxCost=removeFeature(v1,v2,*getCost,subset,removedFeature,verbose);
    if(verbose)
      cout << "removed " << removedFeature << ", " << subset.size() << "/" << minFeatures << " features remain with cost: " << maxCost << endl;
  }//while
  return maxCost;
}

//sequential backward selection
template<class T> double FeatureSelector::backward(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int minFeatures, bool verbose){
  //Select features with least effect on cost when removed (obtain minFeatures eventually)
  double maxCost=0;
  int removedFeature;
  if(subset.empty()){
    for(int iFeature=0;iFeature<v[0][0].size();++iFeature)
      subset.push_back(iFeature);
  }//if
  if(subset.size()==minFeatures)
    maxCost=getCost(v);
  while(subset.size()>minFeatures){
    maxCost=removeFeature(v,*getCost,subset,removedFeature,verbose);
    if(verbose)
      cout << "removed " << removedFeature << ", " << subset.size() << "/" << minFeatures << " features remain with cost: " << maxCost << endl;
  }//while
  return maxCost;
}

template<class T> double FeatureSelector::backward(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int minFeatures, bool verbose){
  //Select features with least effect on cost when removed (obtain minFeatures eventually)
  double maxCost=0;
  int removedFeature;
  if(subset.empty()){
    for(int iFeature=0;iFeature<(v.begin()->second)[0].size();++iFeature)
      subset.push_back(iFeature);
  }//if
  if(subset.size()==minFeatures)
    maxCost=getCost(v,y);
  while(subset.size()>minFeatures){
    maxCost=removeFeature(v,y,*getCost,subset,removedFeature,verbose);
    if(verbose)
      cout << "removed " << removedFeature << ", " << subset.size() << "/" << minFeatures << " features remain with cost: " << maxCost << endl;
  }//while
  return maxCost;
}

//floating search
template<class T> double FeatureSelector::floating(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures, double prior1, double prior2, bool verbose){
  vector<T> cost;
  int maxLevels=v1[0].size();
  if(maxFeatures<1)
    maxFeatures=maxLevels;
//   int k=0;
  int initialK=subset.size();
  int k=initialK;
  int stopK=initialK+1;
  if(initialK>=maxFeatures)
    return -1;
  cost.push_back(-1);//init 0 features as cost -1
  cost.push_back(addFeature(v1,v2,*getCost,subset,verbose));
  ++k;
  assert(subset.size()==k);
  if(verbose)
    cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << cost.back() << endl;
  while(k<maxFeatures){
    cost.push_back(addFeature(v1,v2,*getCost,subset,verbose));
    ++k;
    assert(subset.size()==k);
    if(verbose)
    cout << "added " << subset.back() << ", " << subset.size() << "/" << maxFeatures << " features selected with cost: " << cost.back() << endl;
//     while(k>1){
    while(k>stopK){
      int x_r;
      double cost_r=removeFeature(v1,v2,*getCost,subset,x_r,verbose);
      if(cost_r>cost[k-1]){
	--k;
	assert(subset.size()==k);
	cost[k]=cost_r;
	cost.pop_back();
	if(verbose)
	  cout << "removed " << x_r << ", " << subset.size() << "/" << maxFeatures << " features remain with cost: " << cost_r << endl;
	continue;
      }
      else if(cost_r>=0){
	subset.push_back(x_r);
	break;
      }
      else if(verbose)
	cout << "could not remove any feature" << endl;
      cost.pop_back();
    }
  }
  return cost.back();
}

//floating search
template<class T> double FeatureSelector::floating(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int maxFeatures, bool verbose){
  vector<T> cost;
  int maxLevels=v[0][0].size();
  if(maxFeatures<1)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
  cost.push_back(-1);//init 0 features as cost -1
  cost.push_back(addFeature(v,*getCost,subset,verbose));
  ++k;
  if(verbose)
    cout << "added " << subset.back() << ", " << cost.size()-1 << "/" << maxFeatures << " features selected with cost: " << cost.back() << endl;
  while(k<maxFeatures){
    cost.push_back(addFeature(v,*getCost,subset,verbose));
    ++k;
    if(verbose)
    cout << "added " << subset.back() << ", " << cost.size()-1 << "/" << maxFeatures << " features selected with cost: " << cost.back() << endl;
    while(k>1){
      int x_r;
      double cost_r=removeFeature(v,*getCost,subset,x_r,verbose);
      if(cost_r>cost[k-1]){
	--k;
	cost[k]=cost_r;
	cost.pop_back();
	if(verbose)
	  cout << "removed " << x_r << ", " << cost.size()-1 << "/" << maxFeatures << " features remain with cost: " << cost_r << endl;
	continue;
      }
      else if(cost_r>=0){
	subset.push_back(x_r);
	break;
      }
      else if(verbose)
	cout << "could not remove any feature" << endl;
      cost.pop_back();
    }
  }
  return cost.back();
}

//floating search
template<class T> double FeatureSelector::floating(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int maxFeatures, bool verbose){
  vector<T> cost;
  int maxLevels=(v.begin()->second)[0].size();
  if(maxFeatures<1)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
  cost.push_back(-1);//init 0 features as cost -1
  cost.push_back(addFeature(v,y,*getCost,subset,verbose));
  ++k;
  if(verbose)
    cout << "added " << subset.back() << ", " << cost.size()-1 << "/" << maxFeatures << " features selected with cost: " << cost.back() << endl;
  while(k<maxFeatures){
    cost.push_back(addFeature(v,y,*getCost,subset,verbose));
    ++k;
    if(verbose)
    cout << "added " << subset.back() << ", " << cost.size()-1 << "/" << maxFeatures << " features selected with cost: " << cost.back() << endl;
    while(k>1){
      int x_r;
      double cost_r=removeFeature(v,y,*getCost,subset,x_r,verbose);
      if(cost_r>cost[k-1]){
	--k;
	cost[k]=cost_r;
	cost.pop_back();
	if(verbose)
	  cout << "removed " << x_r << ", " << cost.size()-1 << "/" << maxFeatures << " features remain with cost: " << cost_r << endl;
	continue;
      }
      else if(cost_r>=0){
	subset.push_back(x_r);
	break;
      }
      else if(verbose)
	cout << "could not remove any feature" << endl;
      cost.pop_back();
    }
  }
  return cost.back();
}

//brute force search (search for all possible combinations and select best)
template<class T> double FeatureSelector::bruteForce(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int maxFeatures, double prior1, double prior2, bool verbose){
  int maxLevels=v1[0].size();
  if(maxFeatures<1)
    maxFeatures=maxLevels;
  int k=subset.size();
  if(k>=maxFeatures)
    return -1;
//   gslmm::combination c(v1[0].size(),maxFeatures,false);
  gsl_combination *c;
  c=gsl_combination_calloc(v1[0].size(),maxFeatures);
  
  list<int> tmpset;//temporary set of selected features (levels)
  Vector2d<T> tmp1(v1.size());
  Vector2d<T> tmp2(v2.size());
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
  double cost;
  if(subset.size()>=maxLevels)
    return maxCost;
//   c.first();
  gsl_combination_next(c);
  do{
    for(int ifeature=0;ifeature<maxFeatures;++ifeature)
      tmpset.push_back(c->data[ifeature]);
    v1.selectCols(tmpset,tmp1);
    v2.selectCols(tmpset,tmp2);
    try{
      cost=getCost(tmp1,tmp2,prior1,prior2);
    }
    catch(...){ //singular matrix encountered
      catchset=tmpset;//this tmpset resulted in failure of getCost
      if(verbose){
	cout << "Could not get cost from set: " << endl;
	gsl_combination_fprintf(stdout,c," %u");
	printf("\n");
// 	c.print(stdout);
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::addFeature(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, double prior1, double prior2, bool verbose){
  //select feature with the best value (get maximal cost for 1 feature)
  list<int> tmpset=subset;//temporary set of selected features (levels)
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
//   if(subset.size()){
//     try{
//       v1.selectCols(subset,tmp1);
//       v2.selectCols(subset,tmp2);
//       maxCost=getCost(tmp1,tmp2);
//     }
//     catch(double determinant){
//       return -1;
//     }
//   }
  double cost;
  int maxLevels=v1[0].size();
  if(subset.size()>=maxLevels)
    return maxCost;
  for(int ilevel=0;ilevel<maxLevels;++ilevel){
    Vector2d<T> tmp1(v1.size());
    Vector2d<T> tmp2(v2.size());
    if(find(tmpset.begin(),tmpset.end(),ilevel)!=tmpset.end())
      continue;
    tmpset.push_back(ilevel);
    v1.selectCols(tmpset,tmp1);
    v2.selectCols(tmpset,tmp2);
    try{
      cost=getCost(tmp1,tmp2,prior1,prior2);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
	cout << "Could not add feature " << tmpset.back() << endl;
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::addFeature(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, bool verbose){
  //select feature with the best value (get maximal cost for 1 feature)
  list<int> tmpset=subset;//temporary set of selected features (levels)
  vector< Vector2d<T> > tmp(v.size());
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
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
      cost=getCost(tmp);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
	cout << "Could not add feature " << tmpset.back() << endl;
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::addFeature(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, bool verbose){
  //select feature with the best value (get maximal cost for 1 feature)
  list<int> tmpset=subset;//temporary set of selected features (levels)
  map<string, Vector2d<double> > tmp;
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  map<string, Vector2d<double> >::const_iterator mit;
  double maxCost=-1;
  double cost;
  int maxLevels=(v.begin()->second)[0].size();
  if(subset.size()>=maxLevels)
    return maxCost;
  for(int ilevel=0;ilevel<maxLevels;++ilevel){
    if(find(tmpset.begin(),tmpset.end(),ilevel)!=tmpset.end())
      continue;
    tmpset.push_back(ilevel);
    for(mit=v.begin();mit!=v.end();++mit){
      string theName=mit->first;
      v[theName].selectCols(tmpset,tmp[theName]);
    }
    try{
      cost=getCost(tmp,y);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
	cout << "Could not add feature " << tmpset.back() << endl;
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::removeFeature(Vector2d<T>& v1, Vector2d<T>& v2, double (*getCost)(const Vector2d<T>&, const Vector2d<T>&, double prior1, double prior2), list<int>& subset, int& r, double prior1, double prior2, bool verbose){
  //find the feature that has the least effect on the cost when it is removed from subset
  list<int> tmpset=subset;//temporary set of selected features (levels)
  Vector2d<T> tmp1(v1.size());
  Vector2d<T> tmp2(v2.size());
  int nFeatures=subset.size();
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
//   if(subset.size()){
//     try{
//       v1.selectCols(subset,tmp1);
//       v2.selectCols(subset,tmp2);
//       maxCost=getCost(tmp1,tmp2);
//     }
//     catch(double determinant){
//       return -1;
//     }
//   }
  int last;
  double cost;
  int maxLevels=v1[0].size();
  if(subset.size()>maxLevels||subset.empty()){
    return maxCost;
  }
  list<int>::const_iterator it;
  for(int i=0;i<nFeatures;++i){
    last=tmpset.back();
    tmpset.pop_back();
    v1.selectCols(tmpset,tmp1);
    v2.selectCols(tmpset,tmp2);
    try{
      cost=getCost(tmp1,tmp2,prior1,prior2);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
	cout << "Could not remove feature " << last << endl;
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::removeFeature(vector< Vector2d<T> >& v, double (*getCost)(const vector< Vector2d<T> >&), list<int>& subset, int& r, bool verbose){
  //find the feature that has the least effect on the cost when it is removed from subset
  list<int> tmpset=subset;//temporary set of selected features (levels)
  vector< Vector2d<T> > tmp(v.size());
  int nFeatures=subset.size();
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
  int last;
  double cost;
  int maxLevels=v[0][0].size();
  if(subset.size()>maxLevels||subset.empty()){
    return maxCost;
  }
  list<int>::const_iterator it;
  for(int i=0;i<nFeatures;++i){
    last=tmpset.back();
    tmpset.pop_back();
    for(int iclass=0;iclass<v.size();++iclass){
//       tmp[iclass].resize(v[iclass].size());
      v[iclass].selectCols(tmpset,tmp[iclass]);
    }
    try{
      cost=getCost(tmp);
    }
    catch(...){
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
	cout << "Could not remove feature " << last << endl;
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}

template<class T> double FeatureSelector::removeFeature(map<string, Vector2d<T> > &v, map<string, vector<T> > &y, double (*getCost)(map<string, Vector2d<T> >&, map<string, vector<T> >&), list<int>& subset, int& r, bool verbose){
  //find the feature that has the least effect on the cost when it is removed from subset
  list<int> tmpset=subset;//temporary set of selected features (levels)
  map<string, Vector2d<double> > tmp;
  map<string, Vector2d<double> >::const_iterator mit;
  int nFeatures=subset.size();
  list<int> catchset;//restore set in case of catch all the way to last level (no better cost)
  //initialize maxCost with actual cost for current subset (-1 if empty subset) 
  double maxCost=-1;
  int last;
  double cost;
  int maxLevels=(v.begin()->second)[0].size();
  if(subset.size()>maxLevels||subset.empty()){
    return maxCost;
  }
  list<int>::const_iterator it;
  for(int i=0;i<nFeatures;++i){
    last=tmpset.back();
    tmpset.pop_back();
    for(mit=v.begin();mit!=v.end();++mit){
      string theName=mit->first;
      v[theName].selectCols(tmpset,tmp[theName]);
    }
    try{
      cost=getCost(tmp,y);
    }
    catch(...){ //singular matrix encountered
      catchset=tmpset;//this tmpset resulted in singular matrix
      if(verbose)
	cout << "Could not remove feature " << last << endl;
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
  for(list<int>::iterator lit=subset.begin();lit!=subset.end();++lit){
    list<int>::iterator lit2=lit;//start searching from next element
    assert(find(++lit2,subset.end(),*lit)==subset.end());
  }
  return maxCost;
}
#endif /* _FEATURESELECTOR_H_ */

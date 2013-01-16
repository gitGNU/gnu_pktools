/**********************************************************************
OptFactory.h: factory class for nlopt::opt (selecting algorithm via string)
Copyright (C) 2008-2013 Pieter Kempeneers

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
#ifndef _OPTFACTORY_H_
#define _OPTFACTORY_H_

#include <map>

class OptFactory
{
private:
  static void initMap(std::map<std::string, nlopt::algorithm>& m_algMap){
    //initialize selMap
    m_algMap["GN_DIRECT"]=nlopt::GN_DIRECT;
    m_algMap["GN_DIRECT_L"]=nlopt::GN_DIRECT_L;
    m_algMap["GN_DIRECT_L_RAND"]=nlopt::GN_DIRECT_L_RAND;
    m_algMap["GN_DIRECT_NOSCAL"]=nlopt::GN_DIRECT_NOSCAL;
    m_algMap["GN_DIRECT_L_NOSCAL"]=nlopt::GN_DIRECT_L_NOSCAL;
    m_algMap["GN_DIRECT_L_RAND_NOSCAL"]=nlopt::GN_DIRECT_L_RAND_NOSCAL;
    m_algMap["GN_ORIG_DIRECT"]=nlopt::GN_ORIG_DIRECT;
    m_algMap["GN_ORIG_DIRECT_L"]=nlopt::GN_ORIG_DIRECT_L;
    m_algMap["GD_STOGO"]=nlopt::GD_STOGO;
    m_algMap["GD_STOGO_RAND"]=nlopt::GD_STOGO_RAND;
    m_algMap["LD_LBFGS_NOCEDAL"]=nlopt::LD_LBFGS_NOCEDAL;
    m_algMap["LD_LBFGS"]=nlopt::LD_LBFGS;
    m_algMap["LN_PRAXIS"]=nlopt::LN_PRAXIS;
    m_algMap["LD_VAR1"]=nlopt::LD_VAR1;
    m_algMap["LD_VAR2"]=nlopt::LD_VAR2;
    m_algMap["LD_TNEWTON"]=nlopt::LD_TNEWTON;
    m_algMap["LD_TNEWTON_RESTART"]=nlopt::LD_TNEWTON_RESTART;
    m_algMap["LD_TNEWTON_PRECOND"]=nlopt::LD_TNEWTON_PRECOND;
    m_algMap["LD_TNEWTON_PRECOND_RESTART"]=nlopt::LD_TNEWTON_PRECOND_RESTART;
    m_algMap["GN_CRS2_LM"]=nlopt::GN_CRS2_LM;
    m_algMap["GN_MLSL"]=nlopt::GN_MLSL;
    m_algMap["GD_MLSL"]=nlopt::GD_MLSL;
    m_algMap["GN_MLSL_LDS"]=nlopt::GN_MLSL_LDS;
    m_algMap["GD_MLSL_LDS"]=nlopt::GD_MLSL_LDS;
    m_algMap["LD_MMA"]=nlopt::LD_MMA;
    m_algMap["LN_COBYLA"]=nlopt::LN_COBYLA;
    m_algMap["LN_NEWUOA"]=nlopt::LN_NEWUOA;
    m_algMap["LN_NEWUOA_BOUND"]=nlopt::LN_NEWUOA_BOUND;
    m_algMap["LN_NELDERMEAD"]=nlopt::LN_NELDERMEAD;
    m_algMap["LN_SBPLX"]=nlopt::LN_SBPLX;
    m_algMap["LN_AUGLAG"]=nlopt::LN_AUGLAG;
    m_algMap["LD_AUGLAG"]=nlopt::LD_AUGLAG;
    m_algMap["LN_AUGLAG_EQ"]=nlopt::LN_AUGLAG_EQ;
    m_algMap["LD_AUGLAG_EQ"]=nlopt::LD_AUGLAG_EQ;
    m_algMap["LN_BOBYQA"]=nlopt::LN_BOBYQA;
    m_algMap["GN_ISRES"]=nlopt::GN_ISRES;
    m_algMap["AUGLAG"]=nlopt::AUGLAG;
    m_algMap["AUGLAG_EQ"]=nlopt::AUGLAG_EQ;
    m_algMap["G_MLSL"]=nlopt::G_MLSL;
    m_algMap["G_MLSL_LDS"]=nlopt::G_MLSL_LDS;
    m_algMap["LD_SLSQP "]=nlopt::LD_SLSQP;
  }
public:
  OptFactory(){
  };
  ~OptFactory(){};
  static nlopt::opt getOptimizer(const std::string& algorithmString, unsigned int npar){
    std::map<std::string, nlopt::algorithm> m_algMap;
    initMap(m_algMap);
    switch(m_algMap[algorithmString]){
    case(nlopt::GN_DIRECT):{
      nlopt::opt theOptimizer(nlopt::GN_DIRECT,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_DIRECT_L):{
      nlopt::opt theOptimizer(nlopt::GN_DIRECT_L,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_DIRECT_L_RAND):{
      nlopt::opt theOptimizer(nlopt::GN_DIRECT_L_RAND,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_DIRECT_NOSCAL):{
      nlopt::opt theOptimizer(nlopt::GN_DIRECT_NOSCAL,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_DIRECT_L_NOSCAL):{
      nlopt::opt theOptimizer(nlopt::GN_DIRECT_L_NOSCAL,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_DIRECT_L_RAND_NOSCAL):{
      nlopt::opt theOptimizer(nlopt::GN_DIRECT_L_RAND_NOSCAL,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_ORIG_DIRECT):{
      nlopt::opt theOptimizer(nlopt::GN_ORIG_DIRECT,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_ORIG_DIRECT_L):{
      nlopt::opt theOptimizer(nlopt::GN_ORIG_DIRECT_L,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_PRAXIS):{
      nlopt::opt theOptimizer(nlopt::LN_PRAXIS,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_CRS2_LM):{
      nlopt::opt theOptimizer(nlopt::GN_CRS2_LM,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_MLSL):{
      nlopt::opt theOptimizer(nlopt::GN_MLSL,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_MLSL_LDS):{
      nlopt::opt theOptimizer(nlopt::GN_MLSL_LDS,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_COBYLA):{
      nlopt::opt theOptimizer(nlopt::LN_COBYLA,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_NEWUOA):{
      nlopt::opt theOptimizer(nlopt::LN_NEWUOA,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_NEWUOA_BOUND):{
      nlopt::opt theOptimizer(nlopt::LN_NEWUOA_BOUND,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_NELDERMEAD):{
      nlopt::opt theOptimizer(nlopt::LN_NELDERMEAD,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_SBPLX):{
      nlopt::opt theOptimizer(nlopt::LN_SBPLX,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_AUGLAG):{
      nlopt::opt theOptimizer(nlopt::LN_AUGLAG,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_AUGLAG_EQ):{
      nlopt::opt theOptimizer(nlopt::LN_AUGLAG_EQ,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::LN_BOBYQA):{
      nlopt::opt theOptimizer(nlopt::LN_BOBYQA,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::GN_ISRES):{
      nlopt::opt theOptimizer(nlopt::GN_ISRES,npar);
      return(theOptimizer);
      break;
    }
    case(nlopt::G_MLSL_LDS):
    case(nlopt::AUGLAG):
    case(nlopt::AUGLAG_EQ):
    case(nlopt::G_MLSL):
    case(nlopt::GD_MLSL):
    case(nlopt::GD_MLSL_LDS):
    case(nlopt::GD_STOGO):
    case(nlopt::GD_STOGO_RAND):
    case(nlopt::LD_LBFGS_NOCEDAL):
    case(nlopt::LD_LBFGS):
    case(nlopt::LD_VAR1):
    case(nlopt::LD_VAR2):
    case(nlopt::LD_TNEWTON):
    case(nlopt::LD_TNEWTON_RESTART):
    case(nlopt::LD_TNEWTON_PRECOND):
    case(nlopt::LD_TNEWTON_PRECOND_RESTART):
    case(nlopt::LD_MMA):
    case(nlopt::LD_AUGLAG):
    case(nlopt::LD_AUGLAG_EQ):
    case(nlopt::LD_SLSQP):
    default:{
      std::string errorString="Error: derivative optimization algorithm ";
      errorString+= algorithmString;
      errorString+= " not supported, select derivative-free algorithm ([GL]N_*])";
      throw(errorString);
      break;
    }
    }
  };
};
#endif /* _OPTFACTORY_H_ */

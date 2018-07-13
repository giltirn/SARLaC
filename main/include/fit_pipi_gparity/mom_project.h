#ifndef _PIPI_MOM_PROJECT_H_
#define _PIPI_MOM_PROJECT_H_

#include<set>
#include<regex>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#include "enums.h"

//When reading and projecting the data onto a desired rotation state, we loop over all choices of p1_src and p1_snk, and use a "Selector" to pick out which combinations should be included
//and their coefficient in the projection
struct PiPiCorrelatorSelector{
  virtual bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const = 0;
};


//A basic but quite general selector can be created by combining a "Project" instance for source and sink which control the coefficients of the pipi operator in the projection,
//and an overall filter on (pi1_src, pi1_snk)
struct PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const = 0;
};
struct PiPiMomAllow{
  virtual bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const = 0;
};

//We use factory functions to generate the instance
PiPiProject* getProjector(const PiPiProjector p, const threeMomentum &ptot);
PiPiMomAllow* getMomPairFilter(const PiPiMomAllowed p, const threeMomentum &ptot);

//The selector
struct PiPiCorrelatorBasicSelector: public PiPiCorrelatorSelector{
  PiPiProject* proj_src;
  PiPiProject* proj_snk;
  PiPiMomAllow* allow;
  bool own_ptrs;

  //p_tot is the total momentum of the pipi at the source (p1_src + p2_src)
  PiPiCorrelatorBasicSelector(PiPiProjector proj_src_, PiPiProjector proj_snk_, PiPiMomAllowed allow_, const threeMomentum &p_tot): own_ptrs(true){
    proj_src = getProjector(proj_src_, p_tot);
    proj_snk = getProjector(proj_snk_, -p_tot);
    allow = getMomPairFilter(allow_, p_tot);
  }
  PiPiCorrelatorBasicSelector(PiPiProject* proj_src, PiPiProject* proj_snk, PiPiMomAllow* allow): proj_src(proj_src), proj_snk(proj_snk), allow(allow), own_ptrs(false){}

  ~PiPiCorrelatorBasicSelector(){
    if(own_ptrs){ delete proj_src; delete proj_snk; delete allow; }
  }

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    std::complex<double> csrc, csnk;
    double m;
    if( (*proj_src)(csrc,p1src) && (*proj_snk)(csnk,p1snk) && (*allow)(m, p1src, p1snk) ){
      coeff = std::real(csrc * csnk)*m;
      return true;
    }else return false;
  }
};

//------------------------------------------ Projectors -------------------------------------------------------------

struct PiPiProjectA1: public PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const{ //assumes p_tot = 0 (unchecked)
    if(abs(mom[0])==1 && abs(mom[1])==1 && abs(mom[2])==1){
      coeff = 1./8;
      return true;
    }else{
      return false;
    }
  }
};
struct PiPiProjectAvg4: public PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const{
    coeff = 1./4;
    return true;
  }
};
struct PiPiProjectAvg2: public PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const{
    coeff = 1./2;
    return true;
  }
};
struct PiPiProjectSolo: public PiPiProject{
  virtual bool operator()(std::complex<double> &coeff, const threeMomentum &mom) const{
    coeff = 1.;
    return true;
  }
};


//Compound filter that checks all the pion momenta involved actually exist
struct PiPiProjectAllowOnlyExistingPionMom: public PiPiProject{
  const PiPiProject &base;
  std::set<threeMomentum> pion_momenta;
  threeMomentum ptot;

  PiPiProjectAllowOnlyExistingPionMom(const PiPiProject &base, const threeMomentum &ptot, const std::vector<threeMomentum> pimomv): base(base), ptot(ptot){
    for(int i=0;i<pimomv.size();i++) pion_momenta.insert(pimomv[i]);
  }

  bool operator()(std::complex<double> &coeff, const threeMomentum &p1) const{
    bool c = base(coeff,p1);
    if(!c) return false;

    threeMomentum p2 = ptot - p1;

    return pion_momenta.count(p1) && pion_momenta.count(p2);
  }
};

PiPiProject* getProjector(const PiPiProjector p, const threeMomentum &ptot){
  switch(p){
  case PiPiProjector::A1:
    return (PiPiProject*)(new PiPiProjectA1);
  case PiPiProjector::Avg4:
    return (PiPiProject*)(new PiPiProjectAvg4);
  case PiPiProjector::Avg2:
    return (PiPiProject*)(new PiPiProjectAvg2);
  case PiPiProjector::Solo:
    return (PiPiProject*)(new PiPiProjectSolo);
  default:
    error_exit(std::cout << "getProjector unknown projector " << p << std::endl);
  }
}


//------------------------------------------ Momentum pair filters -------------------------------------------------------------

struct PiPiMomAllowAll: public PiPiMomAllow{
  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    coeff = 1.; return true;
  }
};
struct PiPiMomAllowOrig64: public PiPiMomAllow{
  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    if(abs(p1src[0])==1 && abs(p1src[1])==1 && abs(p1src[2])==1 &&
       abs(p1snk[0])==1 && abs(p1snk[1])==1 && abs(p1snk[2])==1){ 
      coeff = 1.; return true;
    }else{
      return false;
    }
  }
};


struct PiPiMomAllowParityAxisPermSymmReduced: public PiPiMomAllow{
  typedef std::pair<threeMomentum, threeMomentum> momPair;
  std::map<momPair, int> allow;

  PiPiMomAllowParityAxisPermSymmReduced(){
    //Reduced combination set obtained by utilizing parity and axis interchange symmetries
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,-1,1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,1,1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,-1,-1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,-1,1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 12;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,1}) ) ] = 12;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,-1}) ) ] = 6;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,1}) ) ] = 6;
  }    

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    std::map<momPair, int>::const_iterator it = allow.find(momPair(p1src,p1snk));
    if(it == allow.end()) return false;
    else{
      coeff = it->second;
      return true;
    }
  }
};

struct PiPiMomAllowAuxDiagSymmReduced: public PiPiMomAllow{
  typedef std::pair<threeMomentum, threeMomentum> momPair;
  std::map<momPair, int> allow;

  PiPiMomAllowAuxDiagSymmReduced(){
    //Reduced combination set obtained by utilizing auxiliary diagram symmetry
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 1;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({1,1,1}) ) ] = 1;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,-1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,1}) ) ] = 1;
    allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,-1}) ) ] = 1;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 1;
    allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({1,-1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,1,-1}), threeMomentum({1,-1,1}) ) ] = 1;
    allow[ momPair(threeMomentum({-1,1,-1}), threeMomentum({1,-1,-1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,1,-1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({-1,1,1}), threeMomentum({1,-1,-1}) ) ] = 1;
    allow[ momPair(threeMomentum({-1,1,1}), threeMomentum({-1,1,1}) ) ] = 2;
    allow[ momPair(threeMomentum({1,-1,-1}), threeMomentum({-1,1,1}) ) ] = 1;
  }    

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    std::map<momPair, int>::const_iterator it = allow.find(momPair(p1src,p1snk));
    if(it == allow.end()) return false;
    else{
      coeff = it->second;
      return true;
    }
  }
};

struct PiPiMomAllowAuxDiagParityAxisPermSymmReduced: public PiPiMomAllow{
  typedef std::pair<threeMomentum, threeMomentum> momPair;
  std::map<momPair, int> allow;

  PiPiMomAllowAuxDiagParityAxisPermSymmReduced(const threeMomentum &ptot){
    //Reduced combination set obtained by utilizing aux diag, parity and axis interchange symmetries
    if(ptot == threeMomentum({0,0,0})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,1,1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,1}) ) ] = 12;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,1,-1}) ) ] = 12;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,1}) ) ] = 6;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({1,1,-1}) ) ] = 6;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,1,-1}) ) ] = 12;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({1,-1,1}) ) ] = 12;
    }else if(ptot == threeMomentum({2,0,0})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,1}) ) ] = 4;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,1,1}) ) ] = 1;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,-1}) ) ] = 4;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,-1,-1}), threeMomentum({-1,-1,-1}) ) ] = 1;
    }else if(ptot == threeMomentum({0,2,0})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,1}) ) ] = 4;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,-1,1}) ) ] = 1;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,-1}) ) ] = 4;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({-1,1,-1}), threeMomentum({-1,-1,-1}) ) ] = 1;
    }else if(ptot == threeMomentum({0,0,2})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,1,-1}) ) ] = 1;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,1,-1}) ) ] = 4;
      allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,-1,-1}) ) ] = 1;
      allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 4;
      allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({1,-1,-1}) ) ] = 2;
    }else if(ptot == threeMomentum({2,2,0})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,1}) ) ] = 1;
      allow[ momPair(threeMomentum({1,1,-1}), threeMomentum({-1,-1,-1}) ) ] = 1;
    }else if(ptot == threeMomentum({2,0,2})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,1,-1}) ) ] = 1;
      allow[ momPair(threeMomentum({1,-1,1}), threeMomentum({-1,-1,-1}) ) ] = 1;
    }else if(ptot == threeMomentum({0,2,2})){
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({1,1,1}), threeMomentum({1,-1,-1}) ) ] = 1;
      allow[ momPair(threeMomentum({-1,1,1}), threeMomentum({-1,-1,-1}) ) ] = 1;
    }else if(ptot == threeMomentum({-2, 2, 0})){
      allow[ momPair(threeMomentum({-1,1,-1}), threeMomentum({1,-1,1}) ) ] = 2;
      allow[ momPair(threeMomentum({-1,1,-1}), threeMomentum({1,-1,-1}) ) ] = 2;
    }else if(ptot == threeMomentum({-2, 0, 2})){
      allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,-1,-1}) ) ] = 2;      
    }else if(ptot == threeMomentum({0, -2, 2})){
      allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({1,1,-1}) ) ] = 2;
      allow[ momPair(threeMomentum({-1,-1,1}), threeMomentum({-1,1,-1}) ) ] = 2;      
    }else{
      error_exit(std::cout << "PiPiMomAllowAuxDiagParityAxisPermSymmReduced unsupported total momentum " << ptot << std::endl);
    }
  }

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    std::map<momPair, int>::const_iterator it = allow.find(momPair(p1src,p1snk));
    if(it == allow.end()) return false;
    else{
      coeff = it->second;
      return true;
    }
  }
};

//Compound filter that checks all the pion momenta involved actually exist
struct PiPiMomAllowOnlyExistingPionMom: public PiPiMomAllow{
  const PiPiMomAllow &base;
  std::set<threeMomentum> pion_momenta;
  threeMomentum ptot;

  PiPiMomAllowOnlyExistingPionMom(const PiPiMomAllow &base, const threeMomentum &ptot, const std::vector<threeMomentum> pimomv): base(base), ptot(ptot){
    for(int i=0;i<pimomv.size();i++) pion_momenta.insert(pimomv[i]);
  }

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    bool c = base(coeff,p1src,p1snk);
    if(!c) return false;

    threeMomentum p2src = ptot - p1src;
    threeMomentum p2snk = -ptot - p1snk;

    return pion_momenta.count(p1src) && pion_momenta.count(p2src) && pion_momenta.count(p1snk) && pion_momenta.count(p2snk);
  }
};


//The factory function
PiPiMomAllow* getMomPairFilter(const PiPiMomAllowed p, const threeMomentum &ptot){
  switch(p){
  case PiPiMomAllowed::All:
    return (PiPiMomAllow*)(new PiPiMomAllowAll);
  case PiPiMomAllowed::Orig64:
    return (PiPiMomAllow*)(new PiPiMomAllowOrig64);
  case PiPiMomAllowed::ParityAxisPermSymmReduced:
    return (PiPiMomAllow*)(new PiPiMomAllowParityAxisPermSymmReduced);
  case PiPiMomAllowed::AuxDiagSymmReduced:
    return (PiPiMomAllow*)(new PiPiMomAllowAuxDiagSymmReduced);
  case PiPiMomAllowed::AuxDiagParityAxisPermSymmReduced:
    return (PiPiMomAllow*)(new PiPiMomAllowAuxDiagParityAxisPermSymmReduced(ptot));
  default:
    error_exit(std::cout << "getMomPairFilter unknown filter " << p << std::endl);
  }
}



#endif

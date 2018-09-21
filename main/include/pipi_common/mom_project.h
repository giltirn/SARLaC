#ifndef _PIPI_MOM_PROJECT_H_
#define _PIPI_MOM_PROJECT_H_

#include<set>
#include<regex>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#include<config.h>
#include<utils/macros.h>

#include "threemomentum.h"
#include "enums.h"

CPSFIT_START_NAMESPACE

struct PiPiProjectorBase{
  virtual int nMomenta() const = 0;
  virtual threeMomentum momentum(const int i) const = 0;   //momentum of pi1
  virtual std::complex<double> coefficient(const int i) const = 0;
  virtual ~PiPiProjectorBase(){}
};

struct PiPiProjectorA1Basis111: public PiPiProjectorBase{
  int nMomenta() const{ return 8; }
  threeMomentum momentum(const int i) const{ 
    static std::vector<threeMomentum> mom({  {1,1,1}, {-1,1,1}, {1,-1,1}, {1,1,-1},
					     {-1,-1,-1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1} });
    return mom[i];
  }
  std::complex<double> coefficient(const int i) const{ return 1./8; }
};

struct PiPiProjectorA1Basis311: public PiPiProjectorBase{
  int nMomenta() const{ return 24; }
  threeMomentum momentum(const int i) const{ 
    static std::vector<threeMomentum> mom({
	{3,1,1}, {-3,1,1}, {3,-1,1}, {3,1,-1}, {-3,-1,-1}, {3,-1,-1}, {-3,1,-1}, {-3,-1,1}, 
	{1,3,1}, {1,-3,1}, {1,3,-1}, {-1,3,1}, {-1,-3,-1}, {-1,3,-1}, {-1,-3,1}, {1,-3,-1}, 
	{1,1,3}, {1,1,-3}, {-1,1,3}, {1,-1,3}, {-1,-1,-3}, {-1,-1,3}, {1,-1,-3}, {-1,1,-3} });
    return mom[i];
  }
  std::complex<double> coefficient(const int i) const{ return 1./24; }
};



PiPiProjectorBase* getProjector(const PiPiProjector proj, const threeMomentum &ptot){
  switch(proj){
  case PiPiProjector::A1momSet111:
    return new PiPiProjectorA1Basis111;
  case PiPiProjector::A1momSet311:
    return new PiPiProjectorA1Basis311;
  default:
    error_exit(std::cout << "getProjector unknown projector " << proj << std::endl);
  }
}

//Return the list of momenta that appear in either the source or sink projector
std::vector<threeMomentum> getSrcSnkMomentumSet(const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk){
    std::set<threeMomentum> pset;
    for(int i=0;i<proj_src.nMomenta();i++) pset.insert(proj_src.momentum(i));
    for(int i=0;i<proj_snk.nMomenta();i++) pset.insert(proj_snk.momentum(i));
    std::vector<threeMomentum> psetv(pset.size());
    int j=0;
    for(auto it=pset.begin();it!=pset.end();it++) psetv[j++] = *it;
    return psetv;
}

CPSFIT_END_NAMESPACE

#endif

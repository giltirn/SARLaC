#ifndef _PIPI_MOM_PROJECT_H_
#define _PIPI_MOM_PROJECT_H_

#include<set>
#include<regex>
#include <unordered_set>
#include <boost/functional/hash.hpp>

#include<config.h>
#include<utils/macros.h>

#include<parser.h>
#include "threemomentum.h"

SARLAC_START_NAMESPACE

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

//Average over all p1 in the 111 basis for which  ptot - p1 is also in the 111 basis   energy  ~2*E_0
struct PiPiProjectorMovingSwaveGround: public PiPiProjectorBase{
  std::vector<threeMomentum> mom;

  PiPiProjectorMovingSwaveGround(const threeMomentum &ptot){
    static std::vector<threeMomentum> basis111({  {1,1,1}, {-1,1,1}, {1,-1,1}, {1,1,-1},
					       {-1,-1,-1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1} });
    for(int p1i=0; p1i<8; p1i++){
      threeMomentum p2 = ptot - basis111[p1i];
      if(std::find(basis111.begin(), basis111.end(), p2) != basis111.end()) 
	mom.push_back(basis111[p1i]);
    }
  }
  int nMomenta() const{ return mom.size(); }
  inline threeMomentum momentum(const int i) const{ 
    return mom[i];
  }
  std::complex<double> coefficient(const int i) const{ return 1./mom.size(); }
};

//Average over all p1 in the 311/111 basis for which  ptot - p1 is in the 111/311 basis   energy ~E_0 + E_1
struct PiPiProjectorMovingSwaveExc1: public PiPiProjectorBase{
  std::vector<threeMomentum> mom;

  PiPiProjectorMovingSwaveExc1(const threeMomentum &ptot){
    static std::vector<threeMomentum> basis111({  {1,1,1}, {-1,1,1}, {1,-1,1}, {1,1,-1},
					       {-1,-1,-1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1} });
    static std::vector<threeMomentum> basis311({
	{3,1,1}, {-3,1,1}, {3,-1,1}, {3,1,-1}, {-3,-1,-1}, {3,-1,-1}, {-3,1,-1}, {-3,-1,1}, 
	{1,3,1}, {1,-3,1}, {1,3,-1}, {-1,3,1}, {-1,-3,-1}, {-1,3,-1}, {-1,-3,1}, {1,-3,-1}, 
	{1,1,3}, {1,1,-3}, {-1,1,3}, {1,-1,3}, {-1,-1,-3}, {-1,-1,3}, {1,-1,-3}, {-1,1,-3} });

    for(int p1i=0; p1i<8; p1i++){
      threeMomentum p2 = ptot - basis111[p1i];
      if(std::find(basis311.begin(), basis311.end(), p2) != basis311.end()) 
	mom.push_back(basis111[p1i]);
    }
    for(int p1i=0; p1i<24; p1i++){
      threeMomentum p2 = ptot - basis311[p1i];
      if(std::find(basis111.begin(), basis111.end(), p2) != basis111.end()) 
	mom.push_back(basis311[p1i]);
    }
  }
  int nMomenta() const{ return mom.size(); }
  inline threeMomentum momentum(const int i) const{ 
    return mom[i];
  }
  std::complex<double> coefficient(const int i) const{ return 1./mom.size(); }
};

//Average over all p1 in the 311 basis for which  ptot - p1 is also in the 311 basis   energy ~2*E_1
struct PiPiProjectorMovingSwaveExc2: public PiPiProjectorBase{
  std::vector<threeMomentum> mom;

  PiPiProjectorMovingSwaveExc2(const threeMomentum &ptot){
    static std::vector<threeMomentum> basis311({
	{3,1,1}, {-3,1,1}, {3,-1,1}, {3,1,-1}, {-3,-1,-1}, {3,-1,-1}, {-3,1,-1}, {-3,-1,1}, 
	{1,3,1}, {1,-3,1}, {1,3,-1}, {-1,3,1}, {-1,-3,-1}, {-1,3,-1}, {-1,-3,1}, {1,-3,-1}, 
	{1,1,3}, {1,1,-3}, {-1,1,3}, {1,-1,3}, {-1,-1,-3}, {-1,-1,3}, {1,-1,-3}, {-1,1,-3} });

    for(int p1i=0; p1i<24; p1i++){
      threeMomentum p2 = ptot - basis311[p1i];
      if(std::find(basis311.begin(), basis311.end(), p2) != basis311.end()) 
	mom.push_back(basis311[p1i]);
    }
  }
  int nMomenta() const{ return mom.size(); }
  inline threeMomentum momentum(const int i) const{ 
    return mom[i];
  }
  std::complex<double> coefficient(const int i) const{ return 1./mom.size(); }
};

GENERATE_ENUM_AND_PARSER(PiPiProjector, (A1momSet111)(A1momSet311)(MovingSwaveGround)(MovingSwaveExc1)(MovingSwaveExc2) );

PiPiProjectorBase* getProjector(const PiPiProjector proj, const threeMomentum &ptot){
  switch(proj){
  case PiPiProjector::A1momSet111:
    if(ptot != threeMomentum({0,0,0})) error_exit(std::cout << "getProjector projector " << proj << " does not apply to moving pipi states\n");
    return new PiPiProjectorA1Basis111;
  case PiPiProjector::A1momSet311:
    if(ptot != threeMomentum({0,0,0})) error_exit(std::cout << "getProjector projector " << proj << " does not apply to moving pipi states\n");
    return new PiPiProjectorA1Basis311;
  case PiPiProjector::MovingSwaveGround:
    return new PiPiProjectorMovingSwaveGround(ptot);
  case PiPiProjector::MovingSwaveExc1:
    return new PiPiProjectorMovingSwaveExc1(ptot);
  case PiPiProjector::MovingSwaveExc2:
    return new PiPiProjectorMovingSwaveExc2(ptot);
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

SARLAC_END_NAMESPACE

#endif

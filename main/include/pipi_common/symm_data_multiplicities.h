#ifndef SYMM_DATA_MULTIPLICITIES_H_
#define SYMM_DATA_MULTIPLICITIES_H_

//For the extended analysis we have pion momenta from the set (+-1,+-1,+-1) and (+-3,+-1,+-1) + (perms), 32 in total. 
//This leads to a large number of potential pipi diagrams for a given total pipi momentum. In order to reduce this to something
//manageable we utilize the parity, auxiliary diagram and axis permutation symmetries to compute only a subset.

//In order to incorporate this in the analysis we create a filter with the appropriate multiplicities for use in the projection
#include<iostream>
#include<set>
#include<regex>
#include <unordered_set>
#include <unordered_map>
#include <boost/functional/hash.hpp>

#include<config.h>
#include<utils/macros.h>

#include "threemomentum.h"
#include "enums.h"
#include "mom_project.h"

CPSFIT_START_NAMESPACE

struct ConMomentum{
  threeMomentum pi1_src;
  threeMomentum pi1_snk;
  threeMomentum pi2_src;
  threeMomentum pi2_snk;

  threeMomentum p_tot;
  
  inline static int parseMomDir(std::string m){   
    if(m[0] == '_') m[0] = '-';
    int out;
    assert(sscanf(m.c_str(),"%d",&out) == 1);
    return out;
  }

  static threeMomentum parseMom(const std::string &mom){
    //std::cout << "parseMom " << mom << std::endl;
    std::smatch m;
    std::regex r(R"((_?\d+)(_?\d+)(_?\d+))");  
    assert(std::regex_search(mom,m,r));
    
    return threeMomentum({parseMomDir(m[1]),
	                  parseMomDir(m[2]),
			  parseMomDir(m[3])});    
  }
  
  void parse(const std::string &filename){ //should have an enum for selecting filename fmt
    std::smatch m;
    std::regex r( R"(_p1src((?:_?\d+){3})_p2src((?:_?\d+){3})_p1snk((?:_?\d+){3})_p2snk((?:_?\d+){3}))");  
    assert(std::regex_search(filename,m,r));

    /* std::cout << "parse " << filename << std::endl; //(?:(?:_\d+)|\d+){3} */
    /* for(int i=0;i<=4;i++) std::cout << "\"" << m[i] << "\"" << std::endl; */

    pi1_src = parseMom(m[1]);
    pi2_src = parseMom(m[2]);
    pi1_snk = parseMom(m[3]);
    pi2_snk = parseMom(m[4]);

    p_tot = pi1_src + pi2_src;
    assert( pi1_snk + pi2_snk == -p_tot);

    //Daiqian used conventions where the Fourier transforms into momentum space come are \sum_x f(x) exp(ipx) and not with exp(-ipx). 
    //Thus his file format (which I save in) has momentum with flipped sign relative to the actual momentum. Correct for that here
#define USE_DAIQIANS_CONVENTIONS  
#ifdef USE_DAIQIANS_CONVENTIONS
    pi1_src = -pi1_src;
    pi2_src = -pi2_src;
    pi1_snk = -pi1_snk;
    pi2_snk = -pi2_snk;
    p_tot = -p_tot;
#endif

    //std::cout << "Parsed " << filename << " : " << pi1_src << "  " << pi2_src << " " << pi1_snk << " " << pi2_snk << std::endl;
  }

  ConMomentum(): p_tot({0,0,0}){}
  ConMomentum(const std::string &filename){ parse(filename); }
  ConMomentum(const threeMomentum &pi1_src, const threeMomentum &pi1_snk, const threeMomentum &p_tot): pi1_src(pi1_src), pi2_src(p_tot-pi1_src),
												       pi1_snk(pi1_snk), pi2_snk(-p_tot-pi1_snk), p_tot(p_tot){}
      
  void applyParity(){
    pi1_src = -pi1_src;
    pi2_src = -pi2_src;
    pi1_snk = -pi1_snk;
    pi2_snk = -pi2_snk;
    p_tot = -p_tot;
  }
  void applyAxisPerm(const int i){ //i \in {0..5}
    threeMomentum *moms[5] = { &pi1_src, &pi2_src, &pi1_snk, &pi2_snk, &p_tot };   
    for(int pp=0;pp<5;pp++)
      *moms[pp] = axisPerm(i, *moms[pp]);
  }

  void applyAuxDiagSymm(){
    threeMomentum pi1_src_o(pi1_src), pi2_src_o(pi2_src),  pi1_snk_o(pi1_snk), pi2_snk_o(pi2_snk), p_tot_o(p_tot);

    // (pi1_src pi2_src pi2_snk pi1_snk) = (pi2_snk pi1_snk pi1_src pi2_src)
    
    pi1_src = pi2_snk_o;
    pi2_src = pi1_snk_o;
    pi1_snk = pi2_src_o;
    pi2_snk = pi1_src_o;

    p_tot = -p_tot_o;
  }
  
  bool operator==(const ConMomentum &r) const{
    return pi1_src == r.pi1_src && pi2_src == r.pi2_src && pi1_snk == r.pi1_snk && pi2_snk == r.pi2_snk && p_tot == r.p_tot;
  }
};

std::ostream & operator<<(std::ostream &os, const ConMomentum &p){ os << p.pi1_src << " " << p.pi2_src << " " << p.pi1_snk << " " << p.pi2_snk; return os; }

CPSFIT_END_NAMESPACE

namespace boost{
  std::size_t hash_value(const CPSfit::ConMomentum &r){
  std::size_t seed = 0;
  boost::hash_combine(seed, r.pi1_src);
  boost::hash_combine(seed, r.pi2_src);
  boost::hash_combine(seed, r.pi1_snk);
  boost::hash_combine(seed, r.pi2_snk);
  return seed;
}
};

CPSFIT_START_NAMESPACE

struct ConMomentumHasher{
  inline size_t operator()(const ConMomentum &p) const{ return boost::hash_value(p); }
};


struct AvailCorr{
  typedef std::unordered_set<ConMomentum, ConMomentumHasher> hashMapType;
  typename hashMapType::const_iterator it;

  //the following give its relation to the original ConMomentum
  int p; //parity (0,1)
  int a; //aux diag (0,1)
  int r; //axis perm (0..5)
    
  AvailCorr(typename hashMapType::const_iterator it, int p, int a, int r): it(it), p(p), a(a), r(r){}
};

struct PiPiMomSelectSymmetrySubset: public PiPiCorrelatorSelector{
  std::unordered_map<ConMomentum, double, ConMomentumHasher> multiplicities;
  threeMomentum p_tot;

  PiPiMomSelectSymmetrySubset(const threeMomentum &p_tot): p_tot(p_tot){};

  inline void incrementMultiplicity(const ConMomentum &cm, const double by){
    auto it = multiplicities.find(cm);
    if(it == multiplicities.end()) multiplicities[cm] = by;
    else it->second = it->second + by;
  }

  bool operator()(double &coeff, const threeMomentum &p1src, const threeMomentum &p1snk) const{
    ConMomentum cm(p1src, p1snk, p_tot);
    auto it = multiplicities.find(cm);
    if(it != multiplicities.end()){
      coeff = it->second;
      return true;
    }else{
      return false;
    }
  }
};

std::ostream & operator<<(std::ostream &os, const PiPiMomSelectSymmetrySubset &f){
  for(auto it = f.multiplicities.begin(); it != f.multiplicities.end(); it++)
    os << it->first << " : " << it->second << std::endl;
  return os;
}




struct PiPiSymmetrySubset{
  typedef std::unordered_set<ConMomentum, ConMomentumHasher> hashMapType;
  typedef std::map<threeMomentum, hashMapType> PtotMapType;
  PtotMapType corrs_avail;

  //Default regex is for new file format
  void findAvailableCorrs(const std::string &dir){
    std::vector<std::string> files = listFiles(dir, R"(traj_\d+_FigureC)");
    std::cout << "PiPiSymmetrySubset::findAvailableCorrs found " << files.size() << " files with provided format\n";

    for(int i=0;i<files.size();i++){
      ConMomentum momentum(files[i]);
      corrs_avail[momentum.p_tot].insert(momentum);
    }
    for(auto it = corrs_avail.begin(); it != corrs_avail.end(); it++)
      std::cout << it->second.size() << " correlators with total momentum " << it->first << std::endl;
  }

  //Find a partner for the momentum p that exists in the set of data available. We do not allow transformations that change the overall total momentum
  AvailCorr findPartnerInAvailableCorrs(const ConMomentum &cmom) const{
    typename PtotMapType::const_iterator ptot_it = corrs_avail.find(cmom.p_tot);
    assert(ptot_it != corrs_avail.end());

    for(int par = 0; par < 2*2*6; par++){
      int rem = par;
      int p = rem % 2; rem/=2; //parity
      int a = rem % 2; rem/=2; //aux symm
      int r = rem; //axis perms
      
      ConMomentum cmom_new = cmom;
      if(p) cmom_new.applyParity();
      if(a) cmom_new.applyAuxDiagSymm();
      cmom_new.applyAxisPerm(r);

      typename hashMapType::const_iterator corr_it = ptot_it->second.find(cmom_new);
      
      if(corr_it != ptot_it->second.end()){ return AvailCorr(corr_it, p, a, r); }
    }
    error_exit(std::cout << "Could not find partner for correlator " << cmom << " among those available\n");
  }

  //set    0 : (+-1,+-1,+-1)   1: (+-3,+-1,+-1) + perms   2: all 32 momenta
  PiPiMomSelectSymmetrySubset createSelector(const std::vector<threeMomentum> &pimom, const threeMomentum &p_tot, const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow) const{
    int np = pimom.size();

    struct threeMomentumHasher{
      size_t operator()(const threeMomentum &p) const{ return boost::hash_value(p); }
    };    
    std::unordered_set<threeMomentum, threeMomentumHasher> pimoms; 
    for(int i=0; i<np; i++) pimoms.insert(pimom[i]);
    
    PiPiMomSelectSymmetrySubset out(p_tot);
    
    std::complex<double> csrc, csnk;
    double m;

    for(int i=0;i<np;i++){
      if(!proj_src(csrc, pimom[i])) continue;
      for(int j=0;j<np;j++){
	if(!proj_snk(csnk, pimom[j])) continue;

	if(!allow(m,pimom[i],pimom[j])) continue;

	ConMomentum want(pimom[i],  pimom[j], p_tot);
	if(pimoms.count(want.pi2_src) == 0 || pimoms.count(want.pi2_snk) == 0) continue; //all momenta must be in the set 
	AvailCorr avail = findPartnerInAvailableCorrs(want);
	
	std::cout << "Found match to " << want << " : " << *avail.it << " by symms p=" << avail.p << " a=" << avail.a << " r=" << avail.r << std::endl;
	out.incrementMultiplicity(*avail.it, std::real(csrc * csnk)*m );
      }
    }
    
    return out;
  }
  
  //Factory version
  inline PiPiMomSelectSymmetrySubset createSelector(const std::vector<threeMomentum> &pimom, const threeMomentum &p_tot, const PiPiProjector proj_src, const PiPiProjector proj_snk, const PiPiMomAllowed allow){
    std::unique_ptr<PiPiProject> prjsrc( getProjector(proj_src, p_tot) );
    std::unique_ptr<PiPiProject> prjsnk( getProjector(proj_snk, -p_tot) );
    std::unique_ptr<PiPiMomAllow> allw( getMomPairFilter(allow, p_tot) );
    return createSelector(pimom, p_tot, *prjsrc, *prjsnk, *allw);
  }

  PiPiSymmetrySubset(){}
  PiPiSymmetrySubset(const std::string &dir){
    findAvailableCorrs(dir);
  }
       

  //Test the code for all the momenta used in the extended measurement proposal
  void test() const{
    std::vector<threeMomentum> ptot = {
      {0,0,0}, 
      {2,0,0}, {0,2,0}, {0,0,2},
      {2,2,0}, {2,0,2}, {0,2,2},
      {-2,2,0}, {-2,0,2}, {0,-2,2},
      {2,2,2}, 
      {-2,2,2}, {2,-2,2}, {2,2,-2}
    };

    PiPiProjectSolo proj_all;
    PiPiMomAllowAll allow_all;

    for(int set=0;set<2;set++){
      std::vector<threeMomentum> pimom;
      if(set == 0 || set == 2){
	pimom = { {1,1,1}, {-1,1,1}, {1,-1,1}, {1,1,-1},
		  {-1,-1,-1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1} };
      }
      if(set == 1 || set == 2){
	std::vector<threeMomentum> pimom_base = { {3,1,1}, {-3,1,1}, {3,-1,1}, {3,1,-1},
						  {-3,-1,-1}, {3,-1,-1}, {-3,1,-1}, {-3,-1,1} };
	for(int n=0;n<3;n++) //cyclic permutation amount
	  for(int i=0;i<8;i++)
	    pimom.push_back(cyclicPermute(n, pimom_base[i]));      
      }

      for(int p=0;p<ptot.size();p++){
	std::cout << "Doing set " << set << " ptot " << ptot[p] << std::endl;
	auto f = createSelector(pimom, ptot[p], proj_all, proj_all, allow_all);
	std::cout << f << std::endl;
      }
    }
  }


};

CPSFIT_END_NAMESPACE

#endif

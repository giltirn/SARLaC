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
#include<utils/utils/string.h>

#include "../../correlator_utils.h"
#include "../../enums.h"

CPSFIT_START_NAMESPACE

//expect substrings  <TRAJ> <FIG> <TSEP_PIPI> <P1SRC> <P1SNK>   and optionally <P2SRC> <P2SNK>
const std::vector<subStringSpecify> & pipiFileFormatKeys(){ 
#define F(STR) subStringSpecify(STR)
#define FO(STR) subStringSpecify(STR,true)
  static std::vector<subStringSpecify> find = { F("<TRAJ>"), F("<FIG>"), F("<TSEP_PIPI>"), F("<P1SRC>"), F("<P1SNK>"),
						FO("<P2SRC>"),FO("<P2SNK>") };
#undef F
#undef FO
  return find;
}


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

  //Original format assumes ptot == 0
  void parseOrigFmt(const std::string &filename){ 
    std::smatch m;
    std::regex r( R"(_mom((?:_?\d+){3})_mom((?:_?\d+){3}))");  
    assert(std::regex_search(filename,m,r));

    pi1_src = parseMom(m[1]);
    pi2_src = -pi1_src;
    pi1_snk = parseMom(m[2]);
    pi2_snk = -pi2_snk;
    
    p_tot = threeMomentum({0,0,0});
    
    //Daiqian used conventions where the Fourier transforms into momentum space come are \sum_x f(x) exp(ipx) and not with exp(-ipx). 
    //Thus his file format (which I save in) has momentum with flipped sign relative to the actual momentum. Correct for that here
#ifdef USE_DAIQIANS_CONVENTIONS
    pi1_src = -pi1_src;
    pi2_src = -pi2_src;
    pi1_snk = -pi1_snk;
    pi2_snk = -pi2_snk;
    p_tot = -p_tot;
#endif
  }
    
  void parseGeneric(const std::string &filename, const std::string &fmt){ 
    subStringReplace matcher(fmt, pipiFileFormatKeys());

    std::map<std::string,std::string> smatch;
    if(!matcher.match(smatch, filename)) error_exit(std::cout << "Could not match filename " << filename << " to format string " << fmt << std::endl);

    pi1_src = parseMom(smatch["<P1SRC>"]);
    pi1_snk = parseMom(smatch["<P1SNK>"]);
    
    bool foundp2[2] = {false,false};
    std::map<std::string,std::string>::const_iterator it;
    if( (it = smatch.find("<P2SRC>")) != smatch.end()){ 
      pi2_src = parseMom(it->second); foundp2[0] = true;
    }
    if( (it = smatch.find("<P2SNK>")) != smatch.end()){ 
      pi2_snk = parseMom(it->second); foundp2[1] = true;
    }

    if(foundp2[0] && foundp2[1]){
      p_tot = pi1_src + pi2_src;
      if(pi1_snk + pi2_snk != -p_tot) error_exit(std::cout << "Filename " << filename << " source and sink total momentum don't match!\n");
    }else if(foundp2[0] && !foundp2[1]){
      p_tot = pi1_src + pi2_src;
      pi2_snk = -p_tot - pi1_snk;
    }else if(!foundp2[0] && foundp2[1]){
      p_tot = -pi1_snk - pi2_snk;
      pi2_src = p_tot - pi1_src;
    }else{
      //Have to assume zero total momentum (eg old file format)
      pi2_src = -pi1_src;
      pi2_snk = -pi2_snk;
      p_tot = threeMomentum({0,0,0});
    }
    
    //Daiqian used conventions where the Fourier transforms into momentum space come are \sum_x f(x) exp(ipx) and not with exp(-ipx). 
    //Thus his file format (which I save in) has momentum with flipped sign relative to the actual momentum. Correct for that here
#ifdef USE_DAIQIANS_CONVENTIONS
    pi1_src = -pi1_src;
    pi2_src = -pi2_src;
    pi1_snk = -pi1_snk;
    pi2_snk = -pi2_snk;
    p_tot = -p_tot;
#endif
  }

  ConMomentum(): p_tot({0,0,0}){}
  ConMomentum(const std::string &filename){ parse(filename); }
  ConMomentum(const std::string &filename, const std::string &fmt){ parseGeneric(filename,fmt); }
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




struct PiPiSymmetrySubset{
  typedef std::unordered_set<ConMomentum, ConMomentumHasher> hashMapType;
  typedef std::map<threeMomentum, hashMapType> PtotMapType;
  PtotMapType corrs_avail;

  void findAvailableCorrs(const std::string &dir, const std::string &file_fmt, const int traj_start, const int tsep_pipi,
			  const std::vector<threeMomentum> &p_pi, const std::vector<threeMomentum> &p_tot, const MomentumUnit fn_mom_unit = MomentumUnit::PiOverL){
    int pmult = fn_mom_unit == MomentumUnit::PiOverTwoL ? 2 : 1; //allow for pi/L (default) or pi/2L basis in file names. For these pmult = 1 and 2, respectively

    std::cout << "Finding available corrs in directory " << dir << std::endl;
    std::cout << "Searching for total momenta in set " << p_tot << std::endl;
    std::cout << "p1 momentum set is size " << p_pi.size() << " and has contents " << p_pi << std::endl;
    
    subStringReplace repl(file_fmt, pipiFileFormatKeys());
    //<TRAJ> <FIG> <TSEP_PIPI> <P1SRC> <P1SNK>   and optionally <P2SRC> <P2SNK>

    std::vector<std::string> wc = { anyToStr<int>(traj_start), "C", anyToStr<int>(tsep_pipi), "", "", "", "" };
    std::string &p1src_s = wc[3];
    std::string &p1snk_s = wc[4];
    std::string &p2src_s = wc[5];
    std::string &p2snk_s = wc[6];

    //Allow for parity and axis permutation that maintains the same total momentum
    std::set<threeMomentum> p_pi_trans; 
    for(int p=0;p<p_pi.size();p++){
      for(int par=0;par<2;par++){
	threeMomentum mom = par == 0 ? p_pi[p] : -p_pi[p];
	for(int perm=0; perm<6; perm++){
	  p_pi_trans.insert( axisPerm(perm, mom) );
	}
      }
    }

    for(int pt=0; pt< p_tot.size(); pt++){
      std::cout << "Looking for p_tot = " << p_tot[pt] << std::endl;

      for(auto psrc = p_pi_trans.begin(); psrc != p_pi_trans.end(); psrc++){	  
	p1src_s = momStr( (*psrc)*pmult );
	p2src_s = momStr( (p_tot[pt] - (*psrc) )*pmult );

	for(auto psnk = p_pi_trans.begin(); psnk != p_pi_trans.end(); psnk++){	  
	  p1snk_s = momStr( (*psnk)*pmult );
	  p2snk_s = momStr( (-p_tot[pt] - (*psnk) )*pmult );
	      
	  std::ostringstream fn; fn << dir << "/";
	  repl.replace(fn, wc);
	  
	  std::cout << "Checking for " << fn.str();
	      
	  if(fileExists(fn.str())){
	    std::cout << "... found" << std::endl;
	    ConMomentum momentum(*psrc, *psnk, p_tot[pt]);
	    corrs_avail[p_tot[pt]].insert(momentum);
	  }else{
	    std::cout << "... not found" << std::endl;
	  }
	}
      }
    }
    for(auto it = corrs_avail.begin(); it != corrs_avail.end(); it++)
      std::cout << it->second.size() << " correlators with total momentum " << it->first << std::endl;

    if(corrs_avail.size() == 0) error_exit(std::cout << "PiPiSymmetrySubset::findAvailableCorrs found no files in dir " << dir << " matching format " << file_fmt << " with traj = " << traj_start << ", tsep_pipi = " << tsep_pipi << " and fig = C\n");

  }

  //Find a partner for the momentum p that exists in the set of data available.
  //We do not allow transformations that change the overall total momentum unless allow_ptot_parity=true whereby it allows transformations that flip the overall momentum
  AvailCorr findPartnerInAvailableCorrs(const ConMomentum &cmom, bool allow_ptot_parity = false) const{
    typename PtotMapType::const_iterator ptot_it = corrs_avail.find(cmom.p_tot);
    if(ptot_it == corrs_avail.end()){
      if(allow_ptot_parity){
	ptot_it = corrs_avail.find(-cmom.p_tot);
	if(ptot_it == corrs_avail.end()) error_exit(std::cout << "PiPiSymmetrySubset::findPartnerInAvailableCorrs could not find total momentum " << cmom.p_tot << " or its parity partner in list\n");
      }else error_exit(std::cout << "PiPiSymmetrySubset::findPartnerInAvailableCorrs could not find total momentum " << cmom.p_tot << " in list\n");
    }

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
    error_exit(std::cout << "PiPiSymmetrySubset::findPartnerInAvailableCorrs Could not find partner for correlator " << cmom << " among those available\n");
  }

  PiPiSymmetrySubset(){}

  PiPiSymmetrySubset(const std::string &dir, const std::string &file_fmt, const int traj_start, const int tsep_pipi,
		     const std::vector<threeMomentum> &p_pi, const std::vector<threeMomentum> &p_tot, 
		     const MomentumUnit fn_mom_unit = MomentumUnit::PiOverL){
    findAvailableCorrs(dir,file_fmt,traj_start,tsep_pipi,p_pi,p_tot,fn_mom_unit);
  }
};


CPSFIT_END_NAMESPACE

#endif

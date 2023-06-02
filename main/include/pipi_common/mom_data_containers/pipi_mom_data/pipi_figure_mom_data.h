#pragma once

#include<config.h>
#include<utils/macros.h>

#include "../pipi_figure_mom_data_container.h"
#include "../../correlator_utils/mom_project.h"
#include "symm_data_multiplicities.h"

CPSFIT_START_NAMESPACE

////////////////////////
//FILENAME POLICIES
///////////////////////

class figureFilenamePolicyGeneric{
  subStringReplace repl; //expect substrings  <TRAJ> <FIG> <TSEP_PIPI> <P1SRC> <P1SNK>   and optionally <P2SRC> <P2SNK>
  const threeMomentum p_tot;
public:  

  figureFilenamePolicyGeneric(const std::string &fmt, const threeMomentum &p_tot): p_tot(p_tot){
#define F(STR) subStringSpecify(STR)
#define FO(STR) subStringSpecify(STR,true)

    static std::vector<subStringSpecify> find = { F("<TRAJ>"), F("<FIG>"), F("<TSEP_PIPI>"), F("<P1SRC>"), F("<P1SNK>"),
						  FO("<P2SRC>"),FO("<P2SNK>") };
    repl.chunkString(fmt, find);

#undef F
#undef FO
  }


  inline std::string operator()(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi) const{ 
    std::vector<std::string> with = { anyToStr(traj), std::string(1,fig), anyToStr(tsep_pipi), momStr(psrc), momStr(psnk), momStr(p_tot - psrc), momStr(-p_tot-psnk) };
    std::ostringstream os;
    os << data_dir << '/';
    repl.replace(os,with);
    return os.str();
  }
}; 

std::string figureFile(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi, const bool use_symmetric_quark_momenta){
  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_Figure" << fig << "_sep" << tsep_pipi << "_mom" << momStr(psrc) << "_mom" << momStr(psnk);
  if(use_symmetric_quark_momenta) os << "_symm";
  return os.str();
}

struct readFigureStationaryPolicy{
  bool use_symmetric_quark_momenta;
  readFigureStationaryPolicy(bool use_symmetric_quark_momenta): use_symmetric_quark_momenta(use_symmetric_quark_momenta){}
  
  inline std::string operator()(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi) const{ 
    return figureFile(data_dir,fig,traj,psnk,psrc,tsep_pipi,use_symmetric_quark_momenta);
  }
};


std::string figureFileTianleComoving(const std::string &data_dir, const char fig, const int traj, 
				     threeMomentum psnk1, threeMomentum psrc1, threeMomentum ptot,
				     const int tsep_pipi){
  //Tianle's momenta are in units of pi/2L
  psnk1 = psnk1 * 2;
  psrc1 = psrc1 * 2;
  ptot = ptot * 2;

  threeMomentum psrc2 = ptot - psrc1;
  threeMomentum psnk2 = -ptot - psnk1;

  //traj_700_FigureC_sep4_comove_pipi2pt_srcmom1_2_22_srcmom2_2_2_2_snkmom122_2_snkmom2222_v2

  std::ostringstream os;
  os << data_dir << "/traj_" << traj << "_Figure" << fig << "_sep" << tsep_pipi << "_comove_pipi2pt" 
     << "_srcmom1" << momStr(psrc1)
     << "_srcmom2" << momStr(psrc2)
     << "_snkmom1" << momStr(psnk1)
     << "_snkmom2" << momStr(psnk2)
     << "_v2";

  return os.str();
}
struct readFigureTianleComovingPolicy{
  threeMomentum ptot;
  readFigureTianleComovingPolicy(const threeMomentum &ptot): ptot(ptot){}
  
  inline std::string operator()(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi) const{ 
    return figureFileTianleComoving(data_dir,fig,traj,psnk,psrc,ptot,tsep_pipi);
  }
};

//This filename policy uses various symmetries to map particular momentum combinations to other data so as to reduce the total data volume
class PiPiSymmetrySubsetFigureFileMapping: public PiPiSymmetrySubset{
  subStringReplace repl; //expect substrings  <TRAJ> <FIG> <TSEP_PIPI> <P1SRC> <P1SNK>   and optionally <P2SRC> <P2SNK>
  threeMomentum p_tot;
  int pmult; //allow for pi/L (default) or pi/2L basis in file names. For these pmult = 1 and 2, respectively
  bool allow_ptot_parity;

  inline std::string filename(const std::string &data_dir, const char fig, const int traj, 
			      const threeMomentum &p1src, const threeMomentum &p2src, 
			      const threeMomentum &p1snk, const threeMomentum &p2snk, 
			      const int tsep_pipi) const{ 
    std::vector<std::string> with = { anyToStr(traj), std::string(1,fig), anyToStr(tsep_pipi), 
				      momStr(p1src*pmult), momStr(p1snk*pmult), 
				      momStr(p2src*pmult), momStr(p2snk*pmult) 
    };
    std::ostringstream os;
    os << data_dir << '/';
    repl.replace(os,with);
    return os.str();
  }
public:  
  
  //Auxiliary diagram and parity flip source total momentum
  static inline std::vector<threeMomentum> getFileSearchMomenta(const threeMomentum &ptot, bool allow_ptot_parity){
    if(!allow_ptot_parity || ptot == threeMomentum({0,0,0})) return std::vector<threeMomentum>({ptot});
    else return std::vector<threeMomentum>({ptot, -ptot});
  }
  
  PiPiSymmetrySubsetFigureFileMapping(const std::string &dir, const std::string &file_fmt, const int traj_start, const int tsep_pipi,
				      const std::vector<threeMomentum> &p_pi, const threeMomentum &p_tot, 
				      const MomentumUnit fn_mom_unit = MomentumUnit::PiOverL, bool allow_ptot_parity = false): 
    PiPiSymmetrySubset(dir, file_fmt, traj_start, tsep_pipi, p_pi, 
		       getFileSearchMomenta(p_tot, allow_ptot_parity), fn_mom_unit), pmult(fn_mom_unit == MomentumUnit::PiOverTwoL ? 2 : 1),
    p_tot(p_tot), allow_ptot_parity(allow_ptot_parity){

#define F(STR) subStringSpecify(STR)
#define FO(STR) subStringSpecify(STR,true)

    static std::vector<subStringSpecify> find = { F("<TRAJ>"), F("<FIG>"), F("<TSEP_PIPI>"), F("<P1SRC>"), F("<P1SNK>"),
						  FO("<P2SRC>"),FO("<P2SNK>") };
    repl.chunkString(file_fmt, find);

#undef F
#undef FO
  }
																	   
  inline std::string operator()(const std::string &data_dir, const char fig, const int traj, const threeMomentum &psnk, const threeMomentum &psrc, const int tsep_pipi) const{ 
    ConMomentum want(psrc,  psnk, p_tot);
    AvailCorr avail = findPartnerInAvailableCorrs(want, allow_ptot_parity);
    const ConMomentum &found = *avail.it;
    std::cout << "Found match to " << want << " : " << found << " by symms p=" << avail.p << " a=" << avail.a << " r=" << avail.r << std::endl;
    return filename(data_dir, fig, traj, found.pi1_src, found.pi2_src, found.pi1_snk, found.pi2_snk, tsep_pipi);
  }

};


/////////////////////////////
//READ POLICY
//Takes a range and directory with an associated filename policy and returns a filename indexed by a sample index 0...nsample-1
////////////////////////////

template<typename FilenamePolicy>
struct PiPiFigureBasicReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  int tsep_pipi;
  const FilenamePolicy &fn;
  
  PiPiFigureBasicReadPolicy(const FilenamePolicy &fn, const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi):
    fn(fn), traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir), tsep_pipi(tsep_pipi){
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const char fig, const threeMomentum &psnk, const threeMomentum &psrc) const{
    return fn(dir, fig, traj_start + sample * traj_inc, psnk, psrc, tsep_pipi);
  }
};


///////////////////////////////////////////////
//READ PIPI FIGUREDATA FOR ONE OR ALL MOMENTA
///////////////////////////////////////////////

//For a given figure (C,D,R) read the raw data. Projectors are used to only load data that is required
template<typename ReadPolicy>
void readPiPi2ptFigure(figureDataAllMomenta &raw_data, const char fig, const int Lt, 
		       const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk, const ReadPolicy &rp){
  std::cout << "Reading figure " << fig << "\n"; boost::timer::auto_cpu_timer t(std::string("Report: Read figure ") + fig + " in %w s\n");
  int nsample = rp.nsample();

  raw_data.setup(Lt);

  typedef std::map< std::string, std::vector<double> > CacheType;

  int nmompairs = proj_snk.nMomenta() * proj_src.nMomenta();
  figureData* dest_ptrs[nmompairs];
  std::pair<threeMomentum,threeMomentum> dest_mom[nmompairs];
  
  //Initialize containers
  int i=0;
  for(int psnki=0; psnki<proj_snk.nMomenta();psnki++){
    threeMomentum psnk = proj_snk.momentum(psnki);

    for(int psrci=0; psrci<proj_src.nMomenta();psrci++){
      threeMomentum psrc = proj_src.momentum(psrci);

      figureData &into = raw_data(fig, momComb(psnk, psrc));      
      into.initializeElements(rawDataDistributionD(nsample));
      dest_ptrs[i] = &into;
      dest_mom[i] = { psnk, psrc };
      i++;
    }
  }

  //Read in parallel
#pragma omp parallel for
  for(int sample=0; sample < nsample; sample++){
    CacheType sample_cache; //cache data over different momentum combinations as for mapping policy this can be reused
    for(int mp=0; mp<nmompairs; mp++){
      std::string filename = rp.filename(sample, fig, dest_mom[mp].first, dest_mom[mp].second);
      std::cout << "Thread " << omp_get_thread_num() << " parsing " << filename << std::endl;
      dest_ptrs[mp]->parseCDR(filename, sample, &sample_cache);
    }
  }
}

//Call the above with the standard read policy
template<typename FilenamePolicy>
inline void readPiPi2ptFigure(figureDataAllMomenta &raw_data, const char fig, const std::string &data_dir, const int tsep_pipi, const int Lt,
			      const int traj_start, const int traj_inc, const int traj_lessthan, 
			      const PiPiProjectorBase &proj_src, const PiPiProjectorBase &proj_snk, 
			      const FilenamePolicy &fn){
  PiPiFigureBasicReadPolicy<FilenamePolicy> rp(fn, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  readPiPi2ptFigure(raw_data, fig, Lt, proj_src, proj_snk, rp);
}


CPSFIT_END_NAMESPACE

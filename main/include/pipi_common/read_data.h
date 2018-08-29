#ifndef _PIPI_READ_DATA_H_
#define _PIPI_READ_DATA_H_

#include<algorithm>

#include<config.h>
#include<utils/macros.h>

#include "mom_data_containers.h"
#include "mom_project.h"

CPSFIT_START_NAMESPACE

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


class bubbleFilenamePolicyGeneric{
  subStringReplace repl; //expect <TRAJ> <TSEP_PIPI> <PB> and optional <PA>.  Here the pions A and B are the earlier and later pions, respectively
  SourceOrSink src_snk; //whether a source or sink pipi  
  const threeMomentum p_tot;

public:
  bubbleFilenamePolicyGeneric(const std::string &fmt, const threeMomentum &p_tot, const SourceOrSink src_snk): p_tot(p_tot), src_snk(src_snk){
#define F(STR) subStringSpecify(STR)
#define FO(STR) subStringSpecify(STR,true)

    static std::vector<subStringSpecify> find = { F("<TRAJ>"), F("<TSEP_PIPI>"), F("<PB>"), FO("<PA>") };
    repl.chunkString(fmt, find);

#undef F
#undef FO
  }
  
  inline std::string operator()(const std::string &data_dir, const int traj, const threeMomentum &p1, const int tsep_pipi) const{ 
    //Note p1 here is in the standard format, where is is the later pion at source and earlier pion at sink. 
    threeMomentum pA, pB;
    if(src_snk == Source){
      //pB = p1   pA = p2
      pA = p_tot - p1;
      pB = p1;
    }else{ //Sink
      //pB = p2   pA = p1
      pA = p1;
      pB = -p_tot - p1;
    }
    std::vector<std::string> with =   { anyToStr(traj), anyToStr(tsep_pipi), momStr(pB), momStr(pA) };

    std::ostringstream os;
    os << data_dir << '/';
    repl.replace(os, with);
    
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

//For a given figure (C,D,R) read the raw data. Selector is used to only load data that is required
template<typename FilenamePolicy>
void readFigure(figureDataAllMomenta &raw_data, const char fig, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, const FilenamePolicy &fn, const std::vector<threeMomentum> &pion_momenta,
		const PiPiCorrelatorSelector &corr_select){
  std::cout << "Reading figure " << fig << "\n"; boost::timer::auto_cpu_timer t(std::string("Report: Read figure ") + fig + " in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;

  raw_data.setup(Lt,nsample);

  const int nmom = pion_momenta.size();
  
  double dummy;

  for(int psnk=0;psnk<nmom;psnk++){
    for(int psrc=0;psrc<nmom;psrc++){
      if(!corr_select(dummy,pion_momenta[psrc],pion_momenta[psnk])) continue;

      figureData &into = raw_data(fig, momComb(pion_momenta[psnk], pion_momenta[psrc]));      
      
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	int traj = traj_start + sample * traj_inc;
	std::string filename = fn(data_dir, fig, traj, pion_momenta[psnk], pion_momenta[psrc], tsep_pipi);
	std::cout << "Parsing " << filename << std::endl;
	into.parseCDR(filename, sample);
      }
    }
  }
}


struct readBubbleStationaryPolicy{
  SourceOrSink src_snk; //whether a source or sink pipi
  bool use_symmetric_quark_momenta;
  readBubbleStationaryPolicy(bool use_symmetric_quark_momenta, const SourceOrSink src_snk): use_symmetric_quark_momenta(use_symmetric_quark_momenta),src_snk(src_snk){}
  
  inline std::string operator()(const std::string &data_dir, const int traj, const threeMomentum &p1, const int tsep_pipi) const{ 
    std::ostringstream os;
    os << data_dir << "/traj_" << traj << "_FigureVdis_sep" << tsep_pipi << "_mom" << momStr(src_snk == Sink ? -p1 : p1); //file momentum is always for the later pion, which is p2 for the sink vertex
    if(use_symmetric_quark_momenta) os << "_symm";
    return os.str();
  }
};

struct readBubbleTianleComovingPolicy{
  SourceOrSink src_snk; //whether a source or sink pipi
  threeMomentum ptot; //total momentum of this pipi (note factor of 2 in init is because Tianle uses units of pi/2L)

  readBubbleTianleComovingPolicy(const threeMomentum &ptot, const SourceOrSink src_snk): ptot(ptot*2),src_snk(src_snk){}
  
  inline std::string operator()(const std::string &data_dir, const int traj, threeMomentum p1, const int tsep_pipi) const{ 
    //Tianle's momenta are in units of pi/2L
    p1 = p1*2;
    threeMomentum p2 = ptot - p1;

    // In Tianle's conventions the earlier pion is name "pi_src" and the latter "pi_snk"
    const threeMomentum &p_pi_src = src_snk == Sink ? p1 : p2;  //for sink pipi the earlier timeslice is p1
    const threeMomentum &p_pi_snk = src_snk == Sink ? p2 : p1;

    //traj_700_FigureVdis_sep4_comove_pipi2pt_srcmom2_22_snkmom2_2_2_v2
    std::ostringstream os;
    os << data_dir << "/traj_" << traj << "_FigureVdis_sep" << tsep_pipi << "_comove_pipi2pt" 
       << "_srcmom" << momStr(p_pi_src) << "_snkmom" << momStr(p_pi_snk) << "_v2";

    return os.str();
  }
};

template<typename FilenamePolicy>
struct PiPiBubbleBasicReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  int tsep_pipi;
  const FilenamePolicy &fn;
  
  PiPiBubbleBasicReadPolicy(const FilenamePolicy &fn, const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi):
    fn(fn), traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir), tsep_pipi(tsep_pipi){
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const threeMomentum &p_pi1) const{
    return fn(dir, traj_start + sample * traj_inc, p_pi1, tsep_pipi);
  }
};


template<typename ReadPolicy>
void readBubble(bubbleDataAllMomenta &raw_data, const int Lt, const ReadPolicy &rp, const std::vector<threeMomentum> &pion_momenta, const SourceOrSink src_snk){    
  const int nmom = pion_momenta.size();
  const int nsample = rp.nsample();

  for(int p=0;p<nmom;p++){
    bubbleData &into = raw_data(src_snk,pion_momenta[p]);
    
#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      std::string filename = rp.filename(sample,pion_momenta[p]);
      std::cout << "Parsing " << filename << std::endl;
      into.parse(filename, sample);
    }
  }
}

template<typename FilenamePolicy>
inline void readBubble(bubbleDataAllMomenta &raw_data, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, const FilenamePolicy &fn, 
		const std::vector<threeMomentum> &pion_momenta, const SourceOrSink src_snk){   
  PiPiBubbleBasicReadPolicy<FilenamePolicy> rp(fn, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  readBubble(raw_data, Lt, rp, pion_momenta, src_snk);
}

template<typename FilenamePolicy>
void readBubble(bubbleDataAllMomenta &raw_data, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, 
		const FilenamePolicy &fn_src, const FilenamePolicy &fn_snk,
		const std::vector<threeMomentum> &pion_momenta,
		const PiPiCorrelatorSelector &corr_select){
  std::cout << "Reading bubble\n"; boost::timer::auto_cpu_timer t("Report: Read bubble in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  raw_data.setup(Lt,tsep_pipi,nsample);
  
  std::vector<threeMomentum> src_mom_need;
  std::vector<threeMomentum> snk_mom_need;
  double m;
  for(int i=0;i<pion_momenta.size();i++){
    for(int j=0;j<pion_momenta.size();j++){
      if(corr_select(m, pion_momenta[i],pion_momenta[j])){
	if(std::find(src_mom_need.begin(), src_mom_need.end(), pion_momenta[i]) == src_mom_need.end()) src_mom_need.push_back(pion_momenta[i]);
	if(std::find(snk_mom_need.begin(), snk_mom_need.end(), pion_momenta[j]) == snk_mom_need.end()) snk_mom_need.push_back(pion_momenta[j]);
      }
    }
  }

  readBubble(raw_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, fn_src, src_mom_need, Source);
  readBubble(raw_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, fn_snk, snk_mom_need, Sink);
}

CPSFIT_END_NAMESPACE

#endif

#ifndef _PIPI_READ_DATA_H_
#define _PIPI_READ_DATA_H_

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


template<typename FilenamePolicy>
void readFigure(figureDataAllMomenta &raw_data, const char fig, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, const FilenamePolicy &fn, const std::vector<threeMomentum> &pion_momenta,
		const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow){
  std::cout << "Reading figure " << fig << "\n"; boost::timer::auto_cpu_timer t(std::string("Report: Read figure ") + fig + " in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;

  raw_data.setup(Lt,nsample);

  const int nmom = pion_momenta.size();
  
  std::complex<double> dummy; double ddummy;

  for(int psnk=0;psnk<nmom;psnk++){
    if(!proj_snk(dummy,pion_momenta[psnk])) continue;

    for(int psrc=0;psrc<nmom;psrc++){
      if(!proj_src(dummy,pion_momenta[psrc])) continue;
      if(!allow(ddummy,pion_momenta[psrc],pion_momenta[psnk])) continue;

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
void readBubble(bubbleDataAllMomenta &raw_data, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, const FilenamePolicy &fn, const std::vector<threeMomentum> &pion_momenta,
		const PiPiProject &proj, const SourceOrSink src_snk){    
  const int nmom = pion_momenta.size();
  const int nsample = (traj_lessthan - traj_start)/traj_inc;

  std::complex<double> dummy;
  for(int p=0;p<nmom;p++){
    if(!proj(dummy,pion_momenta[p])) continue;

    bubbleData &into = raw_data(src_snk,pion_momenta[p]);
    
#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      int traj = traj_start + sample * traj_inc;
      std::string filename = fn(data_dir, traj, pion_momenta[p], tsep_pipi);
      std::cout << "Parsing " << filename << std::endl;
      into.parse(filename, sample);
    }
  }
}


template<typename FilenamePolicy>
void readBubble(bubbleDataAllMomenta &raw_data, const std::string &data_dir, const int tsep_pipi, const int Lt,
		const int traj_start, const int traj_inc, const int traj_lessthan, 
		const FilenamePolicy &fn_src, const FilenamePolicy &fn_snk,
		const std::vector<threeMomentum> &pion_momenta,
		const PiPiProject &proj_src, const PiPiProject &proj_snk, const PiPiMomAllow &allow){
  std::cout << "Reading bubble\n"; boost::timer::auto_cpu_timer t("Report: Read bubble in %w s\n");
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  raw_data.setup(Lt,tsep_pipi,nsample);
  
  readBubble(raw_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, fn_src, pion_momenta, proj_src, Source);
  readBubble(raw_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, fn_snk, pion_momenta, proj_snk, Sink);
}



#endif

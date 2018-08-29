#ifndef _PIPI_TO_SIGMA_READ_DATA_H_
#define _PIPI_TO_SIGMA_READ_DATA_H_

struct PiPiToSigmaBasicReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  subStringReplace fmt;
  
  PiPiToSigmaBasicReadPolicy(const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir){

    std::string ffmt = "traj_<CONF>_pipitosigma_sigmawdagmom<MOM_QUARK_SIGMA>_pionmom<MOM_PI>_v2";
    fmt.chunkString(ffmt, { subStringSpecify("<CONF>"), subStringSpecify("<MOM_QUARK_SIGMA>"), subStringSpecify("<MOM_PI>") });
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const threeMomentum &mom_quark_sigma, const threeMomentum &mom_pi) const{
    std::ostringstream os;
    os << dir << '/';
    fmt.replace(os, { anyToStr(traj_start + sample * traj_inc), momStr(mom_quark_sigma), momStr(mom_pi) } );
    return os.str();
  }
};
  

template<typename ReadPolicy>
void readPiPiToSigma(figureData &raw_data, const int Lt, const ReadPolicy &rd){
  std::cout << "Reading pipi->sigma data\n"; boost::timer::auto_cpu_timer t("Read pipi->sigma in %w s\n");
  int nsample = rd.nsample();

  raw_data.setup(Lt,nsample);
  raw_data.zero();

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };

  
  std::vector<threeMomentum> pion_mom = { {2,2,2}, {-2,-2,-2},
					  {-2,2,2}, {2,-2,-2},
					  {2,-2,2}, {-2,2,-2},
					  {2,2,-2}, {-2,-2,2} };

  figureData tmp_raw_data(Lt,nsample);

  for(int ppiidx = 0 ; ppiidx < 8 ; ppiidx++){
    for(int psigqidx = 0 ; psigqidx < 8 ; psigqidx++){
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	std::string filename = rd.filename(sample, quark_mom[psigqidx], pion_mom[ppiidx]);
	std::cout << "Parsing " << filename << std::endl;
	tmp_raw_data.parseCDR(filename, sample);
      }

      raw_data = raw_data + tmp_raw_data;
    }
  }
  
  raw_data = raw_data/64.;
}

void readPiPiToSigma(figureData &raw_data, const std::string &data_dir, const int Lt,
		     const int traj_start, const int traj_inc, const int traj_lessthan){ 
  PiPiToSigmaBasicReadPolicy rd(data_dir, traj_start, traj_inc, traj_lessthan);
  readPiPiToSigma(raw_data, Lt, rd);
}


bubbleData A1projectPiPiBubble(const bubbleDataAllMomenta &pipi_self_data, const std::vector<threeMomentum> &pion_mom, const int Lt, const int tsep_pipi){
  int nsample_raw = pipi_self_data.getNsample();
  bubbleData out(Source,Lt,tsep_pipi,nsample_raw);
  out.zero();

  for(int t=0;t<Lt;t++)
    for(int p=0;p<pion_mom.size();p++)
      out(t) = out(t) + pipi_self_data(Source,pion_mom[p])(t)/double(pion_mom.size()); //A1 project
  
  return out;
}

template<typename ReadPolicy>
bubbleData getA1projectedSourcePiPiBubble(const int Lt, const int tsep_pipi, const ReadPolicy &rp){
  bubbleDataAllMomenta pipi_self_data(Lt, tsep_pipi, rp.nsample());
  std::vector<threeMomentum> pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1} };
  readBubble(pipi_self_data, Lt, rp, pion_mom, Source);  
  return A1projectPiPiBubble(pipi_self_data, pion_mom, Lt, tsep_pipi);
}

bubbleData getA1projectedSourcePiPiBubble(const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi, const int Lt){
  readBubbleStationaryPolicy fp(false,Source);
  PiPiBubbleBasicReadPolicy<readBubbleStationaryPolicy> rp(fp, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  return getA1projectedSourcePiPiBubble(Lt, tsep_pipi, rp);
}


#endif

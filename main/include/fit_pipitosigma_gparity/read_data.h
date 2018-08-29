#ifndef _PIPI_TO_SIGMA_READ_DATA_H_
#define _PIPI_TO_SIGMA_READ_DATA_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

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


void reconstructPiPiToSigmaDisconnected(figureData &disconn, const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z){
  const int Lt = pipi_self_data_Z.getLt();
  assert(sigma_self_data_Z.getLt() == Lt);

  const int nsample = pipi_self_data_Z.getNsample();
  assert( sigma_self_data_Z.getNsample() == nsample);

  disconn.setup(Lt, nsample);
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      //note the pipi bubble is computed online as 0.5 tr( mf(t) mf(t-tsep) ), and this is combined with the sigma bubble with a coeff sqrt(6)/2 in the parallel code to form the pipi->sigma disconnected part
      //However the correct formula for the pipi bubble is -0.5 tr( mf(t) mf(t-tsep) )   ; this is corrected for when the bubble is loaded in this analysis code. Hence the coeff here needs to be -sqrt(6)/2
      rawDataDistribution<std::complex<double> > valz = -sqrt(6.)/2 * pipi_self_data_Z(tsrc) * sigma_self_data_Z(tsnk); 
      for(int s=0;s<valz.size();s++) disconn(tsrc,tsep).sample(s) = valz.sample(s).real();
    }
  }
}

//In the current runs we sum the connected and disconnected diagrams online. However for tstep > 1 this means the vacuum diagrams are also only sampled on a subset of timeslices
//This can be rectified by removing the disconnected component and performing the source timeslice average of the connected part, then afterwards constructing the disconnected part
//and source timeslice averaging that over all Lt
void reconstructPiPiToSigmaConnected(figureData &conn, const figureData &full, const figureData &disconn, const int tstep_src_full){
  const int Lt = full.getLt();
  assert(disconn.getLt() == Lt);

  const int nsample = full.getNsample();
  assert( disconn.getNsample() == nsample);

  conn.setup(Lt, nsample);
  conn.zero();
  for(int tsrc=0; tsrc< Lt; tsrc += tstep_src_full)
    for(int tsep=0;tsep<Lt;tsep++)
      conn(tsrc, tsep) = full(tsrc, tsep) - disconn(tsrc, tsep);
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

template<typename ContainerType, typename ReadPolicy>
void getA1projectedSourcePiPiBubble(ContainerType &out, const int Lt, const int tsep_pipi, const ReadPolicy &rp){
  out.setup(Source, Lt, tsep_pipi, rp.nsample());
  out.zero();
  
  ContainerType temp(Source, Lt, tsep_pipi, rp.nsample());
  std::vector<threeMomentum> pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1} };
  for(int p=0;p<pion_mom.size();p++){
    readBubble(temp, Lt, pion_mom[p], rp);
    for(int t=0;t<Lt;t++) out(t) = out(t) + temp(t) * (1./pion_mom.size());
  }
}

template<typename ContainerType>
void getA1projectedSourcePiPiBubble(ContainerType &out, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi, const int Lt){
  readBubbleStationaryPolicy fp(false,Source);
  PiPiBubbleBasicReadPolicy<readBubbleStationaryPolicy> rp(fp, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  getA1projectedSourcePiPiBubble(out, Lt, tsep_pipi, rp);
}

CPSFIT_END_NAMESPACE

#endif

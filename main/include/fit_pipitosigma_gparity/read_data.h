#ifndef _PIPI_TO_SIGMA_READ_DATA_H_
#define _PIPI_TO_SIGMA_READ_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include<pipi_common/pipi_common.h>

CPSFIT_START_NAMESPACE

struct PiPiToSigmaGenericReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  subStringReplace fmt;
  
  PiPiToSigmaGenericReadPolicy(const std::string &file_fmt, const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir){
    fmt.chunkString(file_fmt, { subStringSpecify("<CONF>"), subStringSpecify("<MOM_QUARK_SIGMA>"), subStringSpecify("<MOM_PI>") });
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const threeMomentum &mom_quark_sigma, const threeMomentum &mom_pi) const{
    std::ostringstream os;
    os << dir << '/';
    fmt.replace(os, { anyToStr(traj_start + sample * traj_inc), momStr(mom_quark_sigma), momStr(mom_pi) } );
    return os.str();
  }
};
  
struct PiPiToSigmaBasicReadPolicy: public PiPiToSigmaGenericReadPolicy{
  PiPiToSigmaBasicReadPolicy(const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    PiPiToSigmaGenericReadPolicy("traj_<CONF>_pipitosigma_sigmawdagmom<MOM_QUARK_SIGMA>_pionmom<MOM_PI>_v2", dir, traj_start, traj_inc, traj_lessthan)
  {}
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
void readPiPiToSigma(figureData &raw_data, const std::string &file_fmt, const std::string &data_dir, const int Lt,
		     const int traj_start, const int traj_inc, const int traj_lessthan){ 
  PiPiToSigmaGenericReadPolicy rd(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readPiPiToSigma(raw_data, Lt, rd);
}

//Construct disconnected part from  Re(  pipi_buble * sigma_bubble ) as we did in the parallel calculation
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

//Compute the disconnected part from Re( pipi_bubble ) * Re (sigma_bubble)   - this may have better statistical error
void reconstructPiPiToSigmaDisconnected(figureData &disconn, const bubbleData &pipi_self_data, const sigmaSelfContraction &sigma_self_data){
  const int Lt = pipi_self_data.getLt();
  assert(sigma_self_data.getLt() == Lt);

  const int nsample = pipi_self_data.getNsample();
  assert( sigma_self_data.getNsample() == nsample);

  disconn.setup(Lt, nsample);
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      disconn(tsrc,tsep) = -sqrt(6.)/2 * pipi_self_data(tsrc) * sigma_self_data(tsnk); 
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


struct readReconstructPiPiToSigmaWithDisconnAllTsrcOptions{
  //Choose a tstep for the vacuum diagram that is not 1 for testing
  bool force_disconn_tstep_src;
  int disconn_tstep_src;
  
  //Choose whether to compute the disconnected part using Re ( pipi_bubble * sigma_bubble ) or Re ( pipi_bubble ) * Re ( sigma_bubble )
  //The former is how it was computed in the parallel code, but the latter (default) has very slightly better stat error
  bool compute_disconn_ReRe;

  readReconstructPiPiToSigmaWithDisconnAllTsrcOptions(): force_disconn_tstep_src(false), compute_disconn_ReRe(true){}
};

//Combine above, reading the pipi->sigma data and reconstructing the disconnected component for all tsrc
template<typename ReadPolicy>
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const ReadPolicy &rp,
								    const int Lt, const int tstep_src,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){

  bubbleData pipi_self_data = reIm(pipi_self_data_Z, 0); //real part
  sigmaSelfContraction sigma_self_data = reIm(sigma_self_data_Z, 0); //copy real part

  //Get pipi->sigma data
  figureData pipitosigma_data;
  readPiPiToSigma(pipitosigma_data, Lt, rp);
  
  //Reconstruct disconnected and connected part
  figureData pipitosigma_disconn_data_ReZZ;
  reconstructPiPiToSigmaDisconnected(pipitosigma_disconn_data_ReZZ, pipi_self_data_Z, sigma_self_data_Z); // Re ( pipi_bubble * sigma_bubble )
 
  figureData pipitosigma_disconn_data_ReZReZ;
  reconstructPiPiToSigmaDisconnected(pipitosigma_disconn_data_ReZReZ, pipi_self_data, sigma_self_data); // Re ( pipi_bubble ) * Re ( sigma_bubble )

  figureData pipitosigma_conn_data;
  reconstructPiPiToSigmaConnected(pipitosigma_conn_data, pipitosigma_data, pipitosigma_disconn_data_ReZZ, tstep_src);

  //(Very slightly) better statistics if we use the Re ( pipi_bubble ) * Re ( sigma_bubble ) for the disconnected part, taking advantage of the fact that the bubbles are real under the ensemble avg
  figureData &pipitosigma_disconn_data = opt.compute_disconn_ReRe ? pipitosigma_disconn_data_ReZReZ : pipitosigma_disconn_data_ReZZ;

  //The code computes the disconnected component for all tsrc, but this option can be used to constrain the number of source timeslices to observe the effect
  if(opt.force_disconn_tstep_src){
    for(int tsrc=0;tsrc<Lt;tsrc++){
      if(tsrc % opt.disconn_tstep_src != 0)
	for(int t=0;t<Lt;t++) pipitosigma_disconn_data(tsrc, t).zero();
    }
  }

  //Source avg connected and disconnected parts and sum the contributions
  rawCorrelationFunction correlator_raw_conn = sourceAverage(pipitosigma_conn_data);
  rawCorrelationFunction correlator_raw_disconn = sourceAverage(pipitosigma_disconn_data);

  rawCorrelationFunction correlator_raw = correlator_raw_conn;
  for(int t=0;t<Lt;t++) correlator_raw.value(t) = correlator_raw.value(t) + correlator_raw_disconn.value(t);

  std::cout << "Pipi->sigma raw data connected/disconnected parts:\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << correlator_raw_conn.value(t) << " " << correlator_raw_disconn.value(t) << std::endl;

  return correlator_raw;
}

rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &data_dir, const int Lt, const int tstep_src,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  PiPiToSigmaBasicReadPolicy rd(data_dir, traj_start, traj_inc, traj_lessthan);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,  pipi_self_data_Z, sigma_self_data_Z, opt);
}
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &file_fmt, const std::string &data_dir, const int Lt, const int tstep_src,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  PiPiToSigmaGenericReadPolicy rd(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,  pipi_self_data_Z, sigma_self_data_Z, opt);
}


//Averages over pion momenta to produce a rotationally invariant state
template<typename AllMomentaContainerType>
typename AllMomentaContainerType::ContainerType A1projectSourcePiPiBubble(const AllMomentaContainerType &pipi_self_data, const std::vector<threeMomentum> &pion_mom){
  int nsample_raw = pipi_self_data.getNsample();
  int Lt = pipi_self_data.getLt();
  int tsep_pipi = pipi_self_data.getTsepPiPi();
  typename AllMomentaContainerType::ContainerType out(Source,Lt,tsep_pipi,nsample_raw);
  out.zero();

  for(int t=0;t<Lt;t++)
    for(int p=0;p<pion_mom.size();p++)
      out(t) = out(t) + pipi_self_data(Source,pion_mom[p])(t)/double(pion_mom.size()); //A1 project
  
  return out;
}

rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &data_dir, const int Lt, const int tstep_src,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const std::vector<threeMomentum> &pion_mom,
								    const bubbleDataAllMomentaZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  bubbleDataZ pipi_self_A1_Z = A1projectSourcePiPiBubble(pipi_self_data_Z, pion_mom);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(data_dir, Lt, tstep_src, traj_start, traj_inc, traj_lessthan, pipi_self_A1_Z, sigma_self_data_Z, opt);
}
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &file_fmt, const std::string &data_dir, const int Lt, const int tstep_src,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const std::vector<threeMomentum> &pion_mom,
								    const bubbleDataAllMomentaZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  bubbleDataZ pipi_self_A1_Z = A1projectSourcePiPiBubble(pipi_self_data_Z, pion_mom);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(file_fmt, data_dir, Lt, tstep_src, traj_start, traj_inc, traj_lessthan, pipi_self_A1_Z, sigma_self_data_Z, opt);
}


template<typename ContainerType, typename ReadPolicy>
void getA1projectedSourcePiPiBubble(ContainerType &out, const int Lt, const int tsep_pipi, const std::vector<threeMomentum> &pion_mom, const ReadPolicy &rp){
  out.setup(Source, Lt, tsep_pipi, rp.nsample());
  out.zero();
  
  ContainerType temp(Source, Lt, tsep_pipi, rp.nsample());
  for(int p=0;p<pion_mom.size();p++){
    readBubble(temp, Lt, pion_mom[p], rp);
    for(int t=0;t<Lt;t++) out(t) = out(t) + temp(t) * (1./pion_mom.size());
  }
}


template<typename ContainerType, typename ReadPolicy>
inline void getA1projectedSourcePiPiBubble(ContainerType &out, const int Lt, const int tsep_pipi, const ReadPolicy &rp){
  std::vector<threeMomentum> pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1} };
  getA1projectedSourcePiPiBubble(out, Lt, tsep_pipi, pion_mom, rp);
}


template<typename ContainerType>
void getA1projectedSourcePiPiBubble(ContainerType &out, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int tsep_pipi, const int Lt){
  readBubbleStationaryPolicy fp(false,Source);
  PiPiBubbleBasicReadPolicy<readBubbleStationaryPolicy> rp(fp, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  getA1projectedSourcePiPiBubble(out, Lt, tsep_pipi, rp);
}

CPSFIT_END_NAMESPACE

#endif

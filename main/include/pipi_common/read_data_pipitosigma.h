#ifndef _PIPI_TO_SIGMA_READ_DATA_H_
#define _PIPI_TO_SIGMA_READ_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include "raw_correlator.h"
#include "read_data_pipi.h"
#include "mom_project.h"

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
void readPiPiToSigma(figureData &raw_data, const int Lt, const PiPiProjectorBase &proj_pipi, const ReadPolicy &rd){
  std::cout << "Reading pipi->sigma data\n"; boost::timer::auto_cpu_timer t("Read pipi->sigma in %w s\n");
  int nsample = rd.nsample();

  raw_data.setup(Lt,rawDataDistributionD(nsample));
  raw_data.zero();

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };

  
  figureData tmp_raw_data(Lt,rawDataDistributionD(nsample));

  for(int ppiidx = 0 ; ppiidx < proj_pipi.nMomenta() ; ppiidx++){
    threeMomentum ppi = proj_pipi.momentum(ppiidx) * 2; //Tianle's conventions for the pion energy are in units of pi/2L     
    double pipi_coeff = std::real(proj_pipi.coefficient(ppiidx));    

    for(int psigqidx = 0 ; psigqidx < 8 ; psigqidx++){
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	std::string filename = rd.filename(sample, quark_mom[psigqidx], ppi);
	std::cout << "Parsing " << filename << std::endl;
	tmp_raw_data.parseCDR(filename, sample);
      }

      raw_data = raw_data + 1./8 * pipi_coeff*tmp_raw_data;
    }
  }
}

void readPiPiToSigma(figureData &raw_data, const std::string &data_dir, const int Lt, const PiPiProjectorBase &proj_pipi,
		     const int traj_start, const int traj_inc, const int traj_lessthan){ 
  PiPiToSigmaBasicReadPolicy rd(data_dir, traj_start, traj_inc, traj_lessthan);
  readPiPiToSigma(raw_data, Lt, proj_pipi, rd);
}
void readPiPiToSigma(figureData &raw_data, const std::string &file_fmt, const std::string &data_dir, const int Lt, const PiPiProjectorBase &proj_pipi,
		     const int traj_start, const int traj_inc, const int traj_lessthan){ 
  PiPiToSigmaGenericReadPolicy rd(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readPiPiToSigma(raw_data, Lt, proj_pipi, rd);
}

//Construct disconnected part from  Re(  pipi_buble * sigma_bubble ) as we did in the parallel calculation
//pipi_self_data_Z should be pre-projected
void reconstructPiPiToSigmaDisconnected(figureData &disconn, const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z){
  const int Lt = pipi_self_data_Z.getLt();
  assert(sigma_self_data_Z.getLt() == Lt);

  const int nsample = pipi_self_data_Z(0).size();
  assert( sigma_self_data_Z(0).size() == nsample);

  disconn.setup(Lt, rawDataDistributionD(nsample));
  for(int tsrc=0;tsrc<Lt;tsrc++){
    for(int tsep=0;tsep<Lt;tsep++){
      int tsnk = ( tsrc + tsep ) % Lt;
      //note the pipi bubble is computed online as 0.5 tr( mf(t) mf(t-tsep) ), and this is combined with the sigma bubble with a coeff sqrt(6)/2 in the parallel code to form the pipi->sigma disconnected part
      //However the correct formula for the pipi bubble is -0.5 tr( mf(t) mf(t-tsep) )   ; this is corrected for when the bubble is loaded in this analysis code. Hence the coeff here needs to be -sqrt(6)/2
      rawDataDistribution<std::complex<double> > valz = -sqrt(6.)/2 * pipi_self_data_Z(tsrc) * sigma_self_data_Z(tsnk); 
      for(int s=0;s<nsample;s++) disconn(tsrc,tsep).sample(s) = valz.sample(s).real();
    }
  }
}

//Compute the disconnected part from Re( pipi_bubble ) * Re (sigma_bubble)   - this has a slightly better statistical error than the above
//pipi_self_data should be pre-projected
void reconstructPiPiToSigmaDisconnected(figureData &disconn, const bubbleData &pipi_self_data, const sigmaSelfContraction &sigma_self_data){
  const int Lt = pipi_self_data.getLt();
  assert(sigma_self_data.getLt() == Lt);

  const int nsample = pipi_self_data(0).size();
  assert( sigma_self_data(0).size() == nsample);

  disconn.setup(Lt);
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

  const int nsample = full(0,0).size();
  assert( disconn(0,0).size() == nsample);

  conn.setup(Lt, rawDataDistributionD(nsample,0.));
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
//pipi_self_data_Z should be pre-projected
template<typename ReadPolicy>
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const ReadPolicy &rp,
								    const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){

  bubbleData pipi_self_data = reIm(pipi_self_data_Z, 0); //real part
  sigmaSelfContraction sigma_self_data = reIm(sigma_self_data_Z, 0); //copy real part

  //Get pipi->sigma data
  figureData pipitosigma_data;
  readPiPiToSigma(pipitosigma_data, Lt, proj_pipi, rp);
  
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

rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  PiPiToSigmaBasicReadPolicy rd(data_dir, traj_start, traj_inc, traj_lessthan);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,  proj_pipi, pipi_self_data_Z, sigma_self_data_Z, opt);
}
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &file_fmt, const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,
								    const bubbleDataZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  PiPiToSigmaGenericReadPolicy rd(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(rd, Lt, tstep_src,  proj_pipi, pipi_self_data_Z, sigma_self_data_Z, opt);
}


//Averages over pion momenta to produce a rotationally invariant state
template<typename RealOrComplex>
struct _getReOrZ{
  static inline double get(const std::complex<double> &val){ return std::real(val); }
};
template<typename T>
struct _getReOrZ<std::complex<T> >{
  static inline std::complex<double> get(const std::complex<double> &val){ return val; }
};
template<typename BubbleDataType>
struct projectSourcePiPiBubble_getCoeff{
  typedef typename BubbleDataType::DistributionType::DataType DataType;
  static inline auto getCoeff(const std::complex<double> &val)->decltype(_getReOrZ<DataType>::get(val)){ return _getReOrZ<DataType>::get(val); }
};

template<typename AllMomentaContainerType>
typename AllMomentaContainerType::ContainerType projectSourcePiPiBubble(const AllMomentaContainerType &pipi_self_data, const PiPiProjectorBase &proj_pipi){
  int Lt = pipi_self_data.getLt();
  int tsep_pipi = pipi_self_data.getTsepPiPi();
  typedef typename AllMomentaContainerType::ContainerType OutType;
  OutType out(Source,Lt,tsep_pipi);

  for(int t=0;t<Lt;t++)
    for(int p=0;p<proj_pipi.nMomenta();p++){
      typename AllMomentaContainerType::DistributionType tmp = projectSourcePiPiBubble_getCoeff<OutType>::getCoeff(proj_pipi.coefficient(p)) * pipi_self_data(Source,proj_pipi.momentum(p))(t);
      out(t) = p == 0 ? tmp : out(t) + tmp;
    }
  
  return out;
}

//This version does the pipi projection
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,								    
								    const bubbleDataAllMomentaZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  bubbleDataZ pipi_self_proj_Z = projectSourcePiPiBubble(pipi_self_data_Z, proj_pipi);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(data_dir, Lt, tstep_src, proj_pipi, traj_start, traj_inc, traj_lessthan, pipi_self_proj_Z, sigma_self_data_Z, opt);
}
rawCorrelationFunction readReconstructPiPiToSigmaWithDisconnAllTsrc(const std::string &file_fmt, const std::string &data_dir, const int Lt, const int tstep_src, const PiPiProjectorBase &proj_pipi,
								    const int traj_start, const int traj_inc, const int traj_lessthan,								    
								    const bubbleDataAllMomentaZ &pipi_self_data_Z, const sigmaSelfContractionZ &sigma_self_data_Z,
								    const readReconstructPiPiToSigmaWithDisconnAllTsrcOptions &opt = readReconstructPiPiToSigmaWithDisconnAllTsrcOptions()){
  bubbleDataZ pipi_self_proj_Z = projectSourcePiPiBubble(pipi_self_data_Z, proj_pipi);
  return readReconstructPiPiToSigmaWithDisconnAllTsrc(file_fmt, data_dir, Lt, tstep_src, proj_pipi, traj_start, traj_inc, traj_lessthan, pipi_self_proj_Z, sigma_self_data_Z, opt);
}

//Read and project pipi bubble
template<typename ContainerType, typename ReadPolicy>
void getProjectedSourcePiPiBubble(ContainerType &out, const int Lt, const int tsep_pipi, const PiPiProjectorBase &proj_pipi, const ReadPolicy &rp){
  out.setup(Source, Lt, tsep_pipi);
  
  ContainerType temp(Source, Lt, tsep_pipi);
  for(int p=0;p<proj_pipi.nMomenta();p++){
    readPiPiBubble(temp, Lt, proj_pipi.momentum(p), rp);
    for(int t=0;t<Lt;t++){
      typename ContainerType::DistributionType tmp =  projectSourcePiPiBubble_getCoeff<ContainerType>::getCoeff(proj_pipi.coefficient(p)) * temp(t);
      out(t) = p==0 ? tmp : out(t) + tmp;
    }
  }
}

template<typename ContainerType>
void getProjectedSourcePiPiBubble(ContainerType &out, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, 
				  const int Lt, const int tsep_pipi,  const PiPiProjectorBase &proj_pipi){
  readBubbleStationaryPolicy fp(false,Source);
  PiPiBubbleBasicReadPolicy<readBubbleStationaryPolicy> rp(fp, data_dir, traj_start, traj_inc, traj_lessthan, tsep_pipi);
  getProjectedSourcePiPiBubble(out, Lt, tsep_pipi, proj_pipi, rp);
}

CPSFIT_END_NAMESPACE

#endif

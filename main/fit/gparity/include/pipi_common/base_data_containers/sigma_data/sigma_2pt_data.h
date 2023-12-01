#pragma once

#include "../figure_data_container.h"
#include "../../correlator_utils/threemomentum.h"
#include "sigma_reconstruct_disconn.h"

CPSFIT_START_NAMESPACE

////////////////////
//FILENAME POLICIES
//////////////////
//General read policy producing the filename associated with a particular sample index and momentum combination
struct Sigma2ptGenericReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  subStringReplace fmt;
  
  Sigma2ptGenericReadPolicy(const std::string &fmt_str, const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir){
    fmt.chunkString(fmt_str, { subStringSpecify("<CONF>"), subStringSpecify("<PSRC_QUARK>"), subStringSpecify("<PSNK_QUARK>") });
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const threeMomentum &psrc_quark, const threeMomentum &psnk_quark) const{
    std::ostringstream os;
    os << dir << '/';
    fmt.replace(os, { anyToStr(traj_start + sample * traj_inc), momStr(psrc_quark), momStr(psnk_quark) } );
    return os.str();
  }
};
//Implementation of the above with a particular filename format
struct Sigma2ptBasicReadPolicy: public Sigma2ptGenericReadPolicy{
  Sigma2ptBasicReadPolicy(const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    Sigma2ptGenericReadPolicy("traj_<CONF>_sigmacorr_mompsrc<PSRC_QUARK>psnk<PSNK_QUARK>_v2", dir, traj_start, traj_inc, traj_lessthan)
  {}
};

///////////////////////
//READ SIGMA 2PT DATA
//Note, both the connected and disconnected components are combined on the measurement code
//////////////////////

struct readSigmaSigmaOptions{
  bool include_V_diagram;
  sigmaSelfContractionZ const *sigma_self_data_Z; //required for include_V_diagram

  readSigmaSigmaOptions(): include_V_diagram(true), sigma_self_data_Z(nullptr){}
};

template<typename ReadPolicy>
void readSigmaSigma(figureData &raw_data, const int Lt, const ReadPolicy &rp, const readSigmaSigmaOptions &opts = readSigmaSigmaOptions()){
  std::cout << "Reading sigma 2pt data\n"; boost::timer::auto_cpu_timer t("Read sigma 2pt in %w s\n");
  int nsample = rp.nsample();

  raw_data.setup(Lt,rawDataDistributionD(nsample));
  raw_data.zero();

  static std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
						  {-3,1,1}, {3,-1,-1},
						  {1,-3,1}, {-1,3,-1},
						  {1,1,-3}, {-1,-1,3} };
  figureData tmp_raw_data(Lt,rawDataDistributionD(nsample));

  const int nmom = quark_mom.size();
  
  for(int psnk=0;psnk<nmom;psnk++){
    for(int psrc=0;psrc<nmom;psrc++){
            
#pragma omp parallel for
      for(int sample=0; sample < nsample; sample++){
	std::string filename = rp.filename(sample, quark_mom[psnk], quark_mom[psrc]);
	std::cout << "Parsing " << filename << std::endl;
	tmp_raw_data.parseCDR(filename, sample);
      }

      raw_data = raw_data + tmp_raw_data;
    }
  }
  
  raw_data = raw_data/double(nmom*nmom);

  if(!opts.include_V_diagram){
    assert(opts.sigma_self_data_Z != nullptr);
    figureData disconn;
    reconstructSigma2ptDisconnected(disconn, *opts.sigma_self_data_Z);
    figureData tmp(raw_data);
    reconstructSigma2ptConnected(raw_data,tmp,disconn);
  }
}
//Call the above with the default filename format
inline void readSigmaSigma(figureData &raw_data, const std::string &data_dir, const int Lt,
		    const int traj_start, const int traj_inc, const int traj_lessthan, const readSigmaSigmaOptions &opts = readSigmaSigmaOptions()){
  Sigma2ptBasicReadPolicy rp(data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSigma(raw_data, Lt, rp, opts);
} 
//Call with a user-specified file format
inline void readSigmaSigma(figureData &raw_data, const std::string &file_fmt, const std::string &data_dir, const int Lt,
		    const int traj_start, const int traj_inc, const int traj_lessthan, const readSigmaSigmaOptions &opts = readSigmaSigmaOptions()){
  Sigma2ptGenericReadPolicy rp(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSigma(raw_data, Lt, rp, opts);
} 


CPSFIT_END_NAMESPACE

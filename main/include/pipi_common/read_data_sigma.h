#ifndef _SIGMA_READ_DATA_H_
#define _SIGMA_READ_DATA_H_

#include <boost/timer/timer.hpp>

#include<config.h>
#include<utils/macros.h>

#include "threemomentum.h"
#include "data_containers.h"

CPSFIT_START_NAMESPACE

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

struct Sigma2ptBasicReadPolicy: public Sigma2ptGenericReadPolicy{
  Sigma2ptBasicReadPolicy(const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    Sigma2ptGenericReadPolicy("traj_<CONF>_sigmacorr_mompsrc<PSRC_QUARK>psnk<PSNK_QUARK>_v2", dir, traj_start, traj_inc, traj_lessthan)
  {}
};



template<typename ReadPolicy>
void readSigmaSigma(figureData &raw_data, const int Lt, const ReadPolicy &rp){
  std::cout << "Reading sigma 2pt data\n"; boost::timer::auto_cpu_timer t("Read sigma 2pt in %w s\n");
  int nsample = rp.nsample();

  raw_data.setup(Lt,nsample);
  raw_data.zero();

  static std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
						  {-3,1,1}, {3,-1,-1},
						  {1,-3,1}, {-1,3,-1},
						  {1,1,-3}, {-1,-1,3} };
  figureData tmp_raw_data(Lt,nsample);

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
}

inline void readSigmaSigma(figureData &raw_data, const std::string &data_dir, const int Lt,
		    const int traj_start, const int traj_inc, const int traj_lessthan){
  Sigma2ptBasicReadPolicy rp(data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSigma(raw_data, Lt, rp);
} 
inline void readSigmaSigma(figureData &raw_data, const std::string &file_fmt, const std::string &data_dir, const int Lt,
		    const int traj_start, const int traj_inc, const int traj_lessthan){
  Sigma2ptGenericReadPolicy rp(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSigma(raw_data, Lt, rp);
} 

struct SigmaSelfGenericReadPolicy{
  int traj_start, traj_inc, traj_lessthan;
  std::string dir;
  subStringReplace fmt;
  
  SigmaSelfGenericReadPolicy(const std::string &file_fmt, const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), dir(dir){
    fmt.chunkString(file_fmt, { subStringSpecify("<CONF>"), subStringSpecify("<PQUARK>") });
  }

  int nsample() const{ return (traj_lessthan - traj_start)/traj_inc; }
  
  std::string filename(const int sample, const threeMomentum &pquark) const{
    std::ostringstream os;
    os << dir << '/';
    fmt.replace(os, { anyToStr(traj_start + sample * traj_inc), momStr(pquark)} );
    return os.str();
  }
};

struct SigmaSelfBasicReadPolicy: public SigmaSelfGenericReadPolicy{
  SigmaSelfBasicReadPolicy(const std::string &dir, const int traj_start, const int traj_inc, const int traj_lessthan):
    SigmaSelfGenericReadPolicy("traj_<CONF>_sigmaself_mom<PQUARK>_v2", dir, traj_start, traj_inc, traj_lessthan){}
};



template<typename ContainerType, typename ReadPolicy>
void readSigmaSelf(ContainerType &raw_data, const int Lt, const ReadPolicy &rp){
  const int nsample = rp.nsample();

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };
  const int nmom = quark_mom.size();

  raw_data.setup(Lt,nsample);
  raw_data.zero();

  ContainerType tmp_data(Lt,nsample);

  for(int p=0;p<nmom;p++){    
#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      std::string filename = rp.filename(sample, quark_mom[p]);
      std::cout << "Parsing " << filename << std::endl;
      tmp_data.parse(filename, sample);
    }
    
    for(int t=0;t<Lt;t++) raw_data(t) = raw_data(t) + tmp_data(t) / double(nmom);
  }
}

template<typename ContainerType>
inline void readSigmaSelf(ContainerType &raw_data, const std::string &data_dir, const int Lt,
			  const int traj_start, const int traj_inc, const int traj_lessthan){   
  SigmaSelfBasicReadPolicy rp(data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSelf(raw_data, Lt, rp);
} 
template<typename ContainerType>
inline void readSigmaSelf(ContainerType &raw_data, const std::string &file_fmt, 
			  const std::string &data_dir, const int Lt,
			  const int traj_start, const int traj_inc, const int traj_lessthan){   
  SigmaSelfGenericReadPolicy rp(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSelf(raw_data, Lt, rp);
} 

CPSFIT_END_NAMESPACE

#endif

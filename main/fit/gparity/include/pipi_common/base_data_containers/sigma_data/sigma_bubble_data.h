#pragma once

#include "../sigma_bubble_data_container.h"
#include "../../correlator_utils/threemomentum.h"

CPSFIT_START_NAMESPACE

////////////////////
//FILENAME POLICIES
//////////////////

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

///////////////////////
//READ SIGMA BUBBLE DATA
//////////////////////

//ContainerType should be one of the bubble data containers
template<typename ContainerType, typename ReadPolicy>
void readSigmaSelf(ContainerType &raw_data, const int Lt, const ReadPolicy &rp){
  const int nsample = rp.nsample();

  typedef typename ContainerType::DistributionType DistributionType;

  std::vector<threeMomentum> quark_mom = { {1,1,1}, {-1,-1,-1},
					   {-3,1,1}, {3,-1,-1},
					   {1,-3,1}, {-1,3,-1},
					   {1,1,-3}, {-1,-1,3} };
  const int nmom = quark_mom.size();

  DistributionType zero(nsample); zeroit(zero);

  raw_data.setup(Lt,zero);

  ContainerType tmp_data(Lt,zero);

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
//Call the above with the default filename format
template<typename ContainerType>
inline void readSigmaSelf(ContainerType &raw_data, const std::string &data_dir, const int Lt,
			  const int traj_start, const int traj_inc, const int traj_lessthan){   
  SigmaSelfBasicReadPolicy rp(data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSelf(raw_data, Lt, rp);
} 
//Call the above with a user-specified filename format
template<typename ContainerType>
inline void readSigmaSelf(ContainerType &raw_data, const std::string &file_fmt, 
			  const std::string &data_dir, const int Lt,
			  const int traj_start, const int traj_inc, const int traj_lessthan){   
  SigmaSelfGenericReadPolicy rp(file_fmt, data_dir, traj_start, traj_inc, traj_lessthan);
  readSigmaSelf(raw_data, Lt, rp);
} 



CPSFIT_END_NAMESPACE

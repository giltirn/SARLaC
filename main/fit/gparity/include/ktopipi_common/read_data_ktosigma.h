#ifndef _FIT_KTOSIGMA_READ_DATA_H
#define _FIT_KTOSIGMA_READ_DATA_H

#include<config.h>
#include<utils/macros.h>

#include<algorithm>
#include<parser.h>

#include<pipi_common/correlator_utils/threemomentum.h>
#include<pipi_common/base_data_containers/sigma_bubble_data_container.h>

#include "data_containers.h"

SARLAC_START_NAMESPACE

struct KtoSigmaFilenamePolicyGen{
  subStringReplace repl_type2; //File format expects <TRAJ> <TSEP_K_SIGMA>
  subStringReplace repl_type3; //Same as type 2
  subStringReplace repl_type4; //File format expects <TRAJ>
  
  KtoSigmaFilenamePolicyGen(const std::string &type2_fmt, 
			    const std::string &type3_fmt,
			    const std::string &type4_fmt):
    repl_type2(type2_fmt, { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_SIGMA>") }),
    repl_type3(type3_fmt, { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_SIGMA>") }),
    repl_type4(type4_fmt, { subStringSpecify("<TRAJ>") })
  {}
  
  inline std::string type2filename(const std::string &data_dir, const int traj, const int tsep_k_sigma) const{
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type2.replace(fname, { anyToStr(traj), anyToStr(tsep_k_sigma) });
    return fname.str();
  }
  inline std::string type3filename(const std::string &data_dir, const int traj, const int tsep_k_sigma) const{
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type3.replace(fname, { anyToStr(traj), anyToStr(tsep_k_sigma) });
    return fname.str();
  }
  inline std::string type4filename(const std::string &data_dir, const int traj) const{
    std::ostringstream fname;
    fname << data_dir << '/';
    repl_type4.replace(fname, { anyToStr(traj) });
    return fname.str();
  }
};


template<typename FilenamePolicy>
class BasicKtoSigmaReadPolicy{
  const std::string data_dir;
  const int traj_start;
  const int traj_inc;
  const int traj_lessthan;
  const FilenamePolicy &fp;

  int traj(const int sample) const{ return traj_start + sample*traj_inc; }

public:
  int nSample() const{ return (traj_lessthan - traj_start)/traj_inc; }

  inline std::string type2filename(const int sample, const int tsep_k_sigma) const{
    return fp.type2filename(data_dir, traj(sample), tsep_k_sigma);
  }
  inline std::string type3filename(const int sample, const int tsep_k_sigma) const{
    return fp.type3filename(data_dir, traj(sample), tsep_k_sigma);
  }
  inline std::string type4filename(const int sample) const{
    return fp.type4filename(data_dir, traj(sample));
  }

  BasicKtoSigmaReadPolicy(const std::string &data_dir, int traj_start, int traj_inc, int traj_lessthan, const FilenamePolicy &fp): 
    data_dir(data_dir), traj_start(traj_start), traj_inc(traj_inc), traj_lessthan(traj_lessthan), fp(fp){}
};





//type can be 2,3,4
template<typename ReadPolicy>
type1234Data readKtoSigmaType(const int type, const int tsep_k_sigma, const int Lt, const ReadPolicy &rp){
  int nsample = rp.nSample();
  if(type == 2 || type == 3){
    type1234Data typedata(type,Lt,nsample,Sigma);
#pragma omp parallel for
    for(int sample=0;sample<nsample;sample++){
      std::string filename = type == 2 ? rp.type2filename(sample, tsep_k_sigma) : rp.type3filename(sample, tsep_k_sigma);
      typedata.parse(filename,sample);
    }
    return typedata;
  }else if(type == 4){
    type1234Data typedata(type,Lt,nsample,Sigma);
#pragma omp parallel for
    for(int sample=0;sample<nsample;sample++){
      std::string filename = rp.type4filename(sample);
      typedata.parse(filename, sample);
    }
    return typedata;
  }else{
    error_exit(std::cout << "readKtoSigmaType invalid type " << type << std::endl);
  }
}

//For original and extended data, use bubble_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
//                                               	       { {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
//                      				       { {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
//                      				       { {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       }
NumericTensor<rawDataDistributionD,1> readProjectedSigmaBubble(const std::string &data_dir, const std::string &file_fmt,
							       const int traj_start, const int traj_inc, const int traj_lessthan, 
							       const int Lt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  NumericTensor<rawDataDistributionD,1> out({Lt}, rawDataDistributionD(nsample,0.));
  
  //File format expects <TRAJ>  <MOM>
  static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>"),  subStringSpecify("<MOM>") };
  subStringReplace repl(file_fmt, keys);
  
  for(int p=0;p<bubble_quarkmom_proj.size();p++){
    sigmaSelfContraction self(Lt,rawDataDistributionD(nsample));

#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      int traj = traj_start + sample * traj_inc;
      std::ostringstream filename; filename << data_dir << '/';
      repl.replace(filename, { anyToStr(traj), momStr(bubble_quarkmom_proj[p].first) });

      std::cout << "Parsing " << filename.str() << std::endl;
      self.parse(filename.str(), sample);
    }

    for(int t=0;t<Lt;t++)
      out({t}) = out({t}) + self(t) * bubble_quarkmom_proj[p].second; //project
  }

  return out;
}


SARLAC_END_NAMESPACE

#endif

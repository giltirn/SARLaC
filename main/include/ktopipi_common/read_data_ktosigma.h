#ifndef _FIT_KTOSIGMA_READ_DATA_H
#define _FIT_KTOSIGMA_READ_DATA_H

#include<config.h>
#include<utils/macros.h>

#include<algorithm>
#include<parser.h>

#include<pipi_common/read_data.h>

#include "data_containers.h"

CPSFIT_START_NAMESPACE

//type can be 2,3,4
type1234Data readKtoSigmaType(const int type, const int traj_start, const int traj_inc, const int traj_lessthan,
		      const int tsep_k_sigma, const int Lt, 
		      const std::string &data_dir, const std::string &file_fmt){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  if(type == 2 || type == 3){
    //File format expected <TRAJ> <TSEP_K_SIGMA>
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_SIGMA>") };
    subStringReplace repl(file_fmt, keys);

    type1234Data typedata(type,Lt,nsample,Sigma);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;

      std::ostringstream fname;
      fname << data_dir << '/';
      repl.replace(fname, { anyToStr(traj), anyToStr(tsep_k_sigma) });
      typedata.parse(fname.str(),i);
    }
    return typedata;
  }else if(type == 4){
    //File format expects <TRAJ>
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>") };
    subStringReplace repl(file_fmt, keys);

    type1234Data typedata(type,Lt,nsample,Sigma);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;

      std::ostringstream fname;
      fname << data_dir << '/';
      repl.replace(fname, { anyToStr(traj) });
      typedata.parse(fname.str(),i);
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
    sigmaSelfContraction self(Lt,nsample);

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


CPSFIT_END_NAMESPACE

#endif

#pragma once

#include "../../correlator_utils/raw_correlator.h"
#include "sigma_bubble_data.h"
#include "sigma_2pt_data.h"

SARLAC_START_NAMESPACE

struct readSigma2ptOpts{
  bool include_V_diagram;
  readSigma2ptOpts(): include_V_diagram(true){}
};


inline void readSigma2pt(rawDataCorrelationFunctionD &sigma2pt_raw, sigmaSelfContractionZ &raw_bubble_data,
			 const std::string &data_dir, 
			 const std::string &sigma2pt_file_fmt, const std::string &bubble_file_fmt, 
			 const int Lt,
			 const int traj_start, const int traj_inc, const int traj_lessthan, const readSigma2ptOpts &opts = readSigma2ptOpts()){
  readSigmaSelf(raw_bubble_data, bubble_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
  readSigmaSigmaOptions opt_ss;
  if(!opts.include_V_diagram){
    opt_ss.include_V_diagram = false;
    opt_ss.sigma_self_data_Z = &raw_bubble_data;
  }
  figureData sigma2pt_data;
  readSigmaSigma(sigma2pt_data, sigma2pt_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan, opt_ss);
  sigma2pt_raw = sourceAverage(sigma2pt_data);
}

SARLAC_END_NAMESPACE

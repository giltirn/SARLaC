#pragma once

#include "../../correlator_utils/raw_correlator.h"
#include "sigma_bubble_data.h"
#include "sigma_2pt_data.h"

CPSFIT_START_NAMESPACE

inline void readSigma2pt(rawDataCorrelationFunctionD &sigma2pt_raw, sigmaSelfContractionZ &raw_bubble_data,
			 const std::string &data_dir, 
			 const std::string &sigma2pt_file_fmt, const std::string &bubble_file_fmt, 
			 const int Lt,
			 const int traj_start, const int traj_inc, const int traj_lessthan){
  figureData sigma2pt_data;
  readSigmaSigma(sigma2pt_data, sigma2pt_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
  sigma2pt_raw = sourceAverage(sigma2pt_data);
  readSigmaSelf(raw_bubble_data, bubble_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
}

CPSFIT_END_NAMESPACE

#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_LOAD_FRZPARAMS_H
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_LOAD_FRZPARAMS_H

#include<config.h>
#include<utils/macros.h>

#include "args.h"

SARLAC_START_NAMESPACE

void loadFrozenParameters(jackknifeDistributionD &mK, jackknifeDistributionD &cK,
			  jackknifeDistributionD &E0, jackknifeDistributionD &E1,
			  NumericSquareMatrix<jackknifeDistributionD> &coeffs,
			  const Args &args){
  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc/args.bin_size;	  

  {
    std::vector<jackknifeDistributionD> p;
    readParamsStandard(p,  args.input_params.kaon2pt_fit_result);
    mK = p[args.input_params.idx_mK];
    cK = sqrt( p[args.input_params.idx_cK] );
  }

  {
    double scale = sqrt(args.input_params.pipi_sigma_sim_fit_Ascale);
    std::vector<jackknifeDistributionD> p;
    readParamsStandard(p, args.input_params.pipi_sigma_sim_fit_result);
    for(int i=0;i<p.size();i++) assert(p[i].size() == nsample);
    coeffs(0,0) = p[args.input_params.idx_coeff_pipi_state0] * scale;
    coeffs(0,1) = p[args.input_params.idx_coeff_pipi_state1] * scale;
    coeffs(1,0) = p[args.input_params.idx_coeff_sigma_state0] * scale;
    coeffs(1,1) = p[args.input_params.idx_coeff_sigma_state1] * scale;
    E0 = p[args.input_params.idx_E0];
    E1 = p[args.input_params.idx_E1];
  }

  std::cout << "cK = " << cK << std::endl;
  std::cout << "mK = " << mK << std::endl;
  std::cout << "E0 = " << E0 << std::endl;
  std::cout << "E1 = " << E1 << std::endl;
}

SARLAC_END_NAMESPACE

#endif

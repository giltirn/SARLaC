#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_INPUT_PARAM_ARGS_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_INPUT_PARAM_ARGS_H_

#include<config.h>
#include<utils/macros.h>
#include<parser.h>

SARLAC_START_NAMESPACE


#define INPUT_PARAM_ARGS_MEMBERS			\
    ( std::string, kaon2pt_fit_result )						\
    ( int, idx_cK )							\
    ( int, idx_mK )							\
    ( std::string, pipi_sigma_sim_fit_result )				\
    ( int, idx_E0 )							\
    ( int, idx_E1 )							\
    ( int, idx_E2 )							\
    ( int, idx_coeff_pipi_state0 )					\
    ( int, idx_coeff_pipi_state1 )					\
    ( int, idx_coeff_pipi_state2 )					\
    ( int, idx_coeff_pipi_exc_state0 )					\
    ( int, idx_coeff_pipi_exc_state1 )					\
    ( int, idx_coeff_pipi_exc_state2 )					\
    ( int, idx_coeff_sigma_state0 )					\
    ( int, idx_coeff_sigma_state1 )					\
    ( int, idx_coeff_sigma_state2 )					\
    ( double, pipi_sigma_sim_fit_Ascale )

//State-2 indices are ignored for 2-state fits
struct InputParamArgs{
  GENERATE_MEMBERS(INPUT_PARAM_ARGS_MEMBERS);

  InputParamArgs(): kaon2pt_fit_result("/home/ckelly/projects/32nt64_MDWF_DSDR_GparityXYZ_fixedRNG_fullanalysis/87cfgs_DD2_stream_ext/mK/params.hdf5"),
		    idx_cK(0), idx_mK(1),
		    pipi_sigma_sim_fit_result("/home/ckelly/projects/32nt64_MDWF_DSDR_GparityXYZ_fixedRNG_fullanalysis/87cfgs_DD2_stream_ext/Epipi/I0/sim_gnd_sigma/params.hdf5"),
		    idx_coeff_pipi_state0(0), idx_coeff_pipi_state1(1), 
		    idx_coeff_pipi_exc_state0(2), idx_coeff_pipi_exc_state1(3), 
		    idx_coeff_sigma_state0(4), idx_coeff_sigma_state1(5), 
		    idx_E0(6), idx_E1(7),
		    pipi_sigma_sim_fit_Ascale(1.0){}
};

GENERATE_PARSER(InputParamArgs, INPUT_PARAM_ARGS_MEMBERS);

SARLAC_END_NAMESPACE

#endif

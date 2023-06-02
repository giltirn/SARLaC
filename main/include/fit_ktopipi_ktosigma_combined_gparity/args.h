#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_ARGS_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_ARGS_H_

#include<config.h>
#include<utils/macros.h>

#include <pipi_common/correlator_utils/threemomentum.h>

CPSFIT_START_NAMESPACE

typedef std::pair<threeMomentum, double> momMultiplicityPair;
typedef std::vector<std::string> typeFileFormat;


#define INPUT_PARAM_ARGS_MEMBERS			\
    ( std::string, kaon2pt_fit_result )						\
    ( int, idx_cK )							\
    ( int, idx_mK )							\
    ( std::string, pipi_sigma_sim_fit_result )				\
    ( int, idx_E0 )							\
    ( int, idx_E1 )							\
    ( int, idx_coeff_pipi_state0 )					\
    ( int, idx_coeff_pipi_state1 )					\
    ( int, idx_coeff_sigma_state0 )					\
    ( int, idx_coeff_sigma_state1 )					\
    ( double, pipi_sigma_sim_fit_Ascale )

struct InputParamArgs{
  GENERATE_MEMBERS(INPUT_PARAM_ARGS_MEMBERS);

  InputParamArgs(): kaon2pt_fit_result("/home/ckelly/projects/32nt64_MDWF_DSDR_GparityXYZ_fixedRNG_fullanalysis/87cfgs_DD2_stream_ext/mK/params.hdf5"),
		    idx_cK(0), idx_mK(1),
		    pipi_sigma_sim_fit_result("/home/ckelly/projects/32nt64_MDWF_DSDR_GparityXYZ_fixedRNG_fullanalysis/87cfgs_DD2_stream_ext/Epipi/I0/sim_gnd_sigma/params.hdf5"),
		    idx_coeff_pipi_state0(0), idx_coeff_pipi_state1(1), idx_coeff_sigma_state0(2), idx_coeff_sigma_state1(3), idx_E0(4), idx_E1(5),
		    pipi_sigma_sim_fit_Ascale(1.0){}
};

GENERATE_PARSER(InputParamArgs, INPUT_PARAM_ARGS_MEMBERS);

#define ARGS_MEMBERS							\
  ( std::string, data_dir )						\
  ( typeFileFormat, ktopipi_type_file_fmt )					\
  ( std::vector<momMultiplicityPair>, ktopipi_type1_pimom_proj )		\
  ( std::string, pipi_bubble_file_fmt )					\
  ( std::vector<momMultiplicityPair>, pipi_bubble_pimom_proj )	\
  ( typeFileFormat, ktosigma_type_file_fmt )				\
  ( std::string, sigma_bubble_file_fmt )				\
  ( int, Lt)								\
  ( int, tsep_pipi)							\
  ( std::vector<int>, tsep_k_pi)					\
  ( std::vector<int>, tsep_k_sigma)					\
  ( InputParamArgs, input_params )					\
  ( bool, do_simfit)						\
  ( bool, correlated )							\
  ( int, tmin_k_op)							\
  ( int, tmin_op_snk)							\
  ( int, bin_size )							\
  ( int, traj_start )							\
  ( int, traj_inc )							\
  ( int, traj_lessthan )						
  
struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): data_dir("data"), 
	  
	  ktopipi_type_file_fmt({"traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>_symmpi", 
		"traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_symmpi", 
		"traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_symmpi", 
		"traj_<TRAJ>_type4_symmpi"}),
	  
	  ktopipi_type1_pimom_proj({  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  }),
	  
	  pipi_bubble_file_fmt("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>"),

	  pipi_bubble_pimom_proj( {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
				     { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
				     { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
				     { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } } ),
	  
	  ktosigma_type_file_fmt({"traj_<TRAJ>_ktosigma_type12_deltat_<TSEP_K_SIGMA>", "traj_<TRAJ>_ktosigma_type3_deltat_<TSEP_K_SIGMA>", "traj_<TRAJ>_ktosigma_type4"}),
    
          sigma_bubble_file_fmt("traj_<TRAJ>_sigmaself_mom<MOM>_v2"),

          do_simfit(true), correlated(false), Lt(64), tsep_pipi(4), tsep_k_pi({10,12,14,16,18}), tsep_k_sigma({10,12,14,16,18}),  tmin_k_op(6), tmin_op_snk(4), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1){}

};
GENERATE_PARSER(Args, ARGS_MEMBERS);


CPSFIT_END_NAMESPACE

#endif

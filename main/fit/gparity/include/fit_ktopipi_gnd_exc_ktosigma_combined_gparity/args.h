#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_ARGS_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_ARGS_H_

#include<config.h>
#include<utils/macros.h>

#include <pipi_common/correlator_utils/threemomentum.h>

#include "enums.h"
#include "input_param_args.h"

SARLAC_START_NAMESPACE

typedef std::pair<threeMomentum, double> momMultiplicityPair;
typedef std::vector<std::string> typeFileFormat;

#define ARGS_MEMBERS							\
  ( std::string, data_dir )						\
  ( typeFileFormat, ktopipi_type_file_fmt )					\
  ( typeFileFormat, ktopipi_exc_type_file_fmt )					\
  ( std::string, pipi_bubble_file_fmt )					\
  ( std::vector<momMultiplicityPair>, ktopipi_type1_pimom_proj )		\
  ( std::vector<momMultiplicityPair>, ktopipi_exc_type1_pimom_proj )		\
  ( std::vector<momMultiplicityPair>, pipi_bubble_pimom_proj )	\
  ( std::vector<momMultiplicityPair>, pipi_exc_bubble_pimom_proj )	\
  ( typeFileFormat, ktosigma_type_file_fmt )				\
  ( std::string, sigma_bubble_file_fmt )				\
  ( std::vector<PiPiOperator>, operators )				\
  ( int, Lt)								\
  ( int, tsep_pipi)							\
  ( std::vector<int>, tsep_k_pi)					\
  ( std::vector<int>, tsep_k_sigma)					\
  ( InputParamArgs, input_params )					\
  ( SimFitFunction, fitfunc )						\
  ( int, nstate )							\
  ( Basis, basis )							\
  ( bool, correlated )							\
  ( CovarianceMatrix, covariance_matrix )				\
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

	  ktopipi_exc_type_file_fmt({"traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>_symmpi_ext", 
		"traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_symmpi_ext", 
		"traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_symmpi_ext", 
		"traj_<TRAJ>_type4_symmpi_ext"}),

	  pipi_bubble_file_fmt("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_pi1mom<PB>_pi2mom<PA>_symm"),
	  
	  ktopipi_type1_pimom_proj({  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  }),
	  
	  ktopipi_exc_type1_pimom_proj({  { {3,1,1}, 1./24 }, { {-3,1,1}, 1./24 }, { {3,-1,1}, 1./24 }, { {3,1,-1}, 1./24 }, 
					  { {-3,-1,-1}, 1./24 }, { {3,-1,-1}, 1./24 }, { {-3,1,-1}, 1./24 }, { {-3,-1,1}, 1./24 },
                                          { {1,3,1}, 1./24 }, { {1,-3,1}, 1./24 }, { {1,3,-1}, 1./24 }, { {-1,3,1}, 1./24 }, 
                                          { {-1,-3,-1}, 1./24 }, { {-1,3,-1}, 1./24 }, { {-1,-3,1}, 1./24 }, { {1,-3,-1}, 1./24 }, 
                                          { {1,1,3}, 1./24 }, { {1,1,-3}, 1./24 }, { {-1,1,3}, 1./24 }, { {1,-1,3}, 1./24 }, 
                                          { {-1,-1,-3}, 1./24 }, { {-1,-1,3}, 1./24 }, { {1,-1,-3}, 1./24 }, { {-1,1,-3}, 1./24 }   }),
         
	  pipi_bubble_pimom_proj( {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
				     { {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
				     { {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
				     { {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } } ),

	  pipi_exc_bubble_pimom_proj({  { {3,1,1}, 1./24 }, { {-3,1,1}, 1./24 }, { {3,-1,1}, 1./24 }, { {3,1,-1}, 1./24 }, 
					  { {-3,-1,-1}, 1./24 }, { {3,-1,-1}, 1./24 }, { {-3,1,-1}, 1./24 }, { {-3,-1,1}, 1./24 },
                                          { {1,3,1}, 1./24 }, { {1,-3,1}, 1./24 }, { {1,3,-1}, 1./24 }, { {-1,3,1}, 1./24 }, 
                                          { {-1,-3,-1}, 1./24 }, { {-1,3,-1}, 1./24 }, { {-1,-3,1}, 1./24 }, { {1,-3,-1}, 1./24 }, 
                                          { {1,1,3}, 1./24 }, { {1,1,-3}, 1./24 }, { {-1,1,3}, 1./24 }, { {1,-1,3}, 1./24 }, 
                                          { {-1,-1,-3}, 1./24 }, { {-1,-1,3}, 1./24 }, { {1,-1,-3}, 1./24 }, { {-1,1,-3}, 1./24 }   }),

	  ktosigma_type_file_fmt({"traj_<TRAJ>_ktosigma_type12_deltat_<TSEP_K_SIGMA>", "traj_<TRAJ>_ktosigma_type3_deltat_<TSEP_K_SIGMA>", "traj_<TRAJ>_ktosigma_type4"}),
    
          sigma_bubble_file_fmt("traj_<TRAJ>_sigmaself_mom<MOM>_v2"),
    
    operators({PiPiOperator::PiPiGnd, PiPiOperator::PiPiExc, PiPiOperator::Sigma}),

    fitfunc(SimFitFunction::MultiState),
    
    nstate(3),

    basis(Basis::Basis10),

    correlated(false), Lt(64), tsep_pipi(4), tsep_k_pi({10,12,14,16,18}), tsep_k_sigma({10,12,14,16,18}),  tmin_k_op(6), tmin_op_snk(4), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1), covariance_matrix(CovarianceMatrix::Regular){}

};
GENERATE_PARSER(Args, ARGS_MEMBERS);


SARLAC_END_NAMESPACE

#endif

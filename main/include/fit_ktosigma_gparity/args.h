#ifndef _FIT_KTOSIGMA_GPARITY_ARGS_H_
#define _FIT_KTOSIGMA_GPARITY_ARGS_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

typedef std::vector<std::string> typeFileFormat;

#define ARGS_MEMBERS							\
  ( std::string, data_dir )						\
  ( typeFileFormat, type_file_fmt )				\
  ( std::string, bubble_file_fmt )					\
  ( int, Lt)								\
  ( std::vector<int>, tsep_k_sigma)					\
  ( KtoPiPiFitFunc, fitfunc)						\
  ( bool, correlated )							\
  ( int, tmin_k_op)							\
  ( int, tmin_op_sigma)							\
  ( int, bin_size )							\
  ( int, traj_start )							\
  ( int, traj_inc )							\
  ( int, traj_lessthan )						\
  
struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): data_dir("data"), 

	  type_file_fmt({ "traj_<TRAJ>_ktosigma_type12_deltat_<TSEP_K_SIGMA>",
		"traj_<TRAJ>_ktosigma_type3_deltat_<TSEP_K_SIGMA>",
		"traj_<TRAJ>_ktosigma_type4" }),

	  bubble_file_fmt("traj_<TRAJ>_sigmaself_mom<MOM>_v2"),

	  fitfunc(KtoPiPiFitFunc::FitSeparate), correlated(false), Lt(64), tsep_k_sigma(1,10), tmin_k_op(6), tmin_op_sigma(4), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

CPSFIT_END_NAMESPACE

#endif

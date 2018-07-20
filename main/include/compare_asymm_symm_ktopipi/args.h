#ifndef _COMPARE_ASYMM_SYMM_KTOPIPI_GPARITY_ARGS_H_
#define _COMPARE_ASYMM_SYMM_KTOPIPI_GPARITY_ARGS_H_


#define COMPARE_ARGS_MEMBERS						\
  ( std::string, data_dir_asymm )					\
  ( std::string, data_dir_symm )					\
  ( int, Lt)							\
  ( int, tsep_pipi)						\
  ( std::vector<int>, tsep_k_pi)				\
  ( std::string, mK_file_asymm)					\
  ( int, mK_param_asymm)					\
  ( std::string, mK_file_symm)					\
  ( int, mK_param_symm)					\
  ( int, weighted_avg_tmin_k_op)				\
  ( int, bin_size )						\
  ( int, traj_start )						\
  ( int, traj_inc )						\
  ( int, traj_lessthan )

struct ComparisonArgs{
  GENERATE_MEMBERS(COMPARE_ARGS_MEMBERS);

  ComparisonArgs(): data_dir_asymm("data"), data_dir_symm("data"), Lt(64), tsep_pipi(4), tsep_k_pi(1,10), traj_start(0), traj_inc(1), traj_lessthan(2), mK_file_asymm("kaon_fit_asymm.hdf5"), mK_param_asymm(1), mK_file_symm("kaon_fit_symm.hdf5"), mK_param_symm(1), bin_size(1){}

};
GENERATE_PARSER(ComparisonArgs, COMPARE_ARGS_MEMBERS);


#endif

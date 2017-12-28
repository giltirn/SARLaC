#ifndef _COMPARE_ASYMM_SYMM_KTOPIPI_GPARITY_ARGS_H_
#define _COMPARE_ASYMM_SYMM_KTOPIPI_GPARITY_ARGS_H_


#define COMPARE_ARGS_MEMBERS						\
  ( std::string, data_dir_asymm )					\
  ( std::string, data_dir_symm )					\
  ( int, Lt)							\
  ( int, tsep_pipi)						\
  ( std::vector<int>, tsep_k_pi)				\
  ( int, traj_start )						\
  ( int, traj_inc )						\
  ( int, traj_lessthan )					\

struct ComparisonArgs{
  GENERATE_MEMBERS(COMPARE_ARGS_MEMBERS);

  ComparisonArgs(): data_dir_asymm("data"), data_dir_symm("data"), Lt(64), tsep_pipi(4), tsep_k_pi(1,10), traj_start(0), traj_inc(1), traj_lessthan(2){}

  Args toArgs(const AsymmSymm type) const{
    Args out;
    out.data_dir = type == Asymmetric ? data_dir_asymm : data_dir_symm;
    out.Lt = Lt;
    out.tsep_pipi = tsep_pipi;
    out.tsep_k_pi = tsep_k_pi;
    out.traj_start = traj_start;
    out.traj_inc = traj_inc;
    out.traj_lessthan = traj_lessthan;
    return out;
  }
};
GENERATE_PARSER(ComparisonArgs, COMPARE_ARGS_MEMBERS);


#endif

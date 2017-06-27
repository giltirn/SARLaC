#ifndef _FIT_PIPI_GPARITY_ARGS_H_
#define _FIT_PIPI_GPARITY_ARGS_H_

#define ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( int, Lt) \
  ( int, tsep_pipi) \
  ( int, t_min) \
  ( int, t_max) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)

  Args(): data_dir("data"), Lt(64), tsep_pipi(4), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13) {}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)

#endif

#ifndef _FIT_SIMPLE_ARGS_H_
#define _FIT_SIMPLE_ARGS_H_

#include<fit_simple/data_info.h>

#define ARGS_MEMBERS \
  ( std::vector<DataInfo>, data ) \
  ( Combination, combination ) \
  ( TimeDependence, outer_time_dep ) \
  ( bool, correlated ) \
  ( FitFuncType, fitfunc) \
  ( int, Lt) \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )


struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

Args(): Lt(64), combination(Combination::CombinationAverage), outer_time_dep(TimeDependence::TimeDepNormal), correlated(false), traj_start(0), traj_inc(1), traj_lessthan(2), t_min(0), t_max(32), data(1), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

#endif

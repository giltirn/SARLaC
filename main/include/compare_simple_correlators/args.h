#ifndef _COMPARE_SIMPLE_CORRELATORS_ARGS_H_
#define _COMPARE_SIMPLE_CORRELATORS_ARGS_H_

#include<fit_simple/data_info.h>

#define ARGS_MEMBERS \
  ( std::vector<DataInfo>, data_A ) \
  ( std::vector<DataInfo>, data_B ) \
  ( Combination, combination_A ) \
  ( Combination, combination_B ) \
  ( TimeDependence, outer_time_dep_A ) \
  ( TimeDependence, outer_time_dep_B ) \
  ( int, Lt) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): Lt(64),
    data_A(1), data_B(1),
    combination_A(CombinationAverage), combination_B(CombinationAverage),
    outer_time_dep_A(TimeDepNormal), outer_time_dep_B(TimeDepNormal),
    traj_start(0), traj_inc(1), traj_lessthan(2){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

#endif

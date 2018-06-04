#ifndef FIT_SIMULTANEOUS_ARGS_H
#define FIT_SIMULTANEOUS_ARGS_H

#define CORRINFO_MEMBERS \
  ( std::vector<DataInfo>, data ) \
  ( Combination, combination ) \
  ( TimeDependence, outer_time_dep ) \
  ( FitFuncType, fitfunc)

struct CorrInfo{
  GENERATE_MEMBERS(CORRINFO_MEMBERS);

  CorrInfo(): combination(CombinationAverage), outer_time_dep(TimeDepFold), data(1), fitfunc(FCosh){}
};
GENERATE_PARSER(CorrInfo, CORRINFO_MEMBERS);



#define ARGS_MEMBERS \
  ( std::vector<CorrInfo>, correlators ) \
  ( bool, correlated ) \
  ( int, Lt) \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )


struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

Args(): Lt(32), correlated(false), traj_start(568), traj_inc(8), traj_lessthan(912), t_min(5), t_max(15), correlators(1), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

#endif

#ifndef _FIT_SIMPLE_ARGS_H_
#define _FIT_SIMPLE_ARGS_H_

GENERATE_ENUM_AND_PARSER(ParserType, (ParserStandard)(ParserMultiSourceAverage) );
GENERATE_ENUM_AND_PARSER(TimeDependence, (TimeDepNormal)(TimeDepReflect)(TimeDepFold)(TimeDepAntiFold) );
GENERATE_ENUM_AND_PARSER(Combination, (CombinationAverage) );
GENERATE_ENUM_AND_PARSER(FitFuncType, (FCosh)(FSinh)(FExp) );

#define DATA_INFO_MEMBERS \
  ( ParserType, parser )	       \
  ( std::string, operation )   \
  ( TimeDependence, time_dep ) \
  ( std::string, file_fmt )

//file_fmt should contain a '%d' which is replaced by the trajectory index
//operation is a math expression. Use x to represent the data. Use an empty string to leave as-is

struct DataInfo{
  GENERATE_MEMBERS(DATA_INFO_MEMBERS)
};
GENERATE_PARSER(DataInfo, DATA_INFO_MEMBERS)

#define ARGS_MEMBERS \
  ( std::vector<DataInfo>, data ) \
  ( Combination, combination ) \
  ( TimeDependence, outer_time_dep ) \
  ( bool, correlated ) \
  ( FitFuncType, fitfunc) \
  ( int, Lt) \
  ( int, t_min) \
  ( int, t_max) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )


struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)

  Args(): traj_start(0), traj_inc(1), traj_lessthan(2), t_min(0), t_max(32), data(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)

#endif

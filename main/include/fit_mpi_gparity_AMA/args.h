#ifndef _FIT_MPI_GPARITY_AMA_ARGS_H
#define _FIT_MPI_GPARITY_AMA_ARGS_H

GENERATE_ENUM_AND_PARSER(DataType, (PP_LW_data)(AA_LW_data)(AP_LW_data)(PP_WW_data)(AP_WW_data) )

#define SLOPPY_EXACT_MEMBERS \
  ( std::string, sloppy_fmt )     \
  ( std::string, exact_fmt ) \
  ( bool, include_data )

struct SloppyExact{
  GENERATE_MEMBERS(SLOPPY_EXACT_MEMBERS)
  SloppyExact(): include_data(false){}
};
GENERATE_PARSER(SloppyExact,SLOPPY_EXACT_MEMBERS)

#define TWOPOINTFUNCTION_MEMBERS \
  ( DataType, type )	       \
  ( SloppyExact, FF_data )     \
  ( SloppyExact, BB_data )

struct TwoPointFunction{
  GENERATE_MEMBERS(TWOPOINTFUNCTION_MEMBERS)
};
GENERATE_PARSER(TwoPointFunction, TWOPOINTFUNCTION_MEMBERS)


#define ARGS_MEMBERS \
  ( std::vector<TwoPointFunction>, data ) \
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

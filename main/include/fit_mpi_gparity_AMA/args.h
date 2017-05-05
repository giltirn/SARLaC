#ifndef _FIT_MPI_GPARITY_AMA_ARGS_H
#define _FIT_MPI_GPARITY_AMA_ARGS_H

#define SLOPPY_EXACT_MEMBERS \
  ELEM( std::string, sloppy_fmt )     \
  ELEM( std::string, exact_fmt )      \
  ELEM( bool, include_data )

struct SloppyExact{
#define STRUCT_ARGS SLOPPY_EXACT_MEMBERS
#include<struct_gen.incl>

  SloppyExact(): include_data(false){}
};
#define STRUCT_TYPE SloppyExact
#define STRUCT_ARGS SLOPPY_EXACT_MEMBERS
#include<parser_gen.incl>


#define TWOPOINTFUNCTION_MEMBERS \
  ELEM( SloppyExact, FF_data )     \
  ELEM( SloppyExact, BB_data )

struct TwoPointFunction{
#define STRUCT_ARGS TWOPOINTFUNCTION_MEMBERS
#include<struct_gen.incl>

};
#define STRUCT_TYPE TwoPointFunction
#define STRUCT_ARGS TWOPOINTFUNCTION_MEMBERS
#include<parser_gen.incl>


#define ARGS_MEMBERS \
  ELEM( TwoPointFunction, PP_LW )   \
  ELEM( TwoPointFunction, AP_LW )   \
  ELEM( int, Lt) \
  ELEM( int, t_min) \
  ELEM( int, t_max) \
  ELEM( int, traj_start ) \
  ELEM( int, traj_inc ) \
  ELEM( int, traj_lessthan )


struct Args{
#define STRUCT_ARGS ARGS_MEMBERS
#include<struct_gen.incl>

Args(): traj_start(0), traj_inc(1), traj_lessthan(2), t_min(0), t_max(32){}
};
#define STRUCT_TYPE Args
#define STRUCT_ARGS ARGS_MEMBERS
#include<parser_gen.incl>


#endif

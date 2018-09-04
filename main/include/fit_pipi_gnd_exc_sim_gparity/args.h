#ifndef _FIT_PIPI_GND_EXC_GPARITY_ARGS_H_
#define _FIT_PIPI_GND_EXC_GPARITY_ARGS_H_

#include "enums.h"

//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
#define ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( std::string, figure_file_format )    \
  ( std::string, bubble_file_format )    \
  ( int, Lt) \
  ( int, tsep_pipi) \
  ( int, tstep_pipi) \
  ( bool, do_vacuum_subtraction ) \
  ( FitFuncType, fitfunc) \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)
  Args(): data_dir("data"), 
    figure_file_format("traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_p1src<P1SRC>_p2src<P2SRC>_p1snk<P1SNK>_p2snk<P2SNK>_symm"), 
    bubble_file_format("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_pi1mom<PB>_pi2mom<PA>_symm"),
    Lt(64), tsep_pipi(4), tstep_pipi(8), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), 
    fitfunc(FitFuncType::FSimGenOneState), do_vacuum_subtraction(true), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)


#endif

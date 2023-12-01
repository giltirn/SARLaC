#ifndef _FIT_PIPI_COMOVING_GPARITY_ARGS_H
#define _FIT_PIPI_COMOVING_GPARITY_ARGS_H

#include "enums.h"

typedef std::array<int,3> int_array_3;

//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
//Note nstate applies only for "MultiState" fit func variants
#define ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( std::string, pipi_figure_file_format )    \
  ( std::string, pipi_bubble_file_format )    \
  ( std::vector<Operator>, operators )	      \
  ( std::vector<int_array_3>, p_tot )	      \
  ( int, isospin )			      \
  ( int, Lt)				      \
  ( int, tsep_pipi) \
  ( int, tstep_pipi) \
  ( bool, do_vacuum_subtraction ) \
  ( FitFuncType, fitfunc) \
  ( int, nstate)	  \
  ( int, t_min) \
  ( int, t_max) \
  ( bool, correlated )				\
  ( int, bin_size) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)
  Args(): data_dir("data"), 
    pipi_figure_file_format("traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_p1src<P1SRC>_p2src<P2SRC>_p1snk<P1SNK>_p2snk<P2SNK>_symm"), 
    pipi_bubble_file_format("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_pi1mom<PB>_pi2mom<PA>_symm"),
    operators({Operator::PiPiComoveGnd, Operator::PiPiComoveExc1, Operator::PiPiComoveExc2}),
    p_tot({{2,0,0}}), do_vacuum_subtraction(false),
    isospin(0),Lt(64), tsep_pipi(4), tstep_pipi(8), t_min(0), t_max(32), correlated(true), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), 
    fitfunc(FitFuncType::FSimGenMultiState), nstate(3), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)


#endif

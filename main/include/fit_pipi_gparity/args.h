#ifndef _FIT_PIPI_GPARITY_ARGS_H_
#define _FIT_PIPI_GPARITY_ARGS_H_

#include "enums.h"

typedef std::array<int,3> int_array_3;

//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
#define ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( int, Lt) \
  ( int, tsep_pipi) \
  ( int, tstep_pipi) \
  ( PiPiProjector, proj_src ) \
  ( PiPiProjector, proj_snk ) \
  ( PiPiMomAllowed, allowed_mom ) \
  ( PiPiCorrSelector, corr_selector ) \
  ( std::vector<int_array_3>, pion_momenta) \
  ( std::vector<int_array_3>, total_mom )	    \
  ( int, isospin )		  \
  ( bool, do_vacuum_subtraction ) \
  ( PiPiFitFunction, fitfunc) \
  ( PiPiEffectiveEnergy, effective_energy) \
  ( bool, correlated)	      \
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

  Args(): data_dir("data"), Lt(64), tsep_pipi(4), tstep_pipi(8), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), fitfunc(PiPiFitFunction::FCoshPlusConstant), correlated(true), do_vacuum_subtraction(true), bin_size(1), proj_src(PiPiProjector::A1), proj_snk(PiPiProjector::A1), isospin(0), allowed_mom(PiPiMomAllowed::All), corr_selector(PiPiCorrSelector::Basic), total_mom({ {0,0,0} }), 
    pion_momenta({ {1,1,1}, {-1,-1,-1}, {1,1,-1}, {-1,-1,1}, {1,-1,1}, {-1,1,-1}, {-1,1,1}, {1,-1,-1} }){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)


#endif

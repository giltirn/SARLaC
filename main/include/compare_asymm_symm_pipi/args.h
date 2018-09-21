#ifndef _COMPARE_ASYMM_SYMM_PIPI_GPARITY_ARGS_H_
#define _COMPARE_ASYMM_SYMM_PIPI_GPARITY_ARGS_H_

//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
#define COMPARE_ARGS_MEMBERS \
  ( std::string, data_dir_asymm ) \
  ( std::string, data_dir_symm ) \
  ( int, Lt) \
  ( int, tsep_pipi) \
  ( int, tstep_pipi) \
  ( PiPiProjector, proj_src ) \
  ( PiPiProjector, proj_snk ) \
  ( int, isospin )		  \
  ( bool, do_vacuum_subtraction ) \
  ( int, bin_size ) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \

struct ComparisonArgs{
  GENERATE_MEMBERS(COMPARE_ARGS_MEMBERS)

  ComparisonArgs(): data_dir_asymm("data"), data_dir_symm("data"), Lt(64), tsep_pipi(4), tstep_pipi(8), traj_start(0), traj_inc(1), traj_lessthan(2), do_vacuum_subtraction(true), bin_size(1), proj_src(PiPiProjector::A1momSet111), proj_snk(PiPiProjector::A1momSet111), isospin(0) {}

};
GENERATE_PARSER(ComparisonArgs, COMPARE_ARGS_MEMBERS)

#endif

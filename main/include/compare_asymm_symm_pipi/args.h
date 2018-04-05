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

  ComparisonArgs(): data_dir_asymm("data"), data_dir_symm("data"), Lt(64), tsep_pipi(4), tstep_pipi(8), traj_start(0), traj_inc(1), traj_lessthan(2), do_vacuum_subtraction(true), bin_size(1), proj_src(A1), proj_snk(A1), isospin(0)  {}

  Args toArgs(const AsymmSymm type) const{
    Args out;
    out.data_dir = type == Asymmetric ? data_dir_asymm : data_dir_symm;
    out.Lt = Lt;
    out.tsep_pipi = tsep_pipi;
    out.tstep_pipi = tstep_pipi;
    out.proj_src = proj_src;
    out.proj_snk = proj_snk;
    out.isospin = isospin;
    out.do_vacuum_subtraction = do_vacuum_subtraction;
    out.bin_size = bin_size;
    out.traj_start = traj_start;
    out.traj_inc = traj_inc;
    out.traj_lessthan = traj_lessthan;
    return out;
  }
  
};
GENERATE_PARSER(ComparisonArgs, COMPARE_ARGS_MEMBERS)

#endif

#ifndef _PIPI_SIGMA_SIM_FIT_ARGS_H_
#define _PIPI_SIGMA_SIM_FIT_ARGS_H_

#define PIPI_SIGMA_SIM_ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( int, Lt) \
  ( int, tsep_pipi ) \
  ( int, tstep_pipi2pt )	  \
  ( bool, do_vacuum_subtraction ) \
  ( bool, correlated)	      \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale)

struct PiPiSigmaSimArgs{
  GENERATE_MEMBERS( PIPI_SIGMA_SIM_ARGS_MEMBERS)

  PiPiSigmaSimArgs(): data_dir("data"), Lt(64), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), correlated(true), do_vacuum_subtraction(true), bin_size(1), tsep_pipi(4),tstep_pipi2pt(8){}

  void transfer(SimFitArgs &fargs) const{
    fargs.correlated = correlated;
    fargs.Lt = Lt;
    fargs.tsep_pipi = tsep_pipi;
    fargs.Ascale = Ascale;
    fargs.Cscale = Cscale;
  }    
};
GENERATE_PARSER(PiPiSigmaSimArgs, PIPI_SIGMA_SIM_ARGS_MEMBERS);

#endif

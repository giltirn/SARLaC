#ifndef _PIPI_TO_SIGMA_ARGS_H_
#define _PIPI_TO_SIGMA_ARGS_H_

#define PIPI_TO_SIGMA_ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( int, Lt) \
  ( int, tsep_pipi ) \
  ( bool, do_vacuum_subtraction ) \
  ( PiPiToSigmaFitFunction, fitfunc) \
  ( bool, correlated)	      \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale)

struct PiPiToSigmaArgs{
  GENERATE_MEMBERS(PIPI_TO_SIGMA_ARGS_MEMBERS)

  PiPiToSigmaArgs(): data_dir("data"), Lt(64), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), fitfunc(FCoshPlusConstant), correlated(true), do_vacuum_subtraction(true), bin_size(1), tsep_pipi(4){}

  void transfer(PiPiToSigmaFitArgs &fargs) const{
    fargs.fitfunc = fitfunc;  //FCoshPlusConstant;
    fargs.Lt = Lt;
    fargs.tsep_pipi = tsep_pipi;
    fargs.Ascale = Ascale;
    fargs.Cscale = Cscale;
    fargs.correlated = correlated;    
  }
};
GENERATE_PARSER(PiPiToSigmaArgs, PIPI_TO_SIGMA_ARGS_MEMBERS);


#endif

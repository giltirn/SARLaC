#ifndef _SIGMA_ARGS_H_
#define _SIGMA_ARGS_H_

#define SIGMA_ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( int, Lt) \
  ( bool, do_vacuum_subtraction ) \
  ( SigmaFitFunction, fitfunc) \
  ( bool, correlated)	      \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size) \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale)

struct SigmaArgs{
  GENERATE_MEMBERS(SIGMA_ARGS_MEMBERS)

  SigmaArgs(): data_dir("data"), Lt(64), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), fitfunc(SigmaFitFunction::FCoshPlusConstant), correlated(true), do_vacuum_subtraction(true), bin_size(1){}

  void transfer(SigmaFitArgs &fargs) const{
    fargs.fitfunc = fitfunc;  //FCoshPlusConstant;
    fargs.Lt = Lt;
    fargs.Ascale = Ascale;
    fargs.Cscale = Cscale;
    fargs.correlated = correlated;
  }
};
GENERATE_PARSER(SigmaArgs, SIGMA_ARGS_MEMBERS);

#endif

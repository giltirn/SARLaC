#ifndef _PIPI_SIGMA_SIM_FIT_SAMPLEAMA_ARGS_H_
#define _PIPI_SIGMA_SIM_FIT_SAMPLEAMA_ARGS_H_

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

#define ARGS_MEMBERS \
  ( std::string, config_map )			\
  ( int, Lt) \
  ( int, tsep_pipi ) \
  ( int, tstep_pipi2pt )	  \
  ( int, tstep_pipitosigma )	  \
  ( bool, do_vacuum_subtraction ) \
  ( bool, correlated)	      \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size) \
  ( double, Ascale) \
  ( double, Cscale)

struct Args{
  GENERATE_MEMBERS( ARGS_MEMBERS)

  Args(): config_map("config_map.dat"), Lt(64), t_min(0), t_max(32) , Ascale(1), Cscale(1), 
    correlated(true), do_vacuum_subtraction(true), bin_size(1), tsep_pipi(4), tstep_pipi2pt(8), tstep_pipitosigma(1)
  {}

  void transfer(SimFitArgs &fargs) const{
    fargs.correlated = correlated;
    fargs.Lt = Lt;
    fargs.tsep_pipi = tsep_pipi;
    fargs.Ascale = Ascale;
    fargs.Cscale = Cscale;
  }    
};
GENERATE_PARSER(Args, ARGS_MEMBERS);

SARLAC_END_NAMESPACE

#endif

#ifndef _FIT_PIPI_GPARITY_SAMPLEAMA_ARGS_H_
#define _FIT_PIPI_GPARITY_SAMPLEAMA_ARGS_H_


//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
#define ARGS_SAMPLEAMA_MEMBERS \
  ( std::string, data_dir_S ) \
  ( std::string, data_dir_C ) \
  ( int, Lt) \
  ( int, tsep_pipi) \
  ( int, tstep_pipi) \
  ( PiPiProjector, proj_src ) \
  ( PiPiProjector, proj_snk ) \
  ( PiPiMomAllowed, allowed_mom ) \
  ( int, isospin )		  \
  ( bool, do_vacuum_subtraction ) \
  ( PiPiFitFunction, fitfunc) \
  ( PiPiEffectiveEnergy, effective_energy) \
  ( bool, correlated)	      \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size) \
  ( int, traj_inc ) \
  ( int, traj_start_S )	 \
  ( int, traj_lessthan_S ) \
  ( int, traj_start_C )	 \
  ( int, traj_lessthan_C ) \
  ( double, Ascale) \
  ( double, Cscale)

struct ArgsSampleAMA{
  GENERATE_MEMBERS(ARGS_SAMPLEAMA_MEMBERS)

  ArgsSampleAMA(): data_dir_S("data_S"), data_dir_C("data_C"), Lt(64), tsep_pipi(4), tstep_pipi(8), t_min(0), t_max(32), traj_inc(1), Ascale(1e13), Cscale(1e13), fitfunc(PiPiFitFunction::FCoshPlusConstant), correlated(true), do_vacuum_subtraction(true), bin_size(1), traj_start_S(0), traj_lessthan_S(2), traj_start_C(0), traj_lessthan_C(2), proj_src(PiPiProjector::A1), proj_snk(PiPiProjector::A1), allowed_mom(PiPiMomAllowed::All), isospin(0)  {}

  inline const std::string & data_dir(const char ens) const{ return ens == 'S' ? data_dir_S :  data_dir_C; };

  inline void traj_info(int &traj_start_, int &traj_inc_, int &traj_lessthan_, const char ens) const{
    traj_start_ = ens == 'S' ? traj_start_S : traj_start_C;
    traj_lessthan_ = ens == 'S' ? traj_lessthan_S : traj_lessthan_C;
    traj_inc_ = traj_inc;
  }

};
GENERATE_PARSER(ArgsSampleAMA, ARGS_SAMPLEAMA_MEMBERS)

#endif

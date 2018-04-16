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

  ArgsSampleAMA(): data_dir_S("data_S"), data_dir_C("data_C"), Lt(64), tsep_pipi(4), tstep_pipi(8), t_min(0), t_max(32), traj_inc(1), Ascale(1e13), Cscale(1e13), fitfunc(FCoshPlusConstant), correlated(true), do_vacuum_subtraction(true), bin_size(1), traj_start_S(0), traj_lessthan_S(2), traj_start_C(0), traj_lessthan_C(2), proj_src(A1), proj_snk(A1), allowed_mom(All), isospin(0)  {}

  Args toArgs(const char ens) const{  //'S' or 'C'
    Args out;
    out.Lt = Lt;
    out.tsep_pipi = tsep_pipi;
    out.tstep_pipi = tstep_pipi;
    out.proj_src = proj_src;
    out.proj_snk = proj_snk;
    out.allowed_mom = allowed_mom;
    out.isospin = isospin;
    out.do_vacuum_subtraction = do_vacuum_subtraction;
    out.fitfunc = fitfunc;
    out.effective_energy = effective_energy;
    out.correlated = correlated;
    out.t_min = t_min;
    out.t_max = t_max;
    out.bin_size = bin_size;
    out.traj_inc = traj_inc;
    out.Ascale = Ascale;
    out.Cscale = Cscale;
    if(ens == 'S'){
      out.data_dir = data_dir_S;
      out.traj_start = traj_start_S;
      out.traj_lessthan = traj_lessthan_S;
    }else{
      out.data_dir = data_dir_C;
      out.traj_start = traj_start_C;
      out.traj_lessthan = traj_lessthan_C;
    }
    return out;
  }
};
GENERATE_PARSER(ArgsSampleAMA, ARGS_SAMPLEAMA_MEMBERS)

#endif

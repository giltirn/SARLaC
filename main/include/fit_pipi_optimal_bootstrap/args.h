#ifndef _FIT_PIPI_OPTIMAL_ARGS_BOOTSTRAP_H_
#define _FIT_PIPI_OPTIMAL_ARGS_BOOTSTRAP_H_

#include "fit_pipi_gnd_exc_sigma_sim_gparity/enums.h"

//idx are the indices for the couplings of this operator to each state
#define PARAM_ELEM_MEMBERS \
  (std::string, file)				\
  (std::vector<int>, idx)

struct ParamElem{
  GENERATE_MEMBERS(PARAM_ELEM_MEMBERS)
  ParamElem(): file("file.hdf5"), idx({0,1,2}){}
};
GENERATE_PARSER(ParamElem, PARAM_ELEM_MEMBERS)


//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
#define ARGS_MEMBERS \
  ( std::string, data_dir ) \
  ( std::string, pipi_figure_file_format )    \
  ( std::string, pipi_bubble_file_format )    \
  ( std::string, pipi_to_sigma_file_format )  \
  ( std::string, sigma_bubble_file_format )   \
  ( std::string, sigma2pt_file_format )	      \
  ( std::vector<Operator>, operators )	      \
  ( int, Lt) \
  ( int, tsep_pipi) \
  ( int, tstep_pipi) \
  ( int, tstep_pipi_to_sigma ) \
  ( bool, do_vacuum_subtraction ) \
  ( bool, timeslice_avg_vac_sub ) \
  ( int, t_min) \
  ( int, t_max) \
  ( BootResampleTableType, resample_table_type) \
  ( int, nboot)	      \
  ( int, block_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan ) \
  ( double, Ascale) \
  ( double, Cscale) \
  ( std::vector<ParamElem>, op_amplitudes )

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)
  Args(): data_dir("data"), 
    pipi_figure_file_format("traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_p1src<P1SRC>_p2src<P2SRC>_p1snk<P1SNK>_p2snk<P2SNK>_symm"), 
    pipi_bubble_file_format("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_pi1mom<PB>_pi2mom<PA>_symm"),
    pipi_to_sigma_file_format("traj_<CONF>_pipitosigma_sigmawdagmom<MOM_QUARK_SIGMA>_pionmom<MOM_PI>_v2"),
    sigma_bubble_file_format("traj_<CONF>_sigmaself_mom<PQUARK>_v2"),
    sigma2pt_file_format("traj_<CONF>_sigmacorr_mompsrc<PSRC_QUARK>psnk<PSNK_QUARK>_v2"),
    operators({Operator::PiPiGnd, Operator::PiPiExc, Operator::Sigma}),
    Lt(64), tsep_pipi(4), tstep_pipi(8), tstep_pipi_to_sigma(8), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), 
    do_vacuum_subtraction(true), timeslice_avg_vac_sub(false), 
    block_size(1), nboot(1000), resample_table_type(BootResampleTableType::NonOverlappingBlock), op_amplitudes(3){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)


#endif

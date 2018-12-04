#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_ARGS_H_
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_ARGS_H_

GENERATE_ENUM_AND_PARSER(Operator, (PiPiGnd)(PiPiExc)(Sigma) );
void write(HDF5writer &wr, const Operator op, const std::string &nm){
  int p = (int)op;
  _force_external_lookup<int>::fwrite(wr, p, nm);
}
void read(HDF5reader &rd, Operator &op, const std::string &nm){
  int p;
  _force_external_lookup<int>::fread(rd, p, nm);
  op = (Operator)p;
}


//Note tstep_pipi is the separation between source timeslices that the C, D, R diagrams were measured upon (i.e. every 8 in the 32^3 job)
//Note nstate applies only for "MultiState" fit func variants
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
  ( FitFuncType, fitfunc) \
  ( int, nstate)	  \
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
  Args(): data_dir("data"), 
    pipi_figure_file_format("traj_<TRAJ>_Figure<FIG>_sep<TSEP_PIPI>_p1src<P1SRC>_p2src<P2SRC>_p1snk<P1SNK>_p2snk<P2SNK>_symm"), 
    pipi_bubble_file_format("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_pi1mom<PB>_pi2mom<PA>_symm"),
    pipi_to_sigma_file_format("traj_<CONF>_pipitosigma_sigmawdagmom<MOM_QUARK_SIGMA>_pionmom<MOM_PI>_v2"),
    sigma_bubble_file_format("traj_<CONF>_sigmaself_mom<PQUARK>_v2"),
    sigma2pt_file_format("traj_<CONF>_sigmacorr_mompsrc<PSRC_QUARK>psnk<PSNK_QUARK>_v2"),
    operators({Operator::PiPiGnd, Operator::PiPiExc, Operator::Sigma}),
    Lt(64), tsep_pipi(4), tstep_pipi(8), tstep_pipi_to_sigma(8), t_min(0), t_max(32), traj_start(0), traj_inc(1), traj_lessthan(2), Ascale(1e13), Cscale(1e13), 
    fitfunc(FitFuncType::FSimGenTwoState), nstate(3), do_vacuum_subtraction(true), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)


#endif

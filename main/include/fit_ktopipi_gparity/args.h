#ifndef _FIT_KTOPIPI_GPARITY_ARGS_H_
#define _FIT_KTOPIPI_GPARITY_ARGS_H_

typedef std::pair<threeMomentum, double> momMultiplicityPair;
typedef std::vector<std::string> typeFileFormat;

#define ARGS_MEMBERS							\
  ( std::string, data_dir )						\
  ( typeFileFormat, type_file_fmt )					\
  ( std::vector<momMultiplicityPair>, type1_pimom_proj )		\
  ( std::string, bubble_file_fmt )					\
  ( std::vector<momMultiplicityPair>, bubble_pimom_proj )	\
  ( int, Lt)								\
  ( int, tsep_pipi)							\
  ( std::vector<int>, tsep_k_pi)					\
  ( KtoPiPiFitFunc, fitfunc)						\
  ( bool, correlated )							\
  ( int, tmin_k_op)							\
  ( int, tmin_op_pi)							\
  ( int, bin_size )							\
  ( int, traj_start )							\
  ( int, traj_inc )							\
  ( int, traj_lessthan )						\
  
struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): data_dir("data"), 

	  type_file_fmt({"traj_<TRAJ>_type1_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>_mom<MOM>",
		"traj_<TRAJ>_type2_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
		"traj_<TRAJ>_type3_deltat_<TSEP_K_PI>_sep_<TSEP_PIPI>",
		"traj_<TRAJ>_type4" }),

	  type1_pimom_proj({  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  { {-1,1,1}, 3.0/8.0 }, { {1,-1,-1}, 3.0/8.0 }  }),

	  bubble_file_fmt("traj_<TRAJ>_FigureVdis_sep<TSEP_PIPI>_mom<PB>"),

	  bubble_pimom_proj( {  { {1,1,1}, 1.0/8.0 }, { {-1,-1,-1}, 1.0/8.0 },  
				{ {-1,1,1}, 1.0/8.0 }, { {1,-1,-1}, 1.0/8.0 }, 
				{ {1,-1,1}, 1.0/8.0 }, { {-1,1,-1}, 1.0/8.0 }, 
				{ {1,1,-1}, 1.0/8.0 }, { {-1,-1,1}, 1.0/8.0 } } ),

	  fitfunc(KtoPiPiFitFunc::FitSeparate), correlated(false), Lt(64), tsep_pipi(4), tsep_k_pi(1,10), tmin_k_op(6), tmin_op_pi(4), traj_start(0), traj_inc(1), traj_lessthan(2), bin_size(1){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


#endif

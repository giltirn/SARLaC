#ifndef _FIT_KTOPIPI_GPARITY_ARGS_H_
#define _FIT_KTOPIPI_GPARITY_ARGS_H_

GENERATE_ENUM_AND_PARSER(KtoPiPiFitFunc, (FitSeparate)(FitSimultaneous)(FitSimultaneousChiralBasis) );

#define ARGS_MEMBERS						\
  ( std::string, data_dir )					\
  ( int, Lt)							\
  ( int, tsep_pipi)						\
  ( std::vector<int>, tsep_k_pi)				\
  ( KtoPiPiFitFunc, fitfunc)					\
  ( bool, correlated )						\
  ( int, tmin_k_op)						\
  ( int, tmin_op_pi)						\
  ( int, traj_start )						\
  ( int, traj_inc )						\
  ( int, traj_lessthan )					\

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): data_dir("data"), fitfunc(FitSeparate), correlated(false), Lt(64), tsep_pipi(4), tsep_k_pi(1,10), tmin_k_op(6), tmin_op_pi(4), traj_start(0), traj_inc(1), traj_lessthan(2){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


GENERATE_ENUM_AND_PARSER(FreezeDataReaderType, (UKfitXMLvectorReader) );

#define FREEZE_PARAM_MEMBERS \
  (std::vector<int>, Qlist) \
  (int, param_idx) \
  (FreezeDataReaderType, reader) \
  (std::string, filename) \
  (int, input_idx) \
  (std::string, operation)

struct FreezeParam{
  //If Qlist is left empty it is assumed the freeze is the same for all Q
  //operation can be any math expression. Use the variable 'x' to represent the data.  eg operation = sqrt(1e13 * x).  Use an empty string for no operation to be performed
  GENERATE_MEMBERS(FREEZE_PARAM_MEMBERS);

  FreezeParam(): param_idx(0), reader(UKfitXMLvectorReader), filename("file.dat"), input_idx(0), operation(""), Qlist(0){}
};
GENERATE_PARSER(FreezeParam, FREEZE_PARAM_MEMBERS);


#define FREEZE_PARAMS_MEMBERS			\
  (std::vector<FreezeParam>, sources)

struct FreezeParams{
  GENERATE_MEMBERS(FREEZE_PARAMS_MEMBERS);
  
  FreezeParams(): sources(1){}
};
GENERATE_PARSER(FreezeParams, FREEZE_PARAMS_MEMBERS);


#endif

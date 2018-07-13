#ifndef FIT_PIPI_KTOPIPI_SIM_GPARITY_ARGS_H_
#define FIT_PIPI_KTOPIPI_SIM_GPARITY_ARGS_H_

GENERATE_ENUM_AND_PARSER(CorrelationStatus, (Uncorrelated)(Correlated)(PartiallyCorrelated) );
GENERATE_ENUM_AND_PARSER(SimFitFunc, (OneExp)(TwoExpPiPi) ); // TwoExpPiPi = pipi fit form has excited state, 3pt fit form does not

#define ARGS_MEMBERS \
  ( std::string, pipi_combined_data_file ) \
  ( std::string, ktopipi_amplitude_data_file )\
  ( std::string, freeze_file )	\
  ( int, Lt) \
  ( SimFitFunc, fitfunc) \
  ( CorrelationStatus, correlation_status) \
  ( int, tsep_pipi) \
  ( int, tmin_pipi) \
  ( int, tmax_pipi) \
  ( int, tmin_k_op) \
  ( int, tmin_op_pi) \
  ( double, Ascale) \
  ( double, Cscale)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS)

  Args(): pipi_combined_data_file("pipi_data.hdf5"), ktopipi_amplitude_data_file("ktopipi_data.hdf5"), freeze_file("freeze.hdf5"), Lt(64), fitfunc(SimFitFunc::OneExp), correlation_status(CorrelationStatus::Uncorrelated), tsep_pipi(4), tmin_pipi(6), tmax_pipi(25), tmin_k_op(6), tmin_op_pi(4), Ascale(1e13), Cscale(1e13){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS)



#endif

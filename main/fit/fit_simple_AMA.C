#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>
#include<random.h>

using namespace SARLaC;

#include <fit_simple_AMA/parse_data.h>
#include <fit_simple_AMA/data_manipulations.h>
#include<fit_simple/args.h>
#include<fit_simple/cmdline.h>
#include<fit_simple/read_data.h>
#include<fit_simple/fit.h>
#include<fit_simple/main.h>

#define AMA_DATA_INFO_MEMBERS \
  ( ReIm, reim)		       \
  ( std::string, operation )   \
  ( TimeDependence, time_dep ) \
  ( AMAparserType, parser ) \
  ( std::string, sloppy_file_fmt ) \
  ( std::string, exact_file_fmt ) \
  ( int, nexact_timeslice ) \
  ( bool, have_exact_timeslice_multiplicity_file) \
  ( std::string, exact_timeslice_multiplicity_file_fmt )


//file_fmt : should contain a '%d' which is replaced by the trajectory index
//operation : a math expression. Use x to represent the data. Use an empty string to leave as-is
//nexact_timeslice : the number of exact timeslices
//have_exact_timeslice_multiplicity_file : if more than one exact timeslice is used and are randomly chosen, there is a chance they may overlap. 
//                                         To allow for this, we can specify a timeslice "multiplicity" file containing a single line of space-separated values of the timeslice multiplicity. 
//                                         Their sum should add to nexact_timeslice.
//                                         If a multiplicity file is not included and nexact_timeslice > 1, multiplicities will be randomly assigned such as to add to nexact_timeslice

struct AMAdataInfo{
  GENERATE_MEMBERS(AMA_DATA_INFO_MEMBERS);
  AMAdataInfo(): operation(""), time_dep(TimeDependence::TimeDepNormal), sloppy_file_fmt("sloppy_data.%d"), exact_file_fmt("exact_data.%d"), nexact_timeslice(2), have_exact_timeslice_multiplicity_file(true), exact_timeslice_multiplicity_file_fmt("multiplicity.%d"), parser(AMAparserType::ParserReal){}
};
GENERATE_PARSER(AMAdataInfo, AMA_DATA_INFO_MEMBERS);

void readData(rawDataCorrelationFunctionD &into, const AMAdataInfo &data_info, const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan){
  rawDataDistributionVector v = readData(traj_start, traj_inc, traj_lessthan, data_info.sloppy_file_fmt, data_info.exact_file_fmt, 
					 data_info.have_exact_timeslice_multiplicity_file, data_info.exact_timeslice_multiplicity_file_fmt, data_info.nexact_timeslice,
					 data_info.reim, Lt, data_info.parser);
  into = rawDataCorrelationFunctionD(Lt, [&](const int t){ return rawDataCorrelationFunctionD::ElementType(t, std::move(v[t])); });
  applyOperation(into, data_info.operation);
  applyTimeDep(into, data_info.time_dep,Lt);
}

#define AMA_ARGS_MEMBERS \
  ( std::vector<AMAdataInfo>, data ) \
  ( Combination, combination ) \
  ( TimeDependence, outer_time_dep ) \
  ( CovarianceStrategy, covariance_strategy )	\
  ( FitFuncType, fitfunc) \
  ( int, Lt) \
  ( int, t_min) \
  ( int, t_max) \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )


struct AMAargs{
  GENERATE_MEMBERS(AMA_ARGS_MEMBERS);

  AMAargs(): Lt(64), combination(Combination::CombinationAverage), outer_time_dep(TimeDependence::TimeDepNormal), covariance_strategy(CovarianceStrategy::Uncorrelated), traj_start(0), traj_inc(1), traj_lessthan(2), t_min(0), t_max(32), data(1), bin_size(1){}
};
GENERATE_PARSER(AMAargs, AMA_ARGS_MEMBERS);




int main(const int argc, const char** argv){
  CMDline cmdline(argc,argv,2);

  AMAargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  const int nchannel = args.data.size();
  std::vector<doubleJackknifeCorrelationFunctionD> channels_dj(nchannel);
  std::vector<jackknifeCorrelationFunctionD> channels_j(nchannel);
  for(int i=0;i<nchannel;i++){
    rawDataCorrelationFunctionD channel_raw;
    readData(channel_raw, args.data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    bin(channel_raw, args.bin_size);    

    channels_dj[i] = doubleJackknifeCorrelationFunctionD(args.Lt, [&](const int t){
	return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(channel_raw.value(t)));
      });
    channels_j[i] = jackknifeCorrelationFunctionD(args.Lt, [&](const int t){
	return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(channel_raw.value(t)));
      });
  }
  doubleJackknifeCorrelationFunctionD data_dj;
  applyCombination(data_dj,channels_dj,args.combination);
  applyTimeDep(data_dj, args.outer_time_dep, args.Lt);

  jackknifeCorrelationFunctionD data_j;
  applyCombination(data_j,channels_j,args.combination);
  applyTimeDep(data_j, args.outer_time_dep, args.Lt);

  blockDoubleJackknifeCorrelationFunctionD data_bj;
  assert(args.covariance_strategy != CovarianceStrategy::CorrelatedBlockHybrid &&
	 args.covariance_strategy != CovarianceStrategy::CorrelatedBlock);

  jackknifeDistribution<parameterVectorD> params;
  jackknifeDistributionD chisq;
  int dof;
  fit(params, chisq, dof, data_j,data_dj,data_bj, args, cmdline);
  

  return 0;
}

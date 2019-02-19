#include<common.h>
#include<fit.h>
#include<parser.h>
#include<plot.h>

using namespace CPSfit;

#include <fit_mpi_gparity_AMA/parse_data.h>
#include <fit_mpi_gparity_AMA/data_manipulations.h>
#include<fit_simple/args.h>
#include<fit_simple/cmdline.h>
#include<fit_simple/read_data.h>
#include<fit_simple/fit.h>
#include<fit_simple/main.h>

rawDataDistributionVector readData(const int traj_start, const int traj_inc, const int traj_lessthan,
				   const std::string &sloppy_fmt, const std::string &exact_fmt, const ReIm reim, const int Lt){
  const int ntraj = (traj_lessthan - traj_start)/traj_inc;
  assert(ntraj > 0);

  rawDataDistributionVector corrected;
    
  rawDataDistributionMatrix exact_data(Lt, rawDataDistributionD(ntraj));
  rawDataDistributionMatrix sloppy_data(Lt, rawDataDistributionD(ntraj));

#pragma omp parallel for
  for(int i=0;i<ntraj;i++){
    const int c = traj_start + i*traj_inc;
    read(exact_data, i, exact_fmt, c, reim);
    read(sloppy_data, i, sloppy_fmt, c, reim);
  }
  
  rawDataDistributionVector sloppy_avg = sourceTimeSliceAverage(sloppy_data);
  rawDataDistributionVector correction = computeAMAcorrection(sloppy_data, exact_data);
  
  corrected = sloppy_avg + correction;

  std::cout << "Sloppy data:\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << sloppy_avg[t] << std::endl;
  
  std::cout << "Corrected data:\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << corrected[t] << std::endl;
  
  return corrected;
}

#define AMA_DATA_INFO_MEMBERS \
  ( ReIm, reim)		       \
  ( std::string, operation )   \
  ( TimeDependence, time_dep ) \
  ( std::string, sloppy_file_fmt ) \
  ( std::string, exact_file_fmt )

//file_fmt should contain a '%d' which is replaced by the trajectory index
//operation is a math expression. Use x to represent the data. Use an empty string to leave as-is

struct AMAdataInfo{
  GENERATE_MEMBERS(AMA_DATA_INFO_MEMBERS);
  AMAdataInfo(): operation(""), time_dep(TimeDependence::TimeDepNormal), sloppy_file_fmt("sloppy_data.%d"), exact_file_fmt("exact_data.%d"){}
};
GENERATE_PARSER(AMAdataInfo, AMA_DATA_INFO_MEMBERS);

void readData(rawDataCorrelationFunctionD &into, const AMAdataInfo &data_info, const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan){
  rawDataDistributionVector v = readData(traj_start, traj_inc, traj_lessthan, data_info.sloppy_file_fmt, data_info.exact_file_fmt, data_info.reim, Lt);
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

  fit(data_j,data_dj, args, cmdline);
  
  return 0;
}

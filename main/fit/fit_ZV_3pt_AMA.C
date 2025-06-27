//Fit ZV from the vector 3pt function using the method described in https://arxiv.org/pdf/1411.7017  page 34
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

struct ThreePointDataInfo{
#define THREEPOINTDATAINFO_MEMBERS    \
  (int, tsep_src_snk)		      \
  (int, fit_tmin) \
  (int, fit_tmax) \
  (AMAdataInfo, data) 
  
  GENERATE_MEMBERS(THREEPOINTDATAINFO_MEMBERS);
  ThreePointDataInfo(): tsep_src_snk(10), fit_tmin(2), fit_tmax(8){}
};
GENERATE_PARSER(ThreePointDataInfo, THREEPOINTDATAINFO_MEMBERS);

#define AMA_ARGS_MEMBERS \
  ( std::vector<ThreePointDataInfo>, C_PVP_data ) \
  ( AMAdataInfo, C_PP_WW_data )			  \
  ( std::string, pion_mass_filename )		  \
  ( int, pion_mass_param_idx )			  \
  ( CostType, covariance_strategy )	\
  ( int, Lt) \
  ( int, bin_size)    \
  ( int, traj_start ) \
  ( int, traj_inc ) \
  ( int, traj_lessthan )
  
struct AMAargs{
  GENERATE_MEMBERS(AMA_ARGS_MEMBERS);

    AMAargs(): Lt(64), covariance_strategy(CostType::Uncorrelated), traj_start(0), traj_inc(1), traj_lessthan(2), C_PVP_data(1), pion_mass_filename("pion_mass.hdf5"), pion_mass_param_idx(0), bin_size(1){}
};
GENERATE_PARSER(AMAargs, AMA_ARGS_MEMBERS);

int main(const int argc, const char** argv){
  //CMDline cmdline(argc,argv,2);

  AMAargs args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);

  doubleJackknifeCorrelationFunctionD C_PP_WW_dj;
  jackknifeCorrelationFunctionD C_PP_WW_j;
  {
    rawDataCorrelationFunctionD channel_raw;
    readData(channel_raw, args.C_PP_WW_data, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    bin(channel_raw, args.bin_size);    
    
    C_PP_WW_dj = doubleJackknifeCorrelationFunctionD(args.Lt, [&](const int t){
  	return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(channel_raw.value(t)));
      });
    C_PP_WW_j = jackknifeCorrelationFunctionD(args.Lt, [&](const int t){
  	return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(channel_raw.value(t)));
      });
  }
  
  std::vector<jackknifeDistributionD> pion_params;
  readParamsStandard(pion_params, args.pion_mass_filename);
  jackknifeDistributionD mpi = pion_params[args.pion_mass_param_idx];
  
  std::cout << "Pion mass " << mpi << std::endl;
  
  int nsample = C_PP_WW_j.value(0).size();  
  
  int Lt = args.Lt;
  assert(Lt % 2 == 0);

  jackknifeCorrelationFunctionD C_PP_WW_j_corrected = C_PP_WW_j;
  doubleJackknifeCorrelationFunctionD C_PP_WW_dj_corrected = C_PP_WW_dj;

  //Need a double-jackknife compatible pion mass. We don't consider variations within a resampled ensemble
  doubleJackknifeDistributionD mpi_dj(nsample);
  for(int s=0;s<nsample;s++)
    for(int ss=0;ss<nsample-1;ss++)
      mpi_dj.sample(s).sample(ss) = mpi.sample(s);

  //Perform the ATW correction
  for(int t=0;t<Lt;t++){
    C_PP_WW_j_corrected.value(t) = C_PP_WW_j.value(t) - 0.5 * C_PP_WW_j.value(Lt/2) * exp(-mpi * (Lt/2-t) );
    jackknifeDistributionD rat = C_PP_WW_j_corrected.value(t)/C_PP_WW_j.value(t);
    std::cout << "t=" << t << " C_PP_WW corrected " << C_PP_WW_j_corrected.value(t) << " vs uncorrected " << C_PP_WW_j.value(t) << ", ratio " << rat << std::endl;

    C_PP_WW_dj_corrected.value(t) = C_PP_WW_dj.value(t) - 0.5 * C_PP_WW_dj.value(Lt/2) * exp(-mpi_dj * (Lt/2-t) );    
  }
  
  int ntsep = args.C_PVP_data.size();
  std::vector<doubleJackknifeCorrelationFunctionD> C_PVP_dj(ntsep);
  std::vector<jackknifeCorrelationFunctionD> C_PVP_j(ntsep);;
  for(int i=0;i<ntsep;i++){
    std::cout << "Reading data for tsep = " << args.C_PVP_data[i].tsep_src_snk << std::endl;
    rawDataCorrelationFunctionD channel_raw;
    readData(channel_raw, args.C_PVP_data[i].data, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
    bin(channel_raw, args.bin_size);    
    
    C_PVP_dj[i] = doubleJackknifeCorrelationFunctionD(args.Lt, [&](const int t){
  	return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(channel_raw.value(t)));
      });
    C_PVP_j[i] = jackknifeCorrelationFunctionD(args.Lt, [&](const int t){
  	return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(channel_raw.value(t)));
      });
  }

  //Get the fit data

  doubleJackknifeCorrelationFunctionD fit_data_dj;
  jackknifeCorrelationFunctionD fit_data_j;
  int idx=0;

  for(int i=0;i<ntsep;i++){
    int tsep = args.C_PVP_data[i].tsep_src_snk;
    int tmin = args.C_PVP_data[i].fit_tmin;
    int tmax = args.C_PVP_data[i].fit_tmax;
    
    for(int t=0;t<tsep;t++){
      jackknifeDistributionD rat_cor = C_PP_WW_j_corrected.value(tsep) / C_PVP_j[i].value(t);
      jackknifeDistributionD rat_uncor = C_PP_WW_j.value(tsep) / C_PVP_j[i].value(t);
      std::cout << tsep << (t>=tmin && t<=tmax ? " |" : " ") << t << " " << rat_cor << " " << rat_uncor << std::endl;

      doubleJackknifeDistributionD rat_cor_dj = C_PP_WW_dj_corrected.value(tsep) / C_PVP_dj[i].value(t);

      if(t>=tmin && t<=tmax){
	fit_data_dj.push_back(idx, rat_cor_dj);
	fit_data_j.push_back(idx, rat_cor);
	++idx;
      }
    }
  }

    //Set up the minimizer
  MarquardtLevenbergParameters<double> minparams;
  minparams.verbose = true;
  // if(cmdline.load_mlparams){
  //   parse(minparams, cmdline.mlparams_file);
  //   std::cout << "Loaded minimizer params: " << minparams << std::endl;
  // }else{
  //   minparams.verbose = true;
  // }
  FitConstant ff;
  simpleFitFuncWrapper<FitConstant> fitfunc(ff);

  simpleFitWrapper<jackknifeDistributionD> fitter(fitfunc, MinimizerType::MarquardtLevenberg, minparams);
  fitter.generateCovarianceMatrix(fit_data_dj, args.covariance_strategy);
  
  jackknifeDistribution<parameterVectorD> params(nsample,parameterVectorD(1,1.0));
  jackknifeDistributionD chisq(nsample,0.);
  jackknifeDistributionD chisq_per_dof(nsample,0.);
  int dof;
  assert(fitter.fit(params, chisq, chisq_per_dof, dof, fit_data_j));

  std::cout << "ZV: " << distributionStructPeek(params,0) << std::endl;
  std::cout << "chisq: " << chisq << std::endl;
  std::cout << "chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "dof: " << dof << std::endl;

  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(params, "params.hdf5"); 



  // const int nchannel = args.data.size();
  // std::vector<doubleJackknifeCorrelationFunctionD> channels_dj(nchannel);
  // std::vector<jackknifeCorrelationFunctionD> channels_j(nchannel);
  // for(int i=0;i<nchannel;i++){
  //   rawDataCorrelationFunctionD channel_raw;
  //   readData(channel_raw, args.data[i], args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
  //   bin(channel_raw, args.bin_size);    

  //   channels_dj[i] = doubleJackknifeCorrelationFunctionD(args.Lt, [&](const int t){
  // 	return typename doubleJackknifeCorrelationFunctionD::ElementType(t,  doubleJackknifeDistributionD(channel_raw.value(t)));
  //     });
  //   channels_j[i] = jackknifeCorrelationFunctionD(args.Lt, [&](const int t){
  // 	return typename jackknifeCorrelationFunctionD::ElementType(t,  jackknifeDistributionD(channel_raw.value(t)));
  //     });
  // }
  // doubleJackknifeCorrelationFunctionD data_dj;
  // applyCombination(data_dj,channels_dj,args.combination);
  // applyTimeDep(data_dj, args.outer_time_dep, args.Lt);

  // jackknifeCorrelationFunctionD data_j;
  // applyCombination(data_j,channels_j,args.combination);
  // applyTimeDep(data_j, args.outer_time_dep, args.Lt);

  // blockDoubleJackknifeCorrelationFunctionD data_bj;
  // assert(args.covariance_strategy != CovarianceStrategy::CorrelatedBlockHybrid &&
  // 	 args.covariance_strategy != CovarianceStrategy::CorrelatedBlock);

  // jackknifeDistribution<parameterVectorD> params;
  // jackknifeDistributionD chisq;
  // int dof;
  // fit(params, chisq, dof, data_j,data_dj,data_bj, args, cmdline);
  

  return 0;
}

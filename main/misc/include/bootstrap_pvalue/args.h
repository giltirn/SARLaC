#pragma once


//Note: cov_strat_params_file only required for those strategies that have tunable parameters
struct Args{

#define ARGS_MEM (int, nsample)(int, Lt)(int, ntest)(int, norig_ens)(DataGenStrategy, data_strat)(int, block_size)(FitFuncType, fitfunc)(CovMatStrategy, cov_strat)(std::string, cov_strat_params_file)(BootResampleTableType, bootstrap_strat)(std::vector<preAnalysisType>, preanalysis)(std::vector<std::string>, preanalysis_params_file)(MarquardtLevenbergParameters<double>, MLparams)

  //norig_ens : values larger than 1 will loop over multiple original ensembles and generate statistical errors for the bootstrap estimates
  //preanalysis_params_file : for preanalyses that don't require params files, any string is valid
  GENERATE_MEMBERS(ARGS_MEM);
  
Args(): nsample(200), Lt(30), ntest(5000), norig_ens(1), block_size(1), fitfunc(FitFuncType::FConstant), cov_strat(CovMatStrategy::Correlated), data_strat(DataGenStrategy::NormalUniform), cov_strat_params_file(""), bootstrap_strat(BootResampleTableType::NonOverlappingBlock), preanalysis(1,preAnalysisType::None), preanalysis_params_file(1,""){
    MLparams.verbose = true;
    MLparams.lambda_factor = 1.2;
    MLparams.dampening_matrix = MLdampeningMatrix::Unit;
    MLparams.max_iter = 50000;
  }
};
GENERATE_PARSER( Args, ARGS_MEM );


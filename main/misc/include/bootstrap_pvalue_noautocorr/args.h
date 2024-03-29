#pragma once


//Note: cov_strat_params_file only required for those strategies that have tunable parameters
struct Args{

#define ARGS_MEM (int, nsample)(int, Lt)(int, ntest)(int, norig_ens)(DataGenStrategy, data_strat)(int, block_size)(FitFuncType, fitfunc)(CovMatStrategy, cov_strat)(std::string, cov_strat_params_file)(std::vector<preAnalysisType>, preanalysis)(MarquardtLevenbergParameters<double>, MLparams)

  //norig_ens : values larger than 1 will loop over multiple original ensembles and generate statistical errors for the bootstrap estimates

  GENERATE_MEMBERS(ARGS_MEM);
  
Args(): nsample(200), Lt(30), ntest(5000), norig_ens(1), block_size(1), fitfunc(FitFuncType::FConstant), cov_strat(CovMatStrategy::Correlated), data_strat(DataGenStrategy::NormalUniform), cov_strat_params_file(""), preanalysis(1,preAnalysisType::None){
    MLparams.verbose = true;
    MLparams.lambda_factor = 1.2;
    MLparams.dampening_matrix = MLdampeningMatrix::Unit;
    MLparams.max_iter = 50000;
  }
};
GENERATE_PARSER( Args, ARGS_MEM );


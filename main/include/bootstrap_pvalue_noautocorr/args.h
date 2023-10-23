#pragma once


//Note: cov_strat_params_file only required for those strategies that have tunable parameters
struct Args{
#define ARGS_MEM (int, nsample)(int, Lt)(int, ntest)(DataGenStrategy, data_strat)(FitFuncType, fitfunc)(CovMatStrategy, cov_strat)(std::string, cov_strat_params_file)(preAnalysisType, preanalysis)(MarquardtLevenbergParameters<double>, MLparams)
  GENERATE_MEMBERS(ARGS_MEM);
  
  Args(): nsample(200), Lt(30), ntest(5000), fitfunc(FitFuncType::FConstant), cov_strat(CovMatStrategy::Correlated), data_strat(DataGenStrategy::NormalUniform), cov_strat_params_file(""), preanalysis(preAnalysisType::None){
    MLparams.verbose = true;
    MLparams.lambda_factor = 1.2;
    MLparams.dampening_matrix = MLdampeningMatrix::Unit;
    MLparams.max_iter = 50000;
  }
};
GENERATE_PARSER( Args, ARGS_MEM );


#ifndef _SIMPLE_FIT_WRAPPER_UTILS_H_
#define _SIMPLE_FIT_WRAPPER_UTILS_H_

//Convenience functions and types for loading data for frozen fits

#include<config.h>
#include<utils/macros.h>
#include<fit/simple_fit_wrapper/fitter.h>

CPSFIT_START_NAMESPACE

template<typename T>
generalContainer getMinimizerParamsT(const bool load_minimizer_params, const std::string &minimizer_params_file){
  generalContainer min_params;
  T mp; mp.verbose = true;
  if(load_minimizer_params){
    parse(mp, minimizer_params_file);
    std::cout << "Loaded minimizer params: " << mp << std::endl;
  }
  min_params = mp;
  return min_params;
}
generalContainer getMinimizerParams(const MinimizerType minimizer, const bool load_minimizer_params, const std::string &minimizer_params_file){
  switch(minimizer){
  case MinimizerType::MarquardtLevenberg:
    return getMinimizerParamsT<MarquardtLevenbergParameters<double> >(load_minimizer_params, minimizer_params_file);
  case MinimizerType::GSLtrs:
    return getMinimizerParamsT<GSLtrsMinimizerParams>(load_minimizer_params, minimizer_params_file);
  case MinimizerType::GSLmultimin:
    return getMinimizerParamsT<GSLmultidimMinimizerParams>(load_minimizer_params, minimizer_params_file);
  case MinimizerType::Minuit2:
#ifdef HAVE_MINUIT2
    return getMinimizerParamsT<Minuit2minimizerParams>(load_minimizer_params, minimizer_params_file);
#else
    error_exit(std::cout << "Library not compiled with Minuit2\n");
#endif
  default:
    assert(0);
  }
}

CPSFIT_END_NAMESPACE

#endif

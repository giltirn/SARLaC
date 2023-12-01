#ifndef _SARLAC_SIMPLE_FIT_WRAPPER_COMMON_H
#define _SARLAC_SIMPLE_FIT_WRAPPER_COMMON_H

#include<config.h>
#include<utils/macros.h>

#include<fit/simple_fit_wrapper/fitfunc_wrapper.h>
#include<minimizer/minimizer.h>
#include<minimizer/gsl_trs_minimizer.h>
#include<minimizer/gsl_multidim_minimizer.h>
#include<minimizer/minuit2_minimizer.h>
#include<fit/fitfunc/fitfunc_frozen.h>
#include<fit/fitfunc/fitfunc_bounded.h>
#include<tensors/numeric_square_matrix.h>
#include<data_series/correlationfunction.h>
#include<fit/cost_function/correlated_chisq.h>
#include<fit/cost_function/correlated_chisq_terms.h>

SARLAC_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(MinimizerType, (MarquardtLevenberg)(GSLtrs)(GSLmultimin)(Minuit2) );
GENERATE_ENUM_AND_PARSER(CostType, (Correlated)(Uncorrelated) );

struct simpleFitCommon{
  typedef genericFitFuncBase FitFunc;
  typedef BoundedFitFunc<genericFitFuncBase> FitFuncBounded;
  typedef FrozenFitFunc<FitFuncBounded> FitFuncFrozenBounded;
  typedef FitFunc::ParameterType ParameterType;

  struct Prior{
    double value;
    double weight;
    int param_idx;
    Prior(const double value, const double weight, const int param_idx): value(value), weight(weight), param_idx(param_idx){}
  };

#define TPDF(T) typedef simpleFitCommon::T T
#define INHERIT_COMMON_TYPEDEFS TPDF(FitFunc); TPDF(FitFuncBounded); TPDF(FitFuncFrozenBounded); TPDF(ParameterType); TPDF(Prior);

  template<typename CostFunc, typename FitFuncInternal>
  static void addPriors(CostFunc &cost, const FitFuncInternal &fitfunc_s, const std::vector<Prior> &priors){
    for(int p=0;p<priors.size();p++){
      int subset_pidx = fitfunc_s.getParamsSubsetIndex(priors[p].param_idx);
      if(subset_pidx != -1)
	cost.addPrior(priors[p].value, priors[p].weight, subset_pidx);
    }
  }

  template<typename FitFuncInternal>
  static double fitSampleML(bool &converged, ParameterType &params_s, int &dof,
			    const correlationFunction<generalContainer, double> &data_s,
			    const NumericSquareMatrix<double> &inv_corr_s,
			    const std::vector<double> &sigma_s,
			    const FitFuncInternal &fitfunc_s,
			    const generalContainer &min_params,
			    const std::vector<Prior> &priors,
			    std::pair<double, int> *chisq_dof_nopriors = NULL){ //for frozen fit funcs, the values of the frozen parameters are sample dependent
    typedef CorrelatedChisqCostFunction<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;
    CostFunc cost(fitfunc_s, data_s, sigma_s, inv_corr_s);
    addPriors(cost, fitfunc_s, priors);
    
    dof = cost.Ndof();
    
    typedef MarquardtLevenbergMinimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    double chisq = min.fit(params_s);

    converged = min.hasConverged();

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, inv_corr_s);
      chisq_dof_nopriors->first = cost_nopriors.cost(params_s);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }

    return chisq;
  }

  template<typename FitFuncInternal>
  static double fitSampleGSLtrs(bool &converged, ParameterType &params_s, int &dof,
				const correlationFunction<generalContainer, double> &data_s,
				const NumericSquareMatrix<double> &corr_s,
				const std::vector<double> &sigma_s,
				const FitFuncInternal &fitfunc_s,
				const generalContainer &min_params,
				const std::vector<Prior> &priors,
				std::pair<double, int> *chisq_dof_nopriors = NULL){ //for frozen fit funcs, the values of the frozen parameters are sample dependent
    std::vector<double> corr_evals;
    std::vector<NumericVector<double> > corr_evecs;
    symmetricMatrixEigensolve(corr_evecs, corr_evals, corr_s);

    typedef CorrelatedChisqCostFunctionTerms<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;

    CostFunc cost(fitfunc_s, data_s, sigma_s, corr_evals, corr_evecs);
    addPriors(cost, fitfunc_s, priors);
    
    dof = cost.Ndof();
    
    typedef GSLtrsMinimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    converged = min.hasConverged();
    
    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, corr_evals, corr_evecs);
      NumericVector<double> cost_terms = cost_nopriors.costVector(params_s);
      chisq_dof_nopriors->first = 0.5*dot(cost_terms,cost_terms);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }
    
    return chisq;
  }

  template<typename FitFuncInternal>
  static double fitSampleGSLmultimin(bool &converged, ParameterType &params_s, int &dof,
				     const correlationFunction<generalContainer, double> &data_s,
				     const NumericSquareMatrix<double> &inv_corr_s,
				     const std::vector<double> &sigma_s,
				     const FitFuncInternal &fitfunc_s,
				     const generalContainer &min_params,
				     const std::vector<Prior> &priors,			      
				     std::pair<double, int> *chisq_dof_nopriors = NULL){ //for frozen fit funcs, the values of the frozen parameters are sample dependent
    typedef CorrelatedChisqCostFunction<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;
    CostFunc cost(fitfunc_s, data_s, sigma_s, inv_corr_s);
    addPriors(cost, fitfunc_s, priors);
    
    dof = cost.Ndof();
    
    typedef GSLmultidimMinimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    converged = min.hasConverged();

    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, inv_corr_s);
      chisq_dof_nopriors->first = cost_nopriors.cost(params_s);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }

    return chisq;
  }

  template<typename FitFuncInternal>
  static double fitSampleMinuit2(bool &converged, ParameterType &params_s, int &dof,
				 const correlationFunction<generalContainer, double> &data_s,
				 const NumericSquareMatrix<double> &inv_corr_s,
				 const std::vector<double> &sigma_s,
				 const FitFuncInternal &fitfunc_s,
				 const generalContainer &min_params,
				 const std::vector<Prior> &priors,
				 std::pair<double, int> *chisq_dof_nopriors = NULL){ //for frozen fit funcs, the values of the frozen parameters are sample dependent
#ifndef HAVE_MINUIT2
    error_exit(std::cout << "Library not compiled with Minuit2" << std::endl);
#else
    typedef CorrelatedChisqCostFunction<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;
    CostFunc cost(fitfunc_s, data_s, sigma_s, inv_corr_s);
    addPriors(cost, fitfunc_s, priors);
    
    dof = cost.Ndof();
    
    typedef Minuit2minimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    converged = min.hasConverged();

    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, inv_corr_s);
      chisq_dof_nopriors->first = cost_nopriors.cost(params_s);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }

    return chisq;
#endif
  }

};

SARLAC_END_NAMESPACE

#endif

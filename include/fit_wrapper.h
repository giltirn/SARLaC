#ifndef _FIT_WRAPPER_H_
#define _FIT_WRAPPER_H_

//A convenience framework for performing combinations of correlated/uncorrelated, frozen/unfrozen fits with a common interface hiding much of the boilerplate
//It uses a policy-based approach allowing a significant amount of flexibility

#include<correlationfunction.h>

struct baseFitTypedefs{
  typedef correlationFunction<jackknifeDistributionD> jackknifeCorrelationFunction;
  typedef NumericMatrix<jackknifeDistributionD> jackknifeMatrix;
  typedef NumericVector<jackknifeDistributionD> jackknifeVector;

  typedef sampleSeries<const jackknifeCorrelationFunction> sampleSeriesType;
  typedef NumericMatrixSampleView<const jackknifeMatrix> sampleInvCorrType;
};
#define INHERIT_TYPEDEF(FROM,DEF) typedef typename FROM::DEF DEF

#define INHERIT_BASE_FIT_TYPEDEFS(FROM)			       \
  INHERIT_TYPEDEF(FROM,jackknifeCorrelationFunction); \
  INHERIT_TYPEDEF(FROM,jackknifeMatrix); \
  INHERIT_TYPEDEF(FROM,jackknifeVector); \
  INHERIT_TYPEDEF(FROM,sampleSeriesType); \
  INHERIT_TYPEDEF(FROM,sampleInvCorrType)


template<typename _fitFunc>
class standardFitFuncPolicy: public baseFitTypedefs{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs);

  typedef _fitFunc fitFunc;
  typedef jackknifeDistribution<typename fitFunc::ParameterType> jackknifeFitParameters;

  class fitFuncPolicyState{
    fitFunc const* ff;
  public:
    fitFuncPolicyState(fitFunc const* _ff): ff(_ff){}
    const fitFunc & getFitFunc() const{
      assert(ff != NULL);
      return *ff;
    }
    typename fitFunc::ParameterType getParamsSample(const jackknifeFitParameters &params, const int s){ return params.sample(s); }
    void setParamsSample(jackknifeFitParameters &params, const typename fitFunc::ParameterType &p_s, const int s){ params.sample(s) = p_s; }
  };
private:
  fitFunc const* ff;
public:
  standardFitFuncPolicy(): ff(NULL){}
  
  void importFitFunc(const fitFunc &fitfunc){ ff = &fitfunc; }
  
  inline fitFuncPolicyState generateFitFuncPolicyState(const int s){
    return fitFuncPolicyState(ff);
  }
};

template<typename _fitFunc>
class frozenFitFuncPolicy: public baseFitTypedefs{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs);

  typedef _fitFunc baseFitFunc;
  typedef FrozenFitFunc<_fitFunc> fitFunc;
  typedef jackknifeDistribution<typename baseFitFunc::ParameterType> jackknifeFitParameters; //going to intercept and transform to internal representation

  class fitFuncPolicyState{
    std::unique_ptr<fitFunc> ff_frz;
  public:
    fitFuncPolicyState(baseFitFunc const* ff, const std::vector<int> &freeze_params, const std::unique_ptr<jackknifeFitParameters> &freeze_vals, const int s){
      assert(ff != NULL);
      ff_frz.reset(new fitFunc(*ff));
      if(freeze_params.size() > 0){
	assert(freeze_vals);
	ff_frz->freeze(freeze_params, freeze_vals->sample(s));
      }
    }
    fitFuncPolicyState(fitFuncPolicyState &&r): ff_frz(std::move(r.ff_frz)){}
    
    const fitFunc & getFitFunc() const{
      return *ff_frz;
    }
    typename fitFunc::ParameterType getParamsSample(const jackknifeFitParameters &params, const int s){ return ff_frz->mapParamsSupersetToSubset(params.sample(s)); }
    void setParamsSample(jackknifeFitParameters &params, const typename fitFunc::ParameterType &p_s, const int s){ params.sample(s) = ff_frz->mapParamsSubsetToSuperset(p_s); }    
  };
private:
  baseFitFunc const* ff;
  std::unique_ptr<jackknifeFitParameters> freeze_vals;
  std::vector<int> freeze_params;
public:
  frozenFitFuncPolicy(): ff(NULL){}
  
  void importFitFunc(const baseFitFunc &fitfunc){ ff = &fitfunc; }

  void freeze(const std::vector<int> &_freeze_params, const jackknifeFitParameters &_freeze_vals){
    freeze_vals.reset(new jackknifeFitParameters(_freeze_vals));   freeze_params = _freeze_params;
  }
  
  inline fitFuncPolicyState generateFitFuncPolicyState(const int s){
    return fitFuncPolicyState(ff,freeze_params,freeze_vals,s);
  }
};




#define INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM)		       \
  INHERIT_BASE_FIT_TYPEDEFS(FROM);			       \
  INHERIT_TYPEDEF(FROM,fitFunc); \
  INHERIT_TYPEDEF(FROM,jackknifeFitParameters); \
  INHERIT_TYPEDEF(FROM,fitFuncPolicyState)



  

template<typename FitFuncPolicy>
class correlatedFitPolicy: public FitFuncPolicy{
public:
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FitFuncPolicy);
  typedef CorrelatedChisqCostFunction<fitFunc, sampleSeriesType, sampleInvCorrType> costFunctionType;
  
  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    sampleInvCorrType inv_corr_s;
    std::vector<double> sigma_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const jackknifeCorrelationFunction &data, const jackknifeMatrix &inv_corr, const jackknifeVector &sigma, const int s):
      data_s(data,s), inv_corr_s(inv_corr,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = sigma[d].sample(s);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s, inv_corr_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    const int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  jackknifeMatrix const* inv_corr;
  jackknifeVector const* sigma;
public:
  correlatedFitPolicy(): inv_corr(NULL), sigma(NULL){}
  
  void importCostFunctionParameters(const jackknifeMatrix &_inv_corr, 
				    const jackknifeVector &_sigma){
    inv_corr = &_inv_corr;
    sigma = &_sigma;
  }
protected:
  inline costFunctionState generateCostFunctionState(const jackknifeCorrelationFunction &data, const int s){
    assert(inv_corr != NULL && sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*inv_corr,*sigma,s);
  }
};

template<typename FitFuncPolicy>
class uncorrelatedFitPolicy: public FitFuncPolicy{
public:
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FitFuncPolicy);
  typedef UncorrelatedChisqCostFunction<fitFunc, sampleSeriesType> costFunctionType;
  
  //To ensure each loop iteration independent we can't store global state, instead we have thread local state
  class costFunctionState: public fitFuncPolicyState{
  private:
    sampleSeriesType data_s;
    std::vector<double> sigma_s;
    std::unique_ptr<costFunctionType> cost_func;

  public:
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const jackknifeCorrelationFunction &data, const jackknifeVector &sigma, const int s):
      data_s(data,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = sigma[d].sample(s);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    const int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  jackknifeVector const* sigma;
public:
  uncorrelatedFitPolicy(): sigma(NULL){}
  
  void importCostFunctionParameters(const jackknifeVector &_sigma){
    sigma = &_sigma;
  }
protected:
  inline costFunctionState generateCostFunctionState(const jackknifeCorrelationFunction &data, const int s){
    assert(sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*sigma,s);
  }
};



#define INHERIT_FIT_POLICY_TYPEDEFS(FROM) \
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM); \
  INHERIT_TYPEDEF(FROM,costFunctionType); \
  INHERIT_TYPEDEF(FROM,costFunctionState)


//Convenience wrapper for composing the fit policy
template<typename FitFunc, template<typename> class FitFuncPolicy, template<typename> class CostFunctionPolicy>
struct composeFitPolicy{
  typedef CostFunctionPolicy<FitFuncPolicy<FitFunc> > type;
};


template<typename FitPolicies>
struct fitter: public FitPolicies{
  INHERIT_FIT_POLICY_TYPEDEFS(FitPolicies);
  
  typedef MarquardtLevenbergMinimizer<costFunctionType> minimizerType;
  typedef typename minimizerType::AlgorithmParameterType minimizerParamsType;
  

  void fit(jackknifeFitParameters &params,
	   jackknifeDistributionD &chisq,
	   jackknifeDistributionD &chisq_per_dof,
	   const jackknifeCorrelationFunction &data){

    const int ndata_fit = data.size();
    assert(ndata_fit > 0);
    const int nsample = data.value(0).size();
    assert(params.size() == nsample);
        
    chisq.resize(nsample);
    chisq_per_dof.resize(nsample);
    
    std::cout << "Starting fit with guess " << params << std::endl;

    int nparams;
    minimizerParamsType min_params;
    min_params.verbose = true;
#pragma omp parallel for
    for(int s=0;s<nsample;s++){
      costFunctionState state = this->generateCostFunctionState(data,s);
      const costFunctionType &cost_func = state.getCostFunction();
      minimizerType minimizer(cost_func,min_params);

      typename fitFunc::ParameterType p_s = state.getParamsSample(params,s); //allows for arbitrary mapping from external to internal parameter type
      chisq.sample(s) = minimizer.fit(p_s);
      state.setParamsSample(params, p_s, s);
      
      chisq_per_dof.sample(s) = chisq.sample(s)/state.Ndof();
      assert(minimizer.hasConverged());
    }
  }
};



#endif

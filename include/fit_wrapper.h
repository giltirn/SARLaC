#ifndef _FIT_WRAPPER_H_
#define _FIT_WRAPPER_H_

//A convenience framework for performing combinations of correlated/uncorrelated, frozen/unfrozen fits with a common interface hiding much of the boilerplate
//It uses a policy-based approach allowing a significant amount of flexibility

#include<correlationfunction.h>
#include<common_defs.h>
#include<fitfunc_mapping.h>
#include<cost_function.h>
#include<minimizer.h>
#include<distribution_iterate.h>

#define INHERIT_TYPEDEF(FROM,DEF) typedef typename FROM::DEF DEF

//User can override these basic input types with custom types
template<typename GeneralizedCoordinate>
struct standardInputFitTypes{
  typedef jackknifeDistributionD DistributionType;
  typedef correlationFunction<GeneralizedCoordinate,DistributionType> CorrelationFunctionDistribution;
};

#define INHERIT_INPUT_FIT_TYPEDEFS(FROM)			       \
  INHERIT_TYPEDEF(FROM,DistributionType); \
  INHERIT_TYPEDEF(FROM,CorrelationFunctionDistribution)


template<typename BaseTypes>
struct baseFitTypedefs{
  INHERIT_INPUT_FIT_TYPEDEFS(BaseTypes);
  
  typedef NumericSquareMatrix<DistributionType> MatrixDistribution;
  typedef NumericVector<DistributionType> VectorDistribution;

  typedef sampleSeries<const CorrelationFunctionDistribution> sampleSeriesType;
  typedef NumericSquareMatrixSampleView<const MatrixDistribution> sampleInvCorrType;
};

#define INHERIT_BASE_FIT_TYPEDEFS(FROM)			       \
  INHERIT_INPUT_FIT_TYPEDEFS(FROM);			       \
  INHERIT_TYPEDEF(FROM,MatrixDistribution); \
  INHERIT_TYPEDEF(FROM,VectorDistribution); \
  INHERIT_TYPEDEF(FROM,sampleSeriesType); \
  INHERIT_TYPEDEF(FROM,sampleInvCorrType)


template<typename InputFitTypes, typename _fitFunc>
class standardFitFuncPolicy: public baseFitTypedefs<InputFitTypes>{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs<InputFitTypes>);
  typedef _fitFunc baseFitFunc;
  typedef _fitFunc fitFunc;
  typedef typename DistributionType::template rebase<typename fitFunc::ParameterType> FitParameterDistribution;

  class fitFuncPolicyState{
    fitFunc const* ff;
  public:
    fitFuncPolicyState(fitFunc const* _ff): ff(_ff){}
    const fitFunc & getFitFunc() const{
      assert(ff != NULL);
      return *ff;
    }
    typename fitFunc::ParameterType getParamsSample(const FitParameterDistribution &params, const int s){ return iterate<FitParameterDistribution>::at(s,params); }
    void setParamsSample(FitParameterDistribution &params, const typename fitFunc::ParameterType &p_s, const int s){ iterate<FitParameterDistribution>::at(s,params) = p_s; }
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

template<typename InputFitTypes, typename _fitFunc>
class frozenFitFuncPolicy: public baseFitTypedefs<InputFitTypes>{
public:
  INHERIT_BASE_FIT_TYPEDEFS(baseFitTypedefs<InputFitTypes>);

  typedef _fitFunc baseFitFunc;
  typedef FrozenFitFunc<_fitFunc> fitFunc;
  typedef jackknifeDistribution<typename baseFitFunc::ParameterType> FitParameterDistribution; //going to intercept and transform to internal representation

  class fitFuncPolicyState{
    std::unique_ptr<fitFunc> ff_frz;
  public:
    fitFuncPolicyState(baseFitFunc const* ff, const std::vector<int> &freeze_params, const std::unique_ptr<FitParameterDistribution> &freeze_vals, const int s){
      assert(ff != NULL);
      ff_frz.reset(new fitFunc(*ff));
      if(freeze_params.size() > 0){
	assert(freeze_vals);
	ff_frz->freeze(freeze_params, iterate<FitParameterDistribution>::at(s,*freeze_vals));
      }
    }
    fitFuncPolicyState(fitFuncPolicyState &&r): ff_frz(std::move(r.ff_frz)){}
    
    const fitFunc & getFitFunc() const{
      return *ff_frz;
    }
    typename fitFunc::ParameterType getParamsSample(const FitParameterDistribution &params, const int s){ return ff_frz->mapParamsSupersetToSubset(iterate<FitParameterDistribution>::at(s,params)); }
    void setParamsSample(FitParameterDistribution &params, const typename fitFunc::ParameterType &p_s, const int s){ iterate<FitParameterDistribution>::at(s,params) = ff_frz->mapParamsSubsetToSuperset(p_s); }    
  };
private:
  baseFitFunc const* ff;
  std::unique_ptr<FitParameterDistribution> freeze_vals;
  std::vector<int> freeze_params;
public:
  frozenFitFuncPolicy(): ff(NULL){}
  
  void importFitFunc(const baseFitFunc &fitfunc){ ff = &fitfunc; }

  void freeze(const std::vector<int> &_freeze_params, const FitParameterDistribution &_freeze_vals){
    freeze_vals.reset(new FitParameterDistribution(_freeze_vals));   freeze_params = _freeze_params;
  }
  
  inline fitFuncPolicyState generateFitFuncPolicyState(const int s){
    return fitFuncPolicyState(ff,freeze_params,freeze_vals,s);
  }
};




#define INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM)		       \
  INHERIT_BASE_FIT_TYPEDEFS(FROM);			       \
  INHERIT_TYPEDEF(FROM,fitFunc); \
  INHERIT_TYPEDEF(FROM,baseFitFunc); \
  INHERIT_TYPEDEF(FROM,FitParameterDistribution); \
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
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const MatrixDistribution &inv_corr, const VectorDistribution &sigma, const int s):
      data_s(data,s), inv_corr_s(inv_corr,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = iterate<DistributionType>::at(s,sigma[d]);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s, inv_corr_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    const int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  MatrixDistribution const* inv_corr;
  VectorDistribution const* sigma;
public:
  correlatedFitPolicy(): inv_corr(NULL), sigma(NULL){}
  
  void importCostFunctionParameters(const MatrixDistribution &_inv_corr, 
				    const VectorDistribution &_sigma){
    inv_corr = &_inv_corr;
    sigma = &_sigma;
  }
protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
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
    costFunctionState(fitFuncPolicyState &&_fitfuncstate, const CorrelationFunctionDistribution &data, const VectorDistribution &sigma, const int s):
      data_s(data,s), sigma_s(data.size()), fitFuncPolicyState(std::forward<fitFuncPolicyState>(_fitfuncstate)){
      for(int d=0;d<data.size();d++)
	sigma_s[d] = iterate<DistributionType>::at(s,sigma[d]);
      cost_func.reset(new costFunctionType(this->getFitFunc(),data_s, sigma_s));
    }
    const costFunctionType &getCostFunction() const{ return *cost_func; }

    const int Ndof() const{ return data_s.size() - this->getFitFunc().Nparams(); }
  };

private:
  VectorDistribution const* sigma;
public:
  uncorrelatedFitPolicy(): sigma(NULL){}
  
  void importCostFunctionParameters(const VectorDistribution &_sigma){
    sigma = &_sigma;
  }
protected:
  inline costFunctionState generateCostFunctionState(const CorrelationFunctionDistribution &data, const int s){
    assert(sigma != NULL);
    return costFunctionState(this->generateFitFuncPolicyState(s),data,*sigma,s);
  }
};



#define INHERIT_FIT_POLICY_TYPEDEFS(FROM) \
  INHERIT_FITFUNC_POLICY_TYPEDEFS(FROM); \
  INHERIT_TYPEDEF(FROM,costFunctionType); \
  INHERIT_TYPEDEF(FROM,costFunctionState)


//Convenience wrapper for composing the fit policy
template<typename FitFunc, template<typename,typename> class FitFuncPolicy, template<typename> class CostFunctionPolicy, typename InputFitTypes = standardInputFitTypes<typename FitFunc::GeneralizedCoordinate> >
struct composeFitPolicy{
  typedef CostFunctionPolicy<FitFuncPolicy<InputFitTypes,FitFunc> > type;
};


template<typename FitPolicies>
class fitter: public FitPolicies{
public:
  INHERIT_FIT_POLICY_TYPEDEFS(FitPolicies);
  
  typedef MarquardtLevenbergMinimizer<costFunctionType> minimizerType;
  typedef typename minimizerType::AlgorithmParameterType minimizerParamsType;
private:
  minimizerParamsType min_params;
public:
  fitter(){
    min_params.verbose = true;
  }
  fitter(const minimizerParamsType &_min_params): min_params(_min_params){}
  
  void fit(FitParameterDistribution &params,
	   DistributionType &chisq,
	   DistributionType &chisq_per_dof,
	   const CorrelationFunctionDistribution &data){

    const int ndata_fit = data.size();
    assert(ndata_fit > 0);
    const int nsample = data.value(0).size();
    assert(params.size() == nsample);
        
    chisq.resize(nsample);
    chisq_per_dof.resize(nsample);
    
    std::cout << "Starting fit with guess " << params << std::endl;

    int nosample = iterate<DistributionType>::size(data.value(0)); //some distributions have extra 'samples', eg jackknifeCdistribution has a sample representing the central value
    
    int nparams;
#pragma omp parallel for
    for(int s=0;s<nosample;s++){
      costFunctionState state = this->generateCostFunctionState(data,s);
      const costFunctionType &cost_func = state.getCostFunction();
      minimizerType minimizer(cost_func,min_params);

      typename fitFunc::ParameterType p_s = state.getParamsSample(params,s); //allows for arbitrary mapping from external to internal parameter type
      iterate<DistributionType>::at(s,chisq) = minimizer.fit(p_s);
      state.setParamsSample(params, p_s, s);
      
      iterate<DistributionType>::at(s,chisq_per_dof) = iterate<DistributionType>::at(s,chisq)/state.Ndof();
      assert(minimizer.hasConverged());
    }
  }
};



//Convenience classes to generate correlation matrix / sigma for correlated and uncorrelated fit policy using double-jackknife data
template<template<typename> class corrUncorrFitPolicy, typename FitPolicies>
struct importCostFunctionParameters{};

template<typename FitPolicies>
struct importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistributionD>::value, "Currently only support jackknifeDistributionD");
  
  NumericVector<jackknifeDistributionD> sigma;

  template<typename GeneralizedCoord>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, doubleJackknifeDistributionD> &data): sigma(data.size()){
    for(int d=0;d<data.size();d++)
      sigma(d) = jackknifeDistributionD(sqrt(doubleJackknifeDistributionD::covariance(data.value(d) , data.value(d) ) ) );
    
    fitter.importCostFunctionParameters(sigma);
  }
};
template<typename FitPolicies>
struct importCostFunctionParameters<correlatedFitPolicy,FitPolicies>{
  static_assert(std::is_same<typename FitPolicies::DistributionType,jackknifeDistributionD>::value, "Currently only support jackknifeDistributionD");
  
  NumericSquareMatrix<jackknifeDistributionD> corr;
  NumericSquareMatrix<jackknifeDistributionD> inv_corr;
  NumericVector<jackknifeDistributionD> sigma;

  template<typename GeneralizedCoord>
  importCostFunctionParameters(fitter<FitPolicies> &fitter,
			       const correlationFunction<GeneralizedCoord, doubleJackknifeDistributionD> &data): sigma(data.size()){

    const int nsample = data.value(0).size();
    const int ndata = data.size();    
    NumericSquareMatrix<jackknifeDistributionD> cov(ndata);
    for(int i=0;i<ndata;i++){
      cov(i,i) = doubleJackknifeDistributionD::covariance(data.value(i), data.value(i));
      sigma(i) = sqrt(cov(i,i));

      for(int j=i+1;j<ndata;j++)
	cov(i,j) = cov(j,i) = doubleJackknifeDistributionD::covariance(data.value(i),data.value(j));
    }

    corr =  NumericSquareMatrix<jackknifeDistributionD>(ndata);
    
    for(int i=0;i<ndata;i++){
      corr(i,i) = jackknifeDistributionD(nsample,1.);
      
      for(int j=i+1;j<ndata;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
    }
    
    inv_corr.resize(ndata, jackknifeDistributionD(nsample));
    svd_inverse(inv_corr, corr);

    //Test the quality of the inverse
    NumericSquareMatrix<jackknifeDistributionD> test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistributionD(nsample,1.0);    
    std::cout << "|CorrMat * CorrMat^{-1} - 1|^2 = " << mod2(test) << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
  }
};


#endif

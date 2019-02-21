#ifndef _CPSFIT_COST_FUNCTION_CORRELATED_CHISQ_H_
#define _CPSFIT_COST_FUNCTION_CORRELATED_CHISQ_H_

//The chi^2 cost function with a non-diagonal covariance matrix

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix.h>
#include<fit/cost_function/invert_policy.h>

CPSFIT_START_NAMESPACE

//This version uses the inverse of the correlation matrix and sigma
template<typename FitFunction, typename DataContainer,
	 typename InvCorrMatrixType = NumericSquareMatrix<double>,
	 typename _CostType = double,
	 typename _CostDerivativeType = NumericVector<_CostType>,
	 typename _CostSecondDerivativeMatrixType = NumericSquareMatrix<_CostType>,
	 typename _CostSecondDerivativeInverseMatrixType = NumericSquareMatrix<_CostType>,
	 typename MatrixInvertPolicy = CostFunctionSVDinvert<_CostType>,
	 typename std::enable_if<std::is_floating_point<typename FitFunction::ValueType>::value, int>::type = 0>
class CorrelatedChisqCostFunction: public MatrixInvertPolicy{
  static_assert(std::is_same<typename FitFunction::ValueType, typename DataContainer::DataType>::value, "DataContainer and FitFunction must have same value type");
  static_assert(std::is_same<typename FitFunction::GeneralizedCoordinate, typename DataContainer::GeneralizedCoordinate>::value, "DataContainer and FitFunction must have same coordinate");

public:
  typedef _CostType CostType;
  typedef typename FitFunction::ValueType ValueType;
  typedef typename FitFunction::ParameterType ParameterType;
  typedef typename FitFunction::ValueDerivativeType ValueDerivativeType; //derivative wrt parameters
  typedef typename FitFunction::GeneralizedCoordinate CoordinateType;

  //Types containing the derivative of the cost with respect to the parameters. Use defaults or specify manually using types conforming to spec
  typedef _CostDerivativeType CostDerivativeType;
  typedef _CostSecondDerivativeMatrixType CostSecondDerivativeMatrixType;
  typedef _CostSecondDerivativeInverseMatrixType CostSecondDerivativeInverseMatrixType;

  struct Prior{
    ValueType value;
    ValueType weight;
    int param_idx;
    
    Prior(ValueType value, ValueType weight, int param_idx): value(value), weight(weight), param_idx(param_idx){}

    inline CostType dcost(const ParameterType &params) const{ 
      CostType d = params(param_idx) - value;
      return d*d/weight/weight;
    }
    //Applies change directly to input derivs and second derivs
    inline void dderivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{ 
      derivs(param_idx) = derivs(param_idx) + 2*(params(param_idx) - value)/weight/weight;
      second_derivs(param_idx, param_idx) = second_derivs(param_idx, param_idx) + 2./weight/weight;
    }
  };
private:
  
  const FitFunction &fitfunc;
  const DataContainer &data;
  const std::vector<double> &sigma; //diagonal elements of covariance matrix
  const InvCorrMatrixType &inv_corr; //inverse correlation matrix
  std::vector<Prior> priors;

public:


  CorrelatedChisqCostFunction(const FitFunction &ff, const DataContainer &dd, const std::vector<double> &_sigma, const InvCorrMatrixType &_inv_corr):
    fitfunc(ff), data(dd), sigma(_sigma), inv_corr(_inv_corr){
    assert(inv_corr.size() == data.size());
    assert(sigma.size() == data.size());
  }
  
  //Add a prior for constrained curve fitting (https://arxiv.org/pdf/hep-lat/0110175.pdf). This adds a term to chi^2 for the parameter p with index 'param_idx'
  //of the form d\Chi^2 = ( p - value )^2/weight^2
  void addPrior(ValueType value, ValueType weight, int param_idx){
    priors.push_back(Prior(value,weight,param_idx));
  }
  void addPrior(const Prior &prior){
    priors.push_back(prior);
  }


  CostType cost(const ParameterType &params) const{
    CostType chisq = 0.;
    const int ndata = data.size();

    std::vector<ValueType> dfw(ndata);

    for(int a=0;a<ndata;a++){
      ValueType yfit_a = fitfunc.value(data.coord(a), params);
      dfw[a] = ( data.value(a) - yfit_a ) / sigma[a];

      chisq += dfw[a] * dfw[a] * inv_corr(a,a);

      for(int b=0;b<a;b++)
	chisq += 2*dfw[a] * dfw[b] * inv_corr(a,b); //matrix is symmetric
    } 
    
    //Contribution from priors (if any)
    for(int p=0; p<priors.size(); p++) chisq += priors[p].dcost(params);
     
    return chisq;
  }
  
  void derivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{
    if(params.size() != fitfunc.Nparams()) error_exit(std::cout << "CorrelatedChisqCostFunction::derivatives params size " << params.size() << " doesn't match the number of fit parameters " << fitfunc.Nparams() << std::endl);

    const int nparams = params.size();
    const int ndata = data.size();
    
    derivs.resize(nparams);
    second_derivs.resize(nparams);
    
    derivs.zero();
    second_derivs.zero();

    std::vector<ValueType> dfw(ndata);
    std::vector<std::unique_ptr<ValueDerivativeType> > yderivs(ndata); //storage but without assuming default constructible

    for(int a=0;a<ndata;a++){
      ValueType yfit_a = fitfunc.value(data.coord(a), params);
      yderivs[a] = std::unique_ptr<ValueDerivativeType>(new ValueDerivativeType(fitfunc.parameterDerivatives(data.coord(a), params) ) );
      dfw[a] = ( data.value(a) - yfit_a ) / sigma[a];

      const ValueDerivativeType & yderivs_a = *yderivs[a];
      
      //b==a
      for(int i=0;i<nparams;i++){
	derivs(i) +=  -2.0 * inv_corr(a,a) * dfw[a] * yderivs_a(i)/sigma[a];
	
	for(int j=0;j<=i;j++){
	  ValueType change = 2.0 * inv_corr(a,a) * yderivs_a(i)/sigma[a] * yderivs_a(j)/sigma[a];
	  second_derivs(i,j) += change;
	  if(i!=j) second_derivs(j,i) += change;
	}
      }
      //b<a
      for(int b=0;b<a;b++){	
	if(inv_corr(a,b)==0.0) continue;

	const ValueDerivativeType & yderivs_b = *yderivs[b];
	
	for(int i=0;i<nparams;i++){
	  derivs(i) +=  -2.0*inv_corr(a,b)*( yderivs_a(i)/sigma[a]*dfw[b] + dfw[a]* yderivs_b(i)/sigma[b] );
	
	  for(int j=0;j<=i;j++){
	    ValueType change = 2.0*inv_corr(a,b)*( yderivs_a(i)/sigma[a]* yderivs_b(j)/sigma[b] + yderivs_a(j)/sigma[a]*yderivs_b(i)/sigma[b] );
	    second_derivs(i,j) += change;
	    if(i!=j) second_derivs(j,i) += change;
	  }
	  
	}
      }
    }
    
    //Contribution from priors (if any)
    for(int p=0; p<priors.size(); p++) priors[p].dderivatives(derivs, second_derivs, params);
   
  }

  inline int Ndof() const{
    return data.size() + priors.size() - fitfunc.Nparams();
  }
  inline int Nparams() const{ return fitfunc.Nparams(); }    
};

CPSFIT_END_NAMESPACE
#endif

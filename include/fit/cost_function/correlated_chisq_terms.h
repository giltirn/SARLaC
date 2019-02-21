#ifndef _CPSFIT_COST_FUNCTION_CORRELATED_CHISQ_TERMS_H_
#define _CPSFIT_COST_FUNCTION_CORRELATED_CHISQ_TERMS_H_

//The chi^2 cost function with a non-diagonal covariance matrix
//GSL minimizer (and possibly others) require that we can express chi^2 = 1/2 \sum_i f_i^2
//and ask us to provide the elements of f and the Jacobian   df_i/dp_j  where p_j is a parameter

//Currently  chi^2 = \sum_i,j   d_i (C^-1)_ij d_j
//Thus we need to do an eigenvector decomposition of the inverse correlation matrix
//(C^-1)_ij =  \sum_a l^a  v^a_i v^a_j    where l^a is the eigenvalue of the inverse matrix
//then chi^2 =  \sum_a   (\sum_i d_i v^a_i  sqrt(l^a) )^2
//f_i = sqrt(2)\sum_i d_i v^a_i  sqrt(l^a)
//df_i / dp_j = sqrt(2)\sum_i  d d_i/dp_j  v^a_i sqrt(l^a)
//where   d d_i/dp_j = - (dF_i/dp_j) /sigma_i     with F the fit function



#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_tensor.h>
#include<fit/cost_function/invert_policy.h>

CPSFIT_START_NAMESPACE

//This version uses the inverse of the correlation matrix and sigma
template<typename FitFunction, typename DataContainer,
	 typename _CostType = double,
	 typename _CostTermsType = NumericVector<_CostType>,
	 typename _CostJacobianMatrixType = NumericTensor<_CostType,2>,
	 typename std::enable_if<std::is_floating_point<typename FitFunction::ValueType>::value, int>::type = 0>
class CorrelatedChisqCostFunctionTerms{
  static_assert(std::is_same<typename FitFunction::ValueType, typename DataContainer::DataType>::value, "DataContainer and FitFunction must have same value type");
  static_assert(std::is_same<typename FitFunction::GeneralizedCoordinate, typename DataContainer::GeneralizedCoordinate>::value, "DataContainer and FitFunction must have same coordinate");

public:
  typedef _CostType CostType;
  typedef typename FitFunction::ValueType ValueType;
  typedef typename FitFunction::ParameterType ParameterType;
  typedef typename FitFunction::ValueDerivativeType ValueDerivativeType; //derivative wrt parameters
  typedef typename FitFunction::GeneralizedCoordinate CoordinateType;

  typedef _CostTermsType CostTermsType;
  typedef _CostJacobianMatrixType CostJacobianMatrixType;

  struct Prior{
    ValueType value;
    ValueType weight;
    int param_idx;
    
    Prior(ValueType value, ValueType weight, int param_idx): value(value), weight(weight), param_idx(param_idx){}
  };

private:
  
  const FitFunction &fitfunc;
  const DataContainer &data;
  const std::vector<double> &sigma; //diagonal elements of covariance matrix

  std::vector<Prior> priors;

  //For GSL minimizer we must eigenvalue decompose correlation matrix
  std::vector<double> corr_evals;
  std::vector<NumericVector<double> > corr_evecs;

public:

  //Provide the eigenvectors and eigenvalues of the correlation matrix
  CorrelatedChisqCostFunctionTerms(const FitFunction &ff, const DataContainer &dd, const std::vector<double> &_sigma, 
				   const std::vector<double> &corr_evals, const std::vector<NumericVector<double> > &corr_evecs):
    fitfunc(ff), data(dd), sigma(_sigma), corr_evecs(corr_evecs), corr_evals(corr_evals){
    assert(corr_evals.size() == data.size());
    assert(corr_evecs.size() == data.size());
    assert(sigma.size() == data.size());
  }

  //Add a prior for constrained curve fitting (https://arxiv.org/pdf/hep-lat/0110175.pdf). This adds a term to chi^2 for the parameter p with index 'param_idx'
  //of the form d\Chi^2 = ( p - value )^2/weight^2
  void addPrior(ValueType value, ValueType weight, int param_idx){
    priors.push_back(Prior(value,weight,param_idx));
  }

  
  CostTermsType costVector(const ParameterType &params) const{
    const int nparams = params.size();
    const int ndata = data.size();
    const int nprior = priors.size();

    std::vector<CostType> delta(ndata);
    for(int i=0;i<ndata;i++){
      CostType yfit = fitfunc.value(data.coord(i), params);
      delta[i] = ( data.value(i) - yfit ) / sigma[i];
    }
    
    CostTermsType out(ndata + nprior);
    for(int a=0;a<ndata;a++){
      out(a) = 0.;
      for(int i=0;i<ndata;i++){
	out(a) = out(a) + sqrt(2.)*delta[i]*corr_evecs[a](i)/sqrt(corr_evals[a]);
      }
    }

    for(int i=ndata; i<ndata+nprior;i++){
      const Prior &pr = priors[i-ndata];
      out(i) = sqrt(2.)*(params(pr.param_idx) - pr.value)/pr.weight;
    }

    return out;
  }
  CostJacobianMatrixType costJacobianMatrix(const ParameterType &params) const{
    const int nparams = params.size();
    const int ndata = data.size();
    const int nprior = priors.size();
    
    std::vector<CostType> delta(ndata);
    std::vector<std::vector<CostType> > ddelta(ndata, std::vector<CostType>(nparams));

    for(int i=0;i<ndata;i++){
      CostType yfit_i = fitfunc.value(data.coord(i), params);
      delta[i] = ( data.value(i) - yfit_i ) / sigma[i];

      ValueDerivativeType yderivs_i = fitfunc.parameterDerivatives(data.coord(i), params);      
      for(int j=0;j<nparams;j++)
	ddelta[i][j] = -yderivs_i(j)/sigma[i];
    }

    CostJacobianMatrixType out({ndata+nprior,nparams});

    //df_i / dp_j = sqrt(2)\sum_i  d d_i/dp_j  v^a_i sqrt(l^a)
    for(int a=0;a<ndata;a++){
      for(int j=0;j<nparams;j++){
	auto & o = out(a,j);
	o = 0.;
	for(int i=0;i<ndata;i++)
	  o = o + sqrt(2.) * ddelta[i][j] * corr_evecs[a](i) / sqrt(corr_evals[a]);
      }
    }

    for(int i=ndata; i<ndata+nprior;i++){
      const Prior &pr = priors[i-ndata];
      for(int j=0;j<nparams;j++)
	out(i,j) = j == pr.param_idx ? sqrt(2.)/pr.weight : 0.;
    }

    return out;
  }
  inline int Nterms() const{ return data.size() + priors.size(); }

  inline int Ndof() const{
    return data.size() + priors.size() - fitfunc.Nparams();
  }
  inline int Nparams() const{ return fitfunc.Nparams(); }    

};

CPSFIT_END_NAMESPACE
#endif

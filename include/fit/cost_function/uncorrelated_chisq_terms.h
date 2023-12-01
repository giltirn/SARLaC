#ifndef _SARLAC_COST_FUNCTION_UNCORRELATED_CHISQ_TERMS_H_
#define _SARLAC_COST_FUNCTION_UNCORRELATED_CHISQ_TERMS_H_

//The chi^2 cost function with a diagonal covariance matrix
//GSL minimizer (and possibly others) require that we can express chi^2 = 1/2 \sum_i f_i^2
//and ask us to provide the elements of f and the Jacobian   df_i/dp_j  where p_j is a parameter

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_tensor.h>
#include<fit/cost_function/invert_policy.h>

SARLAC_START_NAMESPACE

template<typename FitFunction, typename DataContainer, typename _CostType = double,
	 typename _CostTermsType = NumericVector<_CostType>,
	 typename _CostJacobianMatrixType = NumericTensor<_CostType,2>,
	 typename std::enable_if<std::is_floating_point<typename FitFunction::ValueType>::value, int>::type = 0>
class UncorrelatedChisqCostFunctionTerms{
  static_assert(std::is_same<typename FitFunction::ValueType, typename DataContainer::DataType>::value, "DataContainer and FitFunction must have same value type");
  static_assert(std::is_same<typename FitFunction::GeneralizedCoordinate, typename DataContainer::GeneralizedCoordinate>::value, "DataContainer and FitFunction must have same coordinate");
  const FitFunction &fitfunc;
  const DataContainer &data;
  const std::vector<double> &sigma;
public:
  typedef _CostType CostType;
  typedef typename FitFunction::ValueType ValueType;
  typedef typename FitFunction::ParameterType ParameterType;
  typedef typename FitFunction::ValueDerivativeType ValueDerivativeType; //derivative wrt parameters
  typedef typename FitFunction::GeneralizedCoordinate CoordinateType;

  typedef _CostTermsType CostTermsType;
  typedef _CostJacobianMatrixType CostJacobianMatrixType;

  UncorrelatedChisqCostFunctionTerms(const FitFunction &ff, const DataContainer &dd, const std::vector<double> &_sigma): fitfunc(ff), data(dd), sigma(_sigma){ assert(sigma.size() == data.size()); }
  
  //GSL minimizer (and possibly others) require that we can express chi^2 = 1/2 \sum_i f_i^2
  //and ask us to provide the elements of f and the Jacobian   df_i/dp_j  where p_j is a parameter
  CostTermsType costVector(const ParameterType &params) const{
    const int nparams = params.size();
    const int ndata = data.size();
    CostTermsType out(ndata);
    for(int i=0;i<ndata;i++){
      ValueType yfit = fitfunc.value(data.coord(i), params);
      out(i) = sqrt(2.) * ( data.value(i) - yfit ) / sigma[i];
    }
    return out;
  }
  CostJacobianMatrixType costJacobianMatrix(const ParameterType &params) const{
    const int nparams = params.size();
    const int ndata = data.size();
    CostJacobianMatrixType out({ndata,nparams});
    for(int i=0;i<ndata;i++){
      ValueDerivativeType yderivs_i = fitfunc.parameterDerivatives(data.coord(i), params);
      for(int j=0;j<nparams;j++)
	out(i,j) = -sqrt(2.0) * yderivs_i(j)/sigma[i];
    }
    return out;
  }
  inline int Nterms() const{ return data.size(); } //number of terms in f

  inline int Nparams() const{ return fitfunc.Nparams(); }

  inline int Ndof() const{
    return data.size() - fitfunc.Nparams();
  }
};

SARLAC_END_NAMESPACE
#endif

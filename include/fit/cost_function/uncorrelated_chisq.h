#ifndef _CPSFIT_COST_FUNCTION_UNCORRELATED_CHISQ_H_
#define _CPSFIT_COST_FUNCTION_UNCORRELATED_CHISQ_H_

//The chi^2 cost function with a diagonal covariance matrix

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix.h>
#include<fit/cost_function/invert_policy.h>

CPSFIT_START_NAMESPACE

template<typename FitFunction, typename DataContainer, typename _CostType = double,
	 typename _CostDerivativeType = NumericVector<_CostType>,
	 typename _CostSecondDerivativeMatrixType = NumericSquareMatrix<_CostType>,
	 typename _CostSecondDerivativeInverseMatrixType = NumericSquareMatrix<_CostType>,
	 typename MatrixInvertPolicy = CostFunctionSVDinvert<_CostType>,
	 typename std::enable_if<std::is_floating_point<typename FitFunction::ValueType>::value, int>::type = 0>
class UncorrelatedChisqCostFunction: public MatrixInvertPolicy{
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

  //Types containing the derivative of the cost with respect to the parameters. Use defaults or specify manually using types conforming to spec
  typedef _CostDerivativeType CostDerivativeType;
  typedef _CostSecondDerivativeMatrixType CostSecondDerivativeMatrixType;
  typedef _CostSecondDerivativeInverseMatrixType CostSecondDerivativeInverseMatrixType;

  UncorrelatedChisqCostFunction(const FitFunction &ff, const DataContainer &dd, const std::vector<double> &_sigma): fitfunc(ff), data(dd), sigma(_sigma){ assert(sigma.size() == data.size()); }
  
  CostType cost(const ParameterType &params) const{
    CostType chisq = 0.;
    for(int i=0;i<data.size();i++){
      ValueType yfit = fitfunc.value(data.coord(i), params);
      ValueType dfw = ( data.value(i) - yfit ) / sigma[i];
      chisq += dfw * dfw;
    }
    return chisq;
  }
  
  void derivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{
    if(params.size() != fitfunc.Nparams()) error_exit(std::cout << "Error: Expected " << printType<ParameterType>() << " of size " << fitfunc.Nparams() << ", but size is " << params.size() << std::endl);
    const int nparams = params.size();
    const int ndata = data.size();
    
    derivs.resize(nparams);
    second_derivs.resize(nparams);
    
    derivs.zero();
    second_derivs.zero();
    
    for(int d=0;d<ndata;d++){
      ValueType yfit = fitfunc.value(data.coord(d), params);
      ValueDerivativeType yderivs = fitfunc.parameterDerivatives(data.coord(d), params);
      
      ValueType dfw = ( data.value(d) - yfit ) / sigma[d];
            
      for(int i=0;i<nparams;i++){
	derivs(i) += -2.0*dfw*yderivs(i)/sigma[d];

	for(int j=0;j<=i;j++){
	  ValueType change = 2.0 * yderivs(i)/sigma[d] * yderivs(j)/sigma[d];  //function treated as having no second derivative
	  second_derivs(i,j) += change;
	  if(i!=j) second_derivs(j,i) += change;
	}
      }
    }
  }

  inline int Ndof() const{
    return data.size() - fitfunc.Nparams();
  }
    

};

CPSFIT_END_NAMESPACE
#endif

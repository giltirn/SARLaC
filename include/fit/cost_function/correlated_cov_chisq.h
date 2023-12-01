#ifndef _SARLAC_COST_FUNCTION_CORRELATED_COV_CHISQ_H_
#define _SARLAC_COST_FUNCTION_CORRELATED_COV_CHISQ_H_

//The chi^2 cost function with a non-diagonal covariance matrix

#include<config.h>
#include<utils/macros.h>
#include<tensors/numeric_vector.h>
#include<tensors/numeric_square_matrix.h>
#include<fit/cost_function/invert_policy.h>

SARLAC_START_NAMESPACE

//This version uses the inverse of the covariance matrix
template<typename FitFunction, typename DataContainer,
	 typename InvCovMatrixType = NumericSquareMatrix<double>,
	 typename _CostType = double,
	 typename _CostDerivativeType = NumericVector<_CostType>,
	 typename _CostSecondDerivativeMatrixType = NumericSquareMatrix<_CostType>,
	 typename _CostSecondDerivativeInverseMatrixType = NumericSquareMatrix<_CostType>,
	 typename MatrixInvertPolicy = CostFunctionSVDinvert<_CostType>,
	 typename std::enable_if<std::is_floating_point<typename FitFunction::ValueType>::value, int>::type = 0>
class CorrelatedCovChisqCostFunction: public MatrixInvertPolicy{
  static_assert(std::is_same<typename FitFunction::ValueType, typename DataContainer::DataType>::value, "DataContainer and FitFunction must have same value type");
  static_assert(std::is_same<typename FitFunction::GeneralizedCoordinate, typename DataContainer::GeneralizedCoordinate>::value, "DataContainer and FitFunction must have same coordinate");
  const FitFunction &fitfunc;
  const DataContainer &data;
  const InvCovMatrixType &inv_cov; //inverse correlation matrix
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

  CorrelatedCovChisqCostFunction(const FitFunction &ff, const DataContainer &dd, const InvCovMatrixType &_inv_cov):
    fitfunc(ff), data(dd), inv_cov(_inv_cov){
    assert(inv_cov.size() == data.size());
  }
  
  CostType cost(const ParameterType &params) const{
    CostType chisq = 0.;
    const int ndata = data.size();

    std::vector<ValueType> dfw(ndata);

    for(int a=0;a<ndata;a++){
      ValueType yfit_a = fitfunc.value(data.coord(a), params);
      dfw[a] = data.value(a) - yfit_a;

      chisq += dfw[a] * dfw[a] * inv_cov(a,a);

      for(int b=0;b<a;b++)
	chisq += 2*dfw[a] * dfw[b] * inv_cov(a,b); //matrix is symmetric
    }      
    return chisq;
  }
  
  void derivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{
    assert(params.size() == fitfunc.Nparams());
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
      dfw[a] = data.value(a) - yfit_a;

      const ValueDerivativeType & yderivs_a = *yderivs[a];
      
      //b==a
      for(int i=0;i<nparams;i++){
	derivs(i) +=  -2.0 * inv_cov(a,a) * dfw[a] * yderivs_a(i);
	
	for(int j=0;j<=i;j++){
	  ValueType change = 2.0 * inv_cov(a,a) * yderivs_a(i) * yderivs_a(j);
	  second_derivs(i,j) += change;
	  if(i!=j) second_derivs(j,i) += change;
	}
      }
      //b<a
      for(int b=0;b<a;b++){	
	if(inv_cov(a,b)==0.0) continue;

	const ValueDerivativeType & yderivs_b = *yderivs[b];
	
	for(int i=0;i<nparams;i++){
	  derivs(i) +=  -2.0*inv_cov(a,b)*( yderivs_a(i)*dfw[b] + dfw[a]* yderivs_b(i) );
	
	  for(int j=0;j<=i;j++){
	    ValueType change = 2.0*inv_cov(a,b)*( yderivs_a(i)* yderivs_b(j) + yderivs_a(j)*yderivs_b(i) );
	    second_derivs(i,j) += change;
	    if(i!=j) second_derivs(j,i) += change;
	  }
	  
	}
      }
    }
  }

  inline int Ndof() const{
    return data.size() - fitfunc.Nparams();
  }
  inline int Nparams() const{ return fitfunc.Nparams(); }    

};

SARLAC_END_NAMESPACE
#endif

#ifndef _COST_FUNCTION_H
#define _COST_FUNCTION_H

#include<cassert>
#include<vector>
#include<memory>
#include<numeric_tensors.h>
#include<utils.h>
/*
  The Marquardt-Levenberg algorithm is a combination of the Gauss-Newton method and gradient descent.
  Gradient descent updates along the direction - grad(f) where f is the fitness function, and thus doesnt need the second derivatives.
  The Gauss-Newton method taylor expands the gradient vector
  grad f(x) = grad f(x0) + (x-x0)^T grad^2 f(x0) + higher order
  and solve for minimum grad f(x) = 0 at some x, then  x_i+1 = x_i - [grad^2 f(x_i)]^-1 grad f(x_i)
  
  as it cuts the grad off at first order, this is equivalent to approximating the function as a quadratic. 

  As  chisq = sum_i ( Y_i - fit( {X}_i ) )^2/sigma_i^2
  then approximating this as quadratic in the parameters {X} means that the fit function is approximated as linear.
  
  As a result we can treat the second derivatives of the fit function as zero without losing anything.
  
  We thus do not bother going to the effort of determining the second derivatives.
  cf. http://www.scribd.com/doc/10093320/Levenberg-Marquardt-Algorithm  (page 2)
*/

template<typename CostType>
struct CostFunctionSVDinvert{
  inline static NumericSquareMatrix<CostType> invert(const NumericSquareMatrix<CostType> &in){
    NumericSquareMatrix<CostType> out(in.size()); svd_inverse(out,in); return out;
  }
};

template<typename FitFunction, typename DataContainer, typename _CostType = double,
	 typename _CostDerivativeType = typename FitFunction::ValueDerivativeType,
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








template<typename FitFunction, typename DataContainer,
	 typename InvCorrMatrixType = NumericSquareMatrix<double>,
	 typename _CostType = double,
	 typename _CostDerivativeType = typename FitFunction::ValueDerivativeType,
	 typename _CostSecondDerivativeMatrixType = NumericSquareMatrix<_CostType>,
	 typename _CostSecondDerivativeInverseMatrixType = NumericSquareMatrix<_CostType>,
	 typename MatrixInvertPolicy = CostFunctionSVDinvert<_CostType>,
	 typename std::enable_if<std::is_floating_point<typename FitFunction::ValueType>::value, int>::type = 0>
class CorrelatedChisqCostFunction: public MatrixInvertPolicy{
  static_assert(std::is_same<typename FitFunction::ValueType, typename DataContainer::DataType>::value, "DataContainer and FitFunction must have same value type");
  static_assert(std::is_same<typename FitFunction::GeneralizedCoordinate, typename DataContainer::GeneralizedCoordinate>::value, "DataContainer and FitFunction must have same coordinate");
  const FitFunction &fitfunc;
  const DataContainer &data;
  const std::vector<double> &sigma; //diagonal elements of covariance matrix
  const InvCorrMatrixType &inv_corr; //inverse correlation matrix
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

  CorrelatedChisqCostFunction(const FitFunction &ff, const DataContainer &dd, const std::vector<double> &_sigma, const InvCorrMatrixType &_inv_corr):
    fitfunc(ff), data(dd), sigma(_sigma), inv_corr(_inv_corr){
    assert(inv_corr.size() == data.size());
    assert(sigma.size() == data.size());
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
  }

  inline int Ndof() const{
    return data.size() - fitfunc.Nparams();
  }
    

};




#endif

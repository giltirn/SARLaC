#include<minimizer/minimizer.h>
#include<containers/single_value_container.h>
#include<tensors/numeric_square_matrix.h>

SARLAC_START_NAMESPACE

//Function should be a lambda or function object with operator(const NumericVector<double> &param) const
//Should return a NumericVector<double> for each of the equations
//Solution is for all f_i(\vec p) = 0
//The number of equations should be equal to the number of parameters
template<typename Function>
class solveMultiDcost{
public:
  typedef double CostType;
  typedef NumericVector<double> ValueType;
  typedef NumericVector<double> ParameterType;
  typedef NumericVector<double> ValueDerivativeType; //derivative wrt parameters

  typedef NumericVector<double> CostDerivativeType;
  typedef NumericSquareMatrix<double> CostSecondDerivativeMatrixType;
  typedef NumericSquareMatrix<double> CostSecondDerivativeInverseMatrixType;

  const Function &func;
  const double deriv_delta;

  solveMultiDcost(const Function &func, const double deriv_delta = 1e-5): func(func), deriv_delta(deriv_delta){
  }

  inline NumericVector<double> value(const ParameterType &p) const{
    NumericVector<double> t = func(p);
    assert(t.size() == p.size());
    return t;
  }
  inline NumericSquareMatrix<double> deriv(const ParameterType &p) const{  
    NumericSquareMatrix<double> dd(p.size());
    for(int i=0;i<p.size();i++){
      ParameterType pplus(p), pminus(p);
      pplus(i) = p(i) + deriv_delta;
      pminus(i) = p(i) - deriv_delta;
      NumericVector<double> vplus = func(pplus);
      NumericVector<double> vminus = func(pminus);
      for(int j=0;j<vplus.size();j++)
	dd(j,i) = (vplus(j) - vminus(j))/ (2*deriv_delta); //[term][param idx]
    }
    return dd;
  }

  inline CostType cost(const ParameterType &params) const{
    NumericVector<double> v = value(params);
    assert(v.size() == params.size());
    double out = 0.;
    for(int t=0;t<v.size();t++) out += v[t]*v[t];
    return out;
  }
  inline void derivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{
    NumericVector<double> v = value(params);
    assert(v.size() == params.size());
    NumericSquareMatrix<double> dterms = deriv(params); //[term][param idx]
    
    derivs.resize(params.size());
    second_derivs.resize(params.size());

    derivs.zero();
    second_derivs.zero();
    double cost = 0.;

    for(int t=0;t<v.size();t++){ //terms
      cost += v[t]*v[t];
      for(int i=0;i<v.size();i++){
	derivs(i) += 2*v[t]*dterms(t,i);
	for(int j=0;j<v.size();j++)
	  second_derivs(i,j) += 2*dterms(t,i)*dterms(t,j);
      }
    }
  }
  inline CostSecondDerivativeMatrixType invert(const CostSecondDerivativeMatrixType &M) const{
    CostSecondDerivativeMatrixType out(M);
    svd_inverse(out, M);
    return out;
  }
  
};

struct solveMultiDparams{
  bool verbose;
  double deriv_delta;
  double tolerance;
  
  solveMultiDparams(): verbose(false), deriv_delta(1e-4), tolerance(1e-4){}
};


template<typename Function>
NumericVector<double> solveMultiD(const Function &func, const NumericVector<double> &guess, const solveMultiDparams &params, double *cost_out = NULL){
  solveMultiDcost<Function> costfunc(func, params.deriv_delta);
  MarquardtLevenbergParameters<double> mlparams;
  mlparams.verbose = params.verbose;
  mlparams.delta_cost_min = params.tolerance * params.tolerance;
  MarquardtLevenbergMinimizer<solveMultiDcost<Function> > min(costfunc, mlparams);
  NumericVector<double> p(guess);
  double cost = min.fit(p);
  if(cost_out) *cost_out = cost;
  return p;
}
template<typename Function>
inline NumericVector<double> solveMultiD(const Function &func, const NumericVector<double> &guess, double *cost_out = NULL){
  solveMultiDparams params;
  return solveMultiD(func, guess, params, cost_out);
}

SARLAC_END_NAMESPACE

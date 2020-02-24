#include<minimizer/minimizer.h>
#include<containers/single_value_container.h>
#include<tensors/numeric_square_matrix.h>

CPSFIT_START_NAMESPACE

//A simple one-dimensional solver for the value of x that satisfies  f(x) = 0

//Function should be a lambda or function object with operator(const double param) const
template<typename Function>
class solveOneDcost{
public:
  typedef double CostType;
  typedef double ValueType;
  typedef singleValueContainer<double> ParameterType;
  typedef singleValueContainer<double> ValueDerivativeType; //derivative wrt parameters

  typedef singleValueContainer<double> CostDerivativeType;
  typedef NumericSquareMatrix<double> CostSecondDerivativeMatrixType;
  typedef NumericSquareMatrix<double> CostSecondDerivativeInverseMatrixType;

  const Function &func;
  const double deriv_delta;

  solveOneDcost(const Function &func, const double deriv_delta = 1e-5): func(func), deriv_delta(deriv_delta){
  }

  inline double value(const ParameterType &p) const{
    return func(p(0));      
  }
  inline double deriv(const ParameterType &p) const{
    double vplus = func(p(0)+deriv_delta);
    double vminus = func(p(0)-deriv_delta);
    return (vplus - vminus)/ (2*deriv_delta);
  }


  inline CostType cost(const ParameterType &params) const{
    double v = value(params);
    return v * v;
  }
  inline void derivatives(CostDerivativeType &derivs, CostSecondDerivativeMatrixType &second_derivs, const ParameterType &params) const{
    double v = value(params);
    double vderiv = deriv(params);
    derivs = singleValueContainer<double>(2.*v*vderiv);
    second_derivs.resize(1);
    second_derivs(0,0) = 2. * vderiv * vderiv;
  }
  inline CostSecondDerivativeMatrixType invert(const CostSecondDerivativeMatrixType &M) const{
    CostSecondDerivativeMatrixType out(1);
    out(0,0) = 1./M(0,0);
    return out;
  }
  
};

struct solveOneDparams{
  bool verbose;
  double deriv_delta;
  double tolerance;
  
  solveOneDparams(): verbose(false), deriv_delta(1e-4), tolerance(1e-4){}
};


template<typename Function>
double solveOneD(const Function &func, const double guess, const solveOneDparams &params, double *cost_out = NULL){
  solveOneDcost<Function> costfunc(func, params.deriv_delta);
  MarquardtLevenbergParameters<double> mlparams;
  mlparams.verbose = params.verbose;
  mlparams.delta_cost_min = params.tolerance * params.tolerance;
  MarquardtLevenbergMinimizer<solveOneDcost<Function> > min(costfunc, mlparams);
  singleValueContainer<double> p(guess);
  double cost = min.fit(p);
  if(cost_out) *cost_out = cost;

  assert(fabs(func(p(0))) <= params.tolerance);
  return p(0);
}
template<typename Function>
inline double solveOneD(const Function &func, const double guess, double *cost_out = NULL){
  solveOneDparams params;
  return solveOneD(func, guess, params, cost_out);
}

CPSFIT_END_NAMESPACE

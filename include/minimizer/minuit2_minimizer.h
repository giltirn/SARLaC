#ifndef _MINUIT2_MINIMIZER_H_
#define _MINUIT2_MINIMIZER_H_

#include<config.h>

#ifdef HAVE_MINUIT2

#include<omp.h>
#include<sstream>

#include<utils/macros.h>
#include<utils/utils.h>
#include<parser/parser.h>
#include<tensors/numeric_square_matrix.h>

#include<Minuit2/FCNGradientBase.h>
#include<Minuit2/FunctionMinimum.h>
#include<Minuit2/MnMinimize.h>
#include<Minuit2/MnSimplex.h>
#include<Minuit2/MnPrint.h>

SARLAC_START_NAMESPACE

//A wrapper around the Minuit2 minimizer
//The cost functions are the usual ones that work with the native minimizer

template<typename CostFunction>
class Minuit2costFunctionWrapper: public ROOT::Minuit2::FCNGradientBase{
private:
  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;
  typedef typename CostFunction::CostDerivativeType CostDerivativeType;
  typedef typename CostFunction::CostSecondDerivativeMatrixType CostSecondDerivativeMatrixType;

  const CostFunction &func;
  ParameterType pbase; //a properly setup instance of the parameter type. Should be resized appropriately
  int nparams;
public:
  Minuit2costFunctionWrapper(const CostFunction &func): func(func), nparams(func.Nparams()){}

  const CostFunction &getFunc() const{ return func; }  

  void setParameterBase(const ParameterType &pb){ assert(pb.size() == nparams); pbase = pb; }

  inline ParameterType getFuncParams(const std::vector<double> &x) const{
    assert(x.size() == nparams);
    ParameterType fparams(pbase);
    for(size_t i=0;i<fparams.size();i++) fparams(i) = x[i];
    return fparams;
  }

  double operator()(const std::vector<double>& x) const{
    return func.cost(getFuncParams(x));
  }

  bool CheckGradient() const{ return false; }

  std::vector<double> Gradient(const std::vector<double>&pin) const{
    ParameterType fparams = getFuncParams(pin);
    CostDerivativeType d;
    CostSecondDerivativeMatrixType dd;
    func.derivatives(d,dd,fparams);    
    std::vector<double> out(nparams);
    for(int j=0;j<nparams;j++)
      out[j] = d(j);
    return out;
  }

   /**
      Error definition of the function. MINUIT defines Parameter errors as the 
      change in Parameter Value required to change the function Value by up. Normally, 
      for chisquared fits it is 1, and for negative log likelihood, its Value is 0.5.
      If the user wants instead the 2-sigma errors for chisquared fits, it becomes 4, 
      as Chi2(x+n*sigma) = Chi2(x) + n*n.
   */

  double Up() const{ return 1.0; }  

  //The error matrix is twice the inverse of the Hessian
  NumericSquareMatrix<double> ErrorMatrix(const std::vector<double>&pin) const{
    ParameterType fparams = getFuncParams(pin);
    CostDerivativeType d;
    CostSecondDerivativeMatrixType dd;
    func.derivatives(d,dd,fparams);    
    NumericSquareMatrix<double> H(nparams);
    for(int i=0;i<nparams;i++)
      for(int j=0;j<nparams;j++)
	H(i,j) = dd(i,j);
    NumericSquareMatrix<double> Hinv(nparams);
    svd_inverse(Hinv,H);
    for(int i=0;i<nparams;i++)
      for(int j=0;j<nparams;j++)
	Hinv(i,j) = Hinv(i,j) * 2.;
    return Hinv;
  }

};


struct Minuit2minimizerParams{  
  int max_iter;
  bool exit_on_convergence_fail;

  //Stopping condition is when estimated distance to minimum is equal to tolerance
  double tolerance;
 
  //Sets the strategy to be used in calculating first and second derivatives and in certain minimization methods. In general, low values of level mean fewer function calls and high values mean more reliable minimization. Currently allowed values are 0 (low), 1 (default), and 2 (high).
  int strategy; 

  bool run_initial_simplex; //perform an initial simplex solve to refine the guess
  double initial_simplex_tolerance;
  int initial_simplex_max_iter;

  bool verbose; 

  Minuit2minimizerParams(): max_iter(10000), exit_on_convergence_fail(true), tolerance(1e-5), strategy(1), verbose(false), run_initial_simplex(false), initial_simplex_tolerance(1e-1), initial_simplex_max_iter(200){}
};

#define MIN2MIN_PARAMS			\
  (double, tolerance)(int, strategy)(bool, run_initial_simplex)(double, initial_simplex_tolerance)(int, initial_simplex_max_iter)(int, max_iter)(bool, verbose)(bool, exit_on_convergence_fail)

GENERATE_PARSER( Minuit2minimizerParams, MIN2MIN_PARAMS );
		    
#undef MIN2MIN_PARAMS


template<typename CostFunction> 
class Minuit2minimizer{  
public:
  typedef Minuit2minimizerParams AlgorithmParameterType;
private:

  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;

  Minuit2minimizerParams min_params;

  typedef Minuit2costFunctionWrapper<CostFunction> WrapperType;
  
  WrapperType func_wrapper;
  int nparams;

  int me;
  std::string prefix;

  bool converged;
  int iter;
  
public:

  Minuit2minimizer(const CostFunction &func, const AlgorithmParameterType &min_params): func_wrapper(func), nparams(func.Nparams()), min_params(min_params){
    me = 0;
    std::ostringstream pf;
    if(omp_in_parallel()){
      me = omp_get_thread_num();
      pf << "Thread " << me << " ";
    }
    pf << "Minuit2minimizer: ";
    prefix = pf.str();
  }
  ~Minuit2minimizer(){
  }

  
  CostType fit(ParameterType &params){
    using namespace ROOT::Minuit2;

    if(min_params.verbose)
      MnPrint::SetLevel(10);

    assert(params.size() == nparams);
    func_wrapper.setParameterBase(params);

    std::vector<double> init_par(nparams);
    for(int i=0;i<nparams;i++) init_par[i] = params(i);

    std::vector<double> init_err(nparams);
    NumericSquareMatrix<double> init_errm = func_wrapper.ErrorMatrix(init_par);

    for(int i=0;i<nparams;i++) init_err[i] = sqrt(init_errm(i,i));

    if(min_params.run_initial_simplex){//perform an initial simplex solve to refine the guess
      double f_init = func_wrapper(init_par);

      MnSimplex min(func_wrapper, init_par, init_err);

      FunctionMinimum sol = min(min_params.initial_simplex_max_iter, min_params.initial_simplex_tolerance);
    
      int simplex_iter = min.NumOfCalls();
      bool simplex_converged = sol.IsValid();

      std::vector<double> sol_params(nparams);
      for(int i=0;i<nparams;i++) sol_params[i] = sol.Parameters().Vec()(i);
      
      double f_final = func_wrapper(sol_params);

      if(simplex_converged || f_final < f_init){
	if(!me)
	  if(simplex_converged) std::cout << prefix << " Initial simplex converged\n";
	  else std::cout << prefix << " Initial simplex did not converge but did improve the guess\n";
	
	init_par = sol_params;
	init_errm = func_wrapper.ErrorMatrix(init_par);
	for(int i=0;i<nparams;i++) init_err[i] = sqrt(init_errm(i,i));
      }else if(!me) std::cout << "Initial simplex did not converge and did not improve guess: reverting to initial guess\n";
    }

    //MnMigrad min(func_wrapper, init_par, init_err, min_params.strategy);
    MnMinimize min(func_wrapper, init_par, init_err, min_params.strategy);
  
    iter = 0;
    converged = false;

    double tolerance = min_params.tolerance/1e-3; //this is multiplied by 1e-3 internally
    
    FunctionMinimum sol = min(min_params.max_iter, tolerance);
    
    iter = min.NumOfCalls();

    converged = sol.IsValid();

    if(!converged && min_params.exit_on_convergence_fail) error_exit(std::cerr << prefix << "Failed to converge\n");    

    for(int i=0;i<nparams;i++)
      params(i) = sol.Parameters().Vec()(i);

    return sol.Fval();
  }
  inline bool hasConverged() const{ return converged; }
  inline int iterations() const{ return iter; }
};


SARLAC_END_NAMESPACE

#endif

#endif

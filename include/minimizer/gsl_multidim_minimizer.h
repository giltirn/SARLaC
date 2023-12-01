#ifndef _GSL_MULTIDIM_MINIMIZER_H_
#define _GSL_MULTIDIM_MINIMIZER_H_

//A wrapper around the GSL multidimensional minimizers (gsl_multimin)
//The cost functions are the usual ones that work with the native minimizer

#include<gsl/gsl_multimin.h>

#include<minimizer/minimizer.h>

SARLAC_START_NAMESPACE

//First line of algorithms use the derivative. The FR and PR conjugate gradient algorithms are Fletcher-Reeves and Polak-Ribiere, respectively
//Second line do not use the derivative. These are variants of the Nedler-Mead simplex algorithm
//For more explanation see 
//https://www.gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-with-Derivatives.html#Multimin-Algorithms-with-Derivatives           [1]
//https://www.gnu.org/software/gsl/manual/html_node/Multimin-Algorithms-without-Derivatives.html#Multimin-Algorithms-without-Derivatives        [2]
GENERATE_ENUM_AND_PARSER(GSLmultiminAlgorithm, 
			 (ConjugateGradientFR)(ConjugateGradientPR)(VectorBFGS)(VectorBFGS2)(SteepestDescent)
			 (NMsimplex)(NMsimplex2)(NMsimplex2rand) );
GENERATE_ENUM_AND_PARSER(GSLmultiminStoppingCondition, (StopGradient)(StopCostDelta)(StopRelativeGradients)(StopSimplexSize) );

typedef std::pair<GSLmultiminStoppingCondition, double> GSLmultiminStoppingConditionAndValue;

//Stopping conditions:
//StopGradient:  Stop when  |g| < value
//StopCostDelta: Stop when absolute change in cost between two successive iterations is < value
//StopRelativeGradients:   Stop when   df/dp_i  * p_i/f  < value  for all parameter indices i
//StopSimplexSize: Only applies to simplex algorithms, related to separation between simplex vertices https://www.gnu.org/software/gsl/doc/html/multimin.html

struct GSLmultidimMinimizerParams{
  GSLmultiminAlgorithm algorithm;
  double line_search_tol; //use only by derivative solvers cf https://www.gnu.org/software/gsl/manual/html_node/Initializing-the-Multidimensional-Minimizer.html#Initializing-the-Multidimensional-Minimizer
  std::vector<double> step_size; //For derivative solvers this should have size 1, otherwise it should have size = nparams
                                 //The exact meaning depends on the algorith, cf [1] and [2] above

  std::vector<GSLmultiminStoppingConditionAndValue> stopping_conditions; //choose one or more stopping conditions and their values

  int max_iter;
  bool verbose;
  std::ostream *output;

  bool exit_on_convergence_fail;

GSLmultidimMinimizerParams(): algorithm(GSLmultiminAlgorithm::ConjugateGradientFR), line_search_tol(0.1), step_size(1, 0.1), stopping_conditions({ {GSLmultiminStoppingCondition::StopCostDelta, 1e-5} }), max_iter(10000), verbose(false), output(&std::cout), exit_on_convergence_fail(true){}
};

#define GSLMDIM_PARAMS			\
  (GSLmultiminAlgorithm, algorithm)(double, line_search_tol)(std::vector<double>, step_size)(std::vector<GSLmultiminStoppingConditionAndValue>, stopping_conditions) \
  (int, max_iter)(bool, verbose)(bool, exit_on_convergence_fail)

GENERATE_PARSER( GSLmultidimMinimizerParams, GSLMDIM_PARAMS );
		    
#undef GSLMDIM_PARAMS


//A function wrapper for interfacing my cost functions with the GSL function definition for use in the multidimensional minimization (gsl_multimin)
template<typename CostFunction>
class GSLmultiminCostFunctionWrapper{ 
private:
  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;
  typedef typename CostFunction::CostDerivativeType CostDerivativeType;
  typedef typename CostFunction::CostSecondDerivativeMatrixType CostSecondDerivativeMatrixType;

  const CostFunction &func;
  ParameterType pbase; //a properly setup instance of the parameter type. Should be resized appropriately
  int nparams;
public:
  GSLmultiminCostFunctionWrapper(const CostFunction &func): func(func), nparams(func.Nparams()){}

  const CostFunction &getFunc() const{ return func; }  

  void setParameterBase(const ParameterType &pb){ assert(pb.size() == nparams); pbase = pb; }

  inline ParameterType getFuncParams(const gsl_vector * x) const{
    assert(x->size == nparams);
    ParameterType fparams(pbase);
    for(size_t i=0;i<fparams.size();i++) fparams(i) = gsl_vector_get(x,i);
    return fparams;
  }

  //For explanations of these functions look at
  //https://www.gnu.org/software/gsl/manual/html_node/Providing-a-function-to-minimize.html#Providing-a-function-to-minimize


  //x - the trial fit parameters
  //params - an arbitrary container for data that the function can use (size p)

  inline double f(const gsl_vector * x, void * params) const{ 
    ParameterType fparams = getFuncParams(x);
    return func.cost(fparams);
  }

  //g_i  = df/dp_i
  inline void df(const gsl_vector * x, void * params,  gsl_vector * g) const{
    assert(g->size == nparams);
    ParameterType fparams = getFuncParams(x);
    CostDerivativeType d;
    CostSecondDerivativeMatrixType dd;
    func.derivatives(d,dd,fparams);    
    for(int j=0;j<nparams;j++)
      gsl_vector_set(g, j, d(j)); 
  }

  inline void fdf(const gsl_vector * x, void * params, double* f, gsl_vector * g) const{
    assert(g->size == nparams);
    ParameterType fparams = getFuncParams(x);
    CostDerivativeType d;
    CostSecondDerivativeMatrixType dd;
    func.derivatives(d,dd,fparams);    
    for(int j=0;j<nparams;j++)
      gsl_vector_set(g, j, d(j)); 
    *f = func.cost(fparams);
  }
  
};

template<typename Wrapper>
inline double evalGSLmultidimCostFunctionWrapper_f(const gsl_vector * x, void * params){ 
  Wrapper *w = reinterpret_cast<Wrapper *>(params);
  return w->f(x,NULL);
}
template<typename Wrapper>
inline void evalGSLmultidimCostFunctionWrapper_df(const gsl_vector * x, void * params, gsl_vector * g){
  Wrapper *w = reinterpret_cast<Wrapper *>(params);
  w->df(x,NULL,g);
}
template<typename Wrapper>
inline void evalGSLmultidimCostFunctionWrapper_fdf(const gsl_vector * x, void * params, double* f, gsl_vector * g){
  Wrapper *w = reinterpret_cast<Wrapper *>(params);
  w->fdf(x,NULL,f,g);
}



template<typename CostFunction> 
class GSLmultidimMinimizer{  
public:
  typedef GSLmultidimMinimizerParams AlgorithmParameterType;
private:

  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;

  GSLmultidimMinimizerParams min_params;

  gsl_multimin_fdfminimizer * min_fdf;
  gsl_multimin_fminimizer * min_f;

  typedef GSLmultiminCostFunctionWrapper<CostFunction> WrapperType;
  
  WrapperType func_wrapper;
  int nparams;
  gsl_multimin_function_fdf func_fdf;
  gsl_multimin_function func_f;

  int iter;
  bool converged;

  int me;
  std::string prefix;

  void setup_functions(){
    func_fdf.f = evalGSLmultidimCostFunctionWrapper_f<WrapperType>;
    func_fdf.df = evalGSLmultidimCostFunctionWrapper_df<WrapperType>;
    func_fdf.fdf = evalGSLmultidimCostFunctionWrapper_fdf<WrapperType>;
    func_fdf.n = nparams;
    func_fdf.params = (void*)&func_wrapper;

    func_f.f = evalGSLmultidimCostFunctionWrapper_f<WrapperType>;
    func_f.n = nparams;
    func_f.params = (void*)&func_wrapper;
  }

  //is an algorithm requiring only the function
  static bool isFalg(GSLmultiminAlgorithm alg){
    switch(alg){
    case GSLmultiminAlgorithm::ConjugateGradientFR:
    case GSLmultiminAlgorithm::ConjugateGradientPR:
    case GSLmultiminAlgorithm::VectorBFGS:
    case GSLmultiminAlgorithm::VectorBFGS2:
    case GSLmultiminAlgorithm::SteepestDescent:
      return false;
    case GSLmultiminAlgorithm::NMsimplex:
    case GSLmultiminAlgorithm::NMsimplex2:
    case GSLmultiminAlgorithm::NMsimplex2rand:
      return true;
    default:
      assert(0);
    }
  }

  void initializeAlgorithm(){
    switch(min_params.algorithm){
    case GSLmultiminAlgorithm::ConjugateGradientFR:
      min_fdf = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_fr, nparams); break;
    case GSLmultiminAlgorithm::ConjugateGradientPR:
      min_fdf = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_conjugate_pr, nparams); break;
    case GSLmultiminAlgorithm::VectorBFGS:
      min_fdf = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs, nparams); break;
    case GSLmultiminAlgorithm::VectorBFGS2:
      min_fdf = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_vector_bfgs2, nparams); break;
    case GSLmultiminAlgorithm::SteepestDescent:
      min_fdf = gsl_multimin_fdfminimizer_alloc(gsl_multimin_fdfminimizer_steepest_descent, nparams); break;
    case GSLmultiminAlgorithm::NMsimplex:
      min_f = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex, nparams); break;
    case GSLmultiminAlgorithm::NMsimplex2:
      min_f = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2, nparams); break;
    case GSLmultiminAlgorithm::NMsimplex2rand:
      min_f = gsl_multimin_fminimizer_alloc(gsl_multimin_fminimizer_nmsimplex2rand, nparams); break;
    default:
      assert(0);
    }
  }

public:

 GSLmultidimMinimizer(const CostFunction &func, const AlgorithmParameterType &min_params): func_wrapper(func), nparams(func.Nparams()), min_params(min_params), min_fdf(NULL), min_f(NULL){
    setup_functions();
    initializeAlgorithm();

    me = 0;
    std::ostringstream pf;
    if(omp_in_parallel()){
      me = omp_get_thread_num();
      pf << "Thread " << me << " ";
    }
    pf << "GSLmultidimMinimizer(" << min_params.algorithm << "): ";
    prefix = pf.str();
  }
  ~GSLmultidimMinimizer(){
    if(min_f != NULL) gsl_multimin_fminimizer_free(min_f);
    if(min_fdf != NULL) gsl_multimin_fdfminimizer_free(min_fdf);
  }


  CostType fit(ParameterType &params){
    assert(params.size() == nparams);
    func_wrapper.setParameterBase(params);

    gsl_vector *guess = gsl_vector_alloc(nparams);
    for(int p=0;p<nparams;p++) gsl_vector_set(guess, p, params(p));

    bool is_falg = isFalg(min_params.algorithm);
    gsl_vector *step_size;
    int err;
    if(is_falg){
      assert(min_params.step_size.size() == nparams);
      step_size = gsl_vector_alloc(nparams);
      for(int i=0;i<nparams;i++) gsl_vector_set(step_size, i, min_params.step_size[i]);

      err = gsl_multimin_fminimizer_set(min_f, &func_f, guess, step_size);
    }else{
      assert(min_params.step_size.size() == 1);
      err = gsl_multimin_fdfminimizer_set(min_fdf, &func_fdf, guess, min_params.step_size[0], min_params.line_search_tol);
    }
    if(err != GSL_SUCCESS) error_exit(std::cout << prefix << "Unable to initialize solver due to error: " << gsl_strerror(err) << std::endl);

    iter = 0;
    converged = false;
    double prev_cost;
     
    //Setup stopping conditions
    assert(min_params.stopping_conditions.size() > 0);
    
    bool stop_dcost = false;
    double stop_cond_dcost;
    bool stop_grad = false;
    double stop_cond_grad;
    bool stop_rel_grad = false;
    double stop_cond_rel_grad;
    bool stop_simplex_size = false;
    double stop_cond_simplex_size;
    for(int i=0;i<min_params.stopping_conditions.size();i++){
      switch(min_params.stopping_conditions[i].first){
      case GSLmultiminStoppingCondition::StopGradient:
	stop_grad = true; stop_cond_grad = min_params.stopping_conditions[i].second; break;
      case GSLmultiminStoppingCondition::StopCostDelta:
	stop_dcost = true; stop_cond_dcost = min_params.stopping_conditions[i].second; break;
      case GSLmultiminStoppingCondition::StopRelativeGradients:
	stop_rel_grad = true; stop_cond_rel_grad = min_params.stopping_conditions[i].second; break;
      case GSLmultiminStoppingCondition::StopSimplexSize:
	stop_simplex_size = true; stop_cond_simplex_size = min_params.stopping_conditions[i].second; break;
      default:
	assert(0);
      }
    }

    if(stop_simplex_size && !is_falg) error_exit(std::cout << "Simplex size stopping condition only applies to simplex algorithms!\n");

    gsl_vector* g_temp = gsl_vector_alloc(nparams);
    double relg[nparams];

    double const* cost_ptr = is_falg ? &min_f->fval : &min_fdf->f;
    gsl_vector const* x =  is_falg ? min_f->x : min_fdf->x;
    gsl_vector const* g = is_falg ? g_temp : min_fdf->gradient;

    //Iterate
    for(iter = 0; iter <= min_params.max_iter; iter++){
      double cost = *cost_ptr;
      double dcost = cost - prev_cost;
      
      if(is_falg) func_wrapper.df(min_f->x, NULL, g_temp); //compute the derivative explicitly as we know it!

      for(int p=0;p<nparams;p++) relg[p] = gsl_vector_get(g,p) * gsl_vector_get(x,p)/ cost;

      double modg = gsl_blas_dnrm2(g);

      if(!me && min_params.verbose){
	auto &op = *min_params.output;
	op << prefix << iter << " Cost=" << cost << " dCost=" << dcost << " |g|=" << modg;
	if(is_falg) op << " simplex size=" << gsl_multimin_fminimizer_size(min_f);
	if(stop_dcost) op << " dCost/tol=" << dcost/stop_cond_dcost;
	if(stop_grad) op << " |g|/tol=" << modg/stop_cond_grad;
	op << " (p_i/f * df/dp_i)=(";
	for(int p=0;p<nparams;p++) op << relg[p] << (p == nparams - 1 ? "" : ", ");
	if(stop_rel_grad){
	  op << ") |p_i/f * df/dp_i|/tol=(";
	  for(int p=0;p<nparams;p++) op << fabs(relg[p])/stop_cond_rel_grad << (p == nparams - 1 ? "" : ", ");
	}
	op << ") Parameters=(";
	for(int p=0;p<nparams;p++) op << gsl_vector_get(x,p) << (p == nparams - 1 ? "" : ", ");
	op << ")" << std::endl;
      }
      
      if(stop_dcost && iter > 0 && dcost < 0 && -dcost < stop_cond_dcost){
	converged=true; break;
      }      
      if(stop_grad && modg < stop_cond_grad){
	converged=true; break;
      }
      if(stop_rel_grad){ // Stop when   df/dp_i  * p_i/f  < value  for all parameter indices i
	converged=true;
	for(int p=0;p<nparams;p++) if( fabs(relg[p]) > stop_cond_rel_grad ){ converged = false; break; }
	if(converged) break;
      }            
      if(stop_simplex_size && gsl_multimin_fminimizer_size(min_f) < stop_cond_simplex_size){
	converged = true; break;
      }
      
      if(iter == min_params.max_iter){
	if(min_params.exit_on_convergence_fail) error_exit(std::cout << prefix << "Reached max iterations\n");
	else break;
      }

      prev_cost = cost;

      //Perform the iteration
      int err = is_falg ? gsl_multimin_fminimizer_iterate(min_f) : gsl_multimin_fdfminimizer_iterate(min_fdf);
      if(err != GSL_SUCCESS) error_exit(std::cout << prefix << "Unable to iterate due to error: " << gsl_strerror(err) << std::endl);

    }
    if(is_falg) gsl_vector_free(step_size);
    gsl_vector_free(guess);
    gsl_vector_free(g_temp);

    for(int i=0;i<nparams;i++) params(i) = gsl_vector_get(x, i);
    return *cost_ptr;
  }
  inline bool hasConverged() const{ return converged; }
  inline int iterations() const{ return iter; }
};

SARLAC_END_NAMESPACE

#endif

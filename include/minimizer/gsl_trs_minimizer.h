#ifndef _GSL_TRS_MINIMIZER_H_
#define _GSL_TRS_MINIMIZER_H_

//A wrapper around the GSL trust region solvers
//Note the cost function is defined somewhat differently here as the GSL minimizer takes the individual terms in the chi^2 and the associated Jacobian

#include <gsl/gsl_multifit_nlinear.h>

#include<minimizer/minimizer.h>

CPSFIT_START_NAMESPACE


  //Stopping conditions
  /*
    gsl_multilarge_nlinear_test()
    Convergence tests for nonlinear least squares minimization
    (1) |dx_i| <= xtol * (1 + |x_i|) for all i   <-StopParamDelta
    (2) || g .* x ||_inf <= gtol ||f||^2   <-StopGradient
    (3) ||f(x+dx) - f(x)|| <= ftol * max(||f(x)||, 1)  //THIS ONE APPEARS TO BE DEPRECATED
    Inputs: xtol - tolerance for step size
    gtol - tolerance for gradient vector
    ftol - tolerance for residual vector
    info - (output)
               1 - stopped by small x step
               2 - stopped by small gradient
               3 - stopped by small residual vector change
	       w    - workspace

    I added StopCostDelta to match the stopping condition used in my ML implementation
  */

GENERATE_ENUM_AND_PARSER(GSLtrsStoppingCondition, (StopParamDelta)(StopGradient)(StopCostDelta) );
GENERATE_ENUM_AND_PARSER(GSLtrsAlgorithm, (MarquardtLevenberg)(MarquardtLevenbergGeodesicAccel)(Dogleg)(DoubleDogleg)(Subspace2D) );

typedef std::pair<GSLtrsStoppingCondition, double> GSLtrsStoppingConditionAndValue;

struct GSLtrsMinimizerParams{
  GSLtrsAlgorithm algorithm;

  std::vector<GSLtrsStoppingConditionAndValue> stopping_conditions; //choose one or more stopping conditions and their values
  
  //Control how the parameter lambda in the Marquardt-Levenberg algorithm is updated
  double factor_up;
  double factor_down; 

  //Choose the Marquardt-Levenberg dampening matrix (cf https://www.gnu.org/software/gsl/doc/html/nls.html)
  MLdampeningMatrix dampening_matrix;

  int max_iter;
  bool verbose;
  std::ostream *output;

  bool exit_on_convergence_fail; //exit with error if failed to converge
  
  GSLtrsMinimizerParams(): algorithm(GSLtrsAlgorithm::MarquardtLevenberg), stopping_conditions({ {GSLtrsStoppingCondition::StopCostDelta, 1e-5} }), 
			   max_iter(10000), exit_on_convergence_fail(true),
			   output(&std::cout), verbose(false), factor_up(10), factor_down(10), dampening_matrix(MLdampeningMatrix::HessianDiag){}
};

#define GSLTRSP_PARAMS			\
  (GSLtrsAlgorithm, algorithm)(std::vector<GSLtrsStoppingConditionAndValue>, stopping_conditions) \
  (double, factor_up)(double, factor_down)(MLdampeningMatrix, dampening_matrix) \
  (int, max_iter)(bool, verbose)(bool, exit_on_convergence_fail)

GENERATE_PARSER( GSLtrsMinimizerParams, GSLTRSP_PARAMS );
		    
#undef GSLTRSP_PARAMS



//A function wrapper for interfacing my cost functions with the GSL function definition for use in the trust region solvers (gsl_multifit_nlinear)
template<typename CostFunction>
class GSLtrsCostFunctionWrapper{ 
private:
  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;

  const CostFunction &func;
  ParameterType pbase; //a properly setup instance of the parameter type. Should be resized appropriately
  int nparams;
  int nterms;
public:
  GSLtrsCostFunctionWrapper(const CostFunction &func): func(func), nparams(func.Nparams()), nterms(func.Nterms()){}

  const CostFunction &getFunc() const{ return func; }  

  void setParameterBase(const ParameterType &pb){ assert(pb.size() == nparams); pbase = pb; }

  inline ParameterType getFuncParams(const gsl_vector * x) const{
    assert(x->size == nparams);
    ParameterType fparams(pbase);
    for(size_t i=0;i<fparams.size();i++) fparams(i) = gsl_vector_get(x,i);
    return fparams;
  }

  //For explanations of these functions look at
  //https://www.gnu.org/software/gsl/manual/html_node/Nonlinear-Least_002dSquares-Function-Definition.html#Nonlinear-Least_002dSquares-Function-Definition

  //x - the trial fit parameters
  //params - an arbitrary container for data that the function can use (size p)
  //f   vector of contributions with  cost = 1/2\sum_i f_i^2

  inline int f(const gsl_vector * x, void * params, gsl_vector * f) const{ 
    ParameterType fparams = getFuncParams(x);
    assert(f->size == nterms);
    auto fv = func.costVector(fparams);
    for(int i=0;i<nterms;i++) gsl_vector_set(f, i, fv(i));
    return GSL_SUCCESS;
  }

  //J is the nterms * nparams jacobian matrix  df_i/dp_j

  int df(const gsl_vector * x, void * params, gsl_matrix * J) const{
    ParameterType fparams = getFuncParams(x);
    assert(J->size1 == nterms);
    assert(J->size2 == nparams);

    auto Jv = func.costJacobianMatrix(fparams);
    for(int i=0;i<nterms;i++)
      for(int j=0;j<nparams;j++)
	gsl_matrix_set(J, i, j, Jv(i,j)); 

    return GSL_SUCCESS; 
  }
   
};

template<typename Wrapper>
inline int evalGSLcostFunctionWrapper_f(const gsl_vector * x, void * params, gsl_vector * f){ 
  Wrapper *w = reinterpret_cast<Wrapper *>(params);
  return w->f(x,NULL,f);
}
template<typename Wrapper>
inline int evalGSLcostFunctionWrapper_df(const gsl_vector * x, void * params, gsl_matrix * J){
  Wrapper *w = reinterpret_cast<Wrapper *>(params);
  return w->df(x,NULL,J);
}




template<typename CostFunction> 
class GSLtrsMinimizer{  
public:
  typedef GSLtrsMinimizerParams AlgorithmParameterType;
private:

  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;

  GSLtrsMinimizerParams min_params;
  gsl_multifit_nlinear_parameters alg_params;
  gsl_multifit_nlinear_workspace *workspace;

  typedef GSLtrsCostFunctionWrapper<CostFunction> WrapperType;
  
  WrapperType func_wrapper;
  int nterms;
  int nparams;
  gsl_multifit_nlinear_fdf fdf;

  int iter;
  bool converged;

  int me;
  std::string prefix;

  void setup_fdf(){
    fdf.f = evalGSLcostFunctionWrapper_f<WrapperType>;
    fdf.df = evalGSLcostFunctionWrapper_df<WrapperType>;
    fdf.fvv = NULL; //used for geodesic acceleration; is computed numerically if set to NULL
    fdf.n = nterms;
    fdf.p = nparams;
    fdf.params = (void*)&func_wrapper;
  }
  void setup_alg_params(){
    alg_params = gsl_multifit_nlinear_default_parameters();
    alg_params.trs = gsl_multifit_nlinear_trs_lm;
    alg_params.factor_up = min_params.factor_up;
    alg_params.factor_down = min_params.factor_down;
    alg_params.solver = gsl_multifit_nlinear_solver_svd; //best quality solver for underlying trust region problem

    switch(min_params.dampening_matrix){
    case(MLdampeningMatrix::HessianDiag):
      alg_params.scale = gsl_multifit_nlinear_scale_marquardt; break;
    case(MLdampeningMatrix::MaxHessianDiag):
      alg_params.scale = gsl_multifit_nlinear_scale_more; break;
    case(MLdampeningMatrix::Unit):
      alg_params.scale = gsl_multifit_nlinear_scale_levenberg; break;
    default:
      assert(0);
    }

    switch(min_params.algorithm){
    case GSLtrsAlgorithm::MarquardtLevenberg:
      alg_params.trs = gsl_multifit_nlinear_trs_lm; break;
    case GSLtrsAlgorithm::MarquardtLevenbergGeodesicAccel:
      alg_params.trs = gsl_multifit_nlinear_trs_lmaccel; break;
    case GSLtrsAlgorithm::Dogleg:
      alg_params.trs = gsl_multifit_nlinear_trs_dogleg; break;
    case GSLtrsAlgorithm::DoubleDogleg:
      alg_params.trs = gsl_multifit_nlinear_trs_ddogleg; break;
    case GSLtrsAlgorithm::Subspace2D:
      alg_params.trs = gsl_multifit_nlinear_trs_subspace2D; break;
    default:
      assert(0);
    }

  }

public:

  GSLtrsMinimizer(const CostFunction &func, const AlgorithmParameterType &min_params): func_wrapper(func), nparams(func.Nparams()), nterms(func.Nterms()), min_params(min_params){
    setup_alg_params();

    workspace = gsl_multifit_nlinear_alloc(gsl_multifit_nlinear_trust, &alg_params, nterms, nparams);
    
    setup_fdf();

    me = 0;
    std::ostringstream pf;
    if(omp_in_parallel()){
      me = omp_get_thread_num();
      pf << "Thread " << me << " ";
    }
    pf << "GSLtrsMinimizer(" << min_params.algorithm << "): ";
    prefix = pf.str();
  }
  ~GSLtrsMinimizer(){
    gsl_multifit_nlinear_free(workspace);
  }


  CostType fit(ParameterType &params){
    assert(params.size() == nparams);
    func_wrapper.setParameterBase(params);

    gsl_vector *guess = gsl_vector_alloc(nparams);
    for(int p=0;p<nparams;p++) gsl_vector_set(guess, p, params(p));

    assert(gsl_multifit_nlinear_init(guess, &fdf, workspace) == GSL_SUCCESS);

    bool stopx=false, stopg=false, stop_dcost=false;
    double xtol=0, gtol=0, dcost_tol=0; //assure it only stops on the conditions we want
    for(auto it=min_params.stopping_conditions.begin(); it != min_params.stopping_conditions.end(); it++){
      switch(it->first){
      case GSLtrsStoppingCondition::StopParamDelta:
	xtol = it->second; stopx=true; break;
      case GSLtrsStoppingCondition::StopGradient:
	gtol = it->second; stopg=true; break;
      case GSLtrsStoppingCondition::StopCostDelta:
	dcost_tol = it->second; stop_dcost=true; break;
      default:
	assert(0);
      }
    }

    iter = 0;
    converged = false;
    double prev_cost;

    for(iter = 0; iter <= min_params.max_iter; iter++){
      CostType cost = 0.5 * pow(gsl_blas_dnrm2(workspace->f),2);

      double dcost = cost - prev_cost;

      //Test for convergence
      int info;
      gsl_multifit_nlinear_test(xtol,gtol,0., &info, workspace);

      if(!me && min_params.verbose){
	auto &op = *min_params.output;
	op << prefix << iter << " Cost=" << cost << " dCost=" << dcost;
	if(stop_dcost) op << " dCost/tol=" << dcost/dcost_tol;
	op << " Parameters=";
	for(int p=0;p<nparams;p++) op << gsl_vector_get(workspace->x,p) << (p == nparams - 1 ? "" : ", ");
	op << ")" << std::endl;
      }

      
      if(info == 1 && stopx){
	converged=true; break;
      }else if(info == 2 && stopg){
	converged=true; break;
      }else if(iter > 0 && dcost < 0 && -dcost < dcost_tol && stop_dcost){ //the algorithm internally will always reduce the cost so no need for fabs
	converged=true; break;
      }else if(iter == min_params.max_iter){
	if(min_params.exit_on_convergence_fail) error_exit(std::cout << prefix << "Reached max iterations\n");
	else break;
      }

      prev_cost = cost;

      //Perform the iteration
      int err = gsl_multifit_nlinear_iterate(workspace);
      if(err != GSL_SUCCESS) error_exit(std::cout << prefix << "Unable to iterate due to error: " << gsl_strerror(err) << std::endl);

    }
    gsl_vector_free(guess);

    for(int i=0;i<nparams;i++) params(i) = gsl_vector_get(workspace->x, i);
    return  0.5 * pow(gsl_blas_dnrm2(workspace->f),2);
  }
  inline bool hasConverged() const{ return converged; }
  inline int iterations() const{ return iter; }
};


CPSFIT_END_NAMESPACE

#endif

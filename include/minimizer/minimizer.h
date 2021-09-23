#ifndef _MINIMIZER_H
#define _MINIMIZER_H

//An implementation of the Maquardt-Levenberg minimizer for a generic cost function

#include<omp.h>
#include<iostream>
#include<sstream>
#include<cstdlib>
#include<type_traits>
#include<complex>
#include<memory>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<parser/parser.h>
#include<tensors/gsl_svdinverse.h>

CPSFIT_START_NAMESPACE

//The choice of the dampening matrix can have an effect on the fit's stability. cf. section 2.2 of arXiv:1201.5885

//Levenberg originally suggested using the unit matrix (Option: Unit). This option however is not invariant under rescaling of the parameters and as a result can converge slowly in a canyon

//Using the diagonal elements of the (approximate) Hessian (Option: HessianDiag) makes the algorithm invariant under this rescaling but at the cost of increased susceptibility to parameter evaporation where the algorithm wanders off to a minimum at infinity

//Minpack uses a compromise method (Option: MaxHessianDiag) where the damping matrix is chosen to be diagonal with entries corresponding to the maximum absolute values observed previously. This improves resistance to parameter evaporation but it still susceptible to poor guesses that do not provide enough damping. (https://hwbdocuments.env.nm.gov/Los%20Alamos%20National%20Labs/General/32166.pdf pg 112)

GENERATE_ENUM_AND_PARSER(MLdampeningMatrix, (Unit)(HessianDiag)(MaxHessianDiag) );

//Parameters of the minimizer
template<typename CostType>
struct MarquardtLevenbergParameters{
  double lambda_factor;
  double lambda_min;
  double lambda_max;
  double lambda_start;
  CostType delta_cost_min;
  
  MLdampeningMatrix dampening_matrix;

  int max_iter;
  bool verbose;
  std::ostream *output;

  bool exit_on_convergence_fail; //exit with error if failed to converge
  
  MarquardtLevenbergParameters(){
    lambda_factor = 10.;
    lambda_min = 1e-05;
    lambda_max = 1e6;
    lambda_start = 0.001;
    delta_cost_min = 1e-05;
    dampening_matrix = MLdampeningMatrix::HessianDiag;
    max_iter = 10000;
    output = &std::cout;
    verbose = false;
    exit_on_convergence_fail = true;
  }
};

#define MLP_DOUBLE_MEM \
  (double, lambda_factor)(double, lambda_min)(double, lambda_max)(double, lambda_start)(double, delta_cost_min) \
  (MLdampeningMatrix, dampening_matrix)(int, max_iter)			\
  (bool, verbose)(bool, exit_on_convergence_fail)

GENERATE_PARSER_GM( MarquardtLevenbergParameters<double>, MarquardtLevenbergParameters_double_grammar, MLP_DOUBLE_MEM)
		    
#undef MLP_DOUBLE_MEM

//Run on a single thread
template<typename CostFunction> 
class MarquardtLevenbergMinimizer{  
  typedef typename CostFunction::CostType CostType;
  typedef typename CostFunction::ParameterType ParameterType;
  typedef typename CostFunction::CostDerivativeType CostDerivativeType;
  typedef typename CostFunction::CostSecondDerivativeMatrixType CostSecondDerivativeMatrixType;
  typedef typename CostFunction::CostSecondDerivativeInverseMatrixType CostSecondDerivativeInverseMatrixType;
public:
  typedef MarquardtLevenbergParameters<CostType> AlgorithmParameterType;
private:
  AlgorithmParameterType mlparams;
  
  const CostFunction &function;

  //Environment
  int me; //thread index
  std::string prefix; //output prefix
  std::streamsize init_precision;
  
  //Algorithm state
  int iter;
  double lambda;
  bool converged;
  bool fail;

  //Cost function previous state
  CostType last_cost;
  std::unique_ptr<ParameterType> last_params;
  CostDerivativeType last_derivs;
  CostSecondDerivativeMatrixType last_second_derivs;
  
  //Cost function current state
  std::unique_ptr<ParameterType> parameters;
  CostType cost;
  bool use_derivs_from_last_step;
  CostDerivativeType derivs;
  CostSecondDerivativeMatrixType second_derivs;

  //Temp data
  CostSecondDerivativeMatrixType M;
  CostSecondDerivativeInverseMatrixType inv_M;
  CostDerivativeType c;
  std::unique_ptr<ParameterType> delta;
  
  CostSecondDerivativeMatrixType D; //Dampening matrix

  void initialize(const ParameterType &init_parameters){
    parameters = std::unique_ptr<ParameterType>(new ParameterType(init_parameters));
    delta = std::unique_ptr<ParameterType>(new ParameterType(init_parameters));
    last_params = std::unique_ptr<ParameterType>(new ParameterType(init_parameters));

    iter = 0;
    converged = false;
    fail = false;
    use_derivs_from_last_step = false;

    lambda = mlparams.lambda_start;
    cost = function.cost(init_parameters);  
    
    prefix = "";
    me = 0;
    if(omp_in_parallel()){
      me = omp_get_thread_num();
      std::ostringstream os; os << "Thread " << me << " ";
      prefix = os.str();
    }
    prefix = prefix + "MarquardtLevenbergMinimizer: ";

#define MINPRINT_NOPREFIX if(mlparams.verbose && !me) (*mlparams.output)
#define MINPRINT MINPRINT_NOPREFIX << prefix

    init_precision = mlparams.output->precision();
    if(!me) mlparams.output->precision(14);

    if(cost < 0.){
      if(mlparams.exit_on_convergence_fail) error_exit(std::cout << prefix << "Cost " << cost << " cannot be negative!\n");
      else{
	fail = true; //will prevent iteration
      }
    }
  }

  void restoreEnvironment(){
    if(!me) mlparams.output->precision(init_precision);
  }


  void updateAlgorithmState(){
    use_derivs_from_last_step = false; //no point recalculating the derivatives again if we are restoring parameters to what they were before
    
    cost = function.cost(*parameters);   
    if(cost < 0.){
      if(mlparams.exit_on_convergence_fail) error_exit(std::cout << prefix << "Cost " << cost << " cannot be negative!\n");
      else{
	fail = true;
	return;
      }
    }

    MINPRINT << "current cost " << cost << std::endl;
    
    CostType cost_change = cost - last_cost; 
    MINPRINT << "previous cost " << last_cost << std::endl;
    MINPRINT << "d(cost) = "<<cost_change << " d(cost)/tolerance = "<< cost_change/mlparams.delta_cost_min << std::endl;

    if(cost_change >=0.0){
      MINPRINT << "cost increase " << last_cost << " -> " << cost << ", restoring parameters and increasing lambda\n";
      *parameters = *last_params;
      cost = last_cost;
      
      if(lambda >= mlparams.lambda_max){
	if(  fabs(cost_change/cost)<= 1e-12 ){
	  //At this point the fractional change in cost is likely to be heavily influenced by the finite precision of the
	  //SVD matrix inversion. There is no point trying to get closer to the minimum from here, so we may as well assume
	  //that it has been found.
	  MINPRINT << "WARNING: cost increased but at a level within the bounds of the error on the fit algorithm ( delta(cost)/cost <= 1e-12 ). Assuming minima has been found.\n";
	  converged=true;
	  return;
	}else{
	  //This is the end. If we let it continue it will just get stuck in a loop repeatedly trying this step
	  //with ever increasing lambda. If we let lambda continue to grow, the changes in the parameters 
	  //will get smaller and smaller without finding the true minima (unless lambda_max is just too small)

	  if(mlparams.exit_on_convergence_fail){
	    error_exit(std::cout << prefix << "lambda value of " << lambda << " has reached the maximum of " << mlparams.lambda_max << " on iteration " << iter << ", could not find minima\n"); 
	  }else{
	    fail = true;
	    return;
	  }
 
	}
      }else{
	lambda *= mlparams.lambda_factor;
      }
      use_derivs_from_last_step = true; //restoring so use stored derivatives
    }else{
      //test whether algorithm has converged
      if(fabs(cost_change) <= mlparams.delta_cost_min){
	MINPRINT << "change within converge tolerance\n";
	MINPRINT << "Final parameters : " << parameters->print() << std::endl;
	converged=true;
	return;
      }

      //otherwise modify lambda for parameter update
      if(lambda <= mlparams.lambda_min){
	MINPRINT << "lambda has reached its minimum, not reducing\n";
      }else{
	MINPRINT << "cost decrease, reducing lambda\n";
	lambda = lambda/mlparams.lambda_factor;
      }
    }
  }

  //Suggest parameters for next iteration
  void updateParameters(){
    const int nparam = derivs.size(); //requires a size()
     
    M = second_derivs; //requires operator()(int,int) accessor
    c = derivs; //requires operator()(int) accessor
    
    for(int i=0;i<nparam;i++){
      c(i) = -0.5*c(i);

      for(int j=0;j<nparam;j++)
	M(i,j) = 0.5*M(i,j);

      M(i,i) = M(i,i) + 0.5 * lambda*D(i,i);
    }

    MINPRINT << "inverting update matrix\n";
    MINPRINT_NOPREFIX << M.print() << std::endl;
    inv_M = function.invert(M);

    MINPRINT << "inverse matrix : " << std::endl;
    MINPRINT_NOPREFIX << inv_M.print() << std::endl;
    
    delta->zero();
    
    for(int i=0; i<nparam; i++)
      for(int j=0;j<nparam;j++)
	(*delta)(i) = (*delta)(i) + inv_M(i,j)*c(j);

    MINPRINT << "delta : " << delta->print() << std::endl;
    
    *parameters = *parameters + *delta; //must have + operator

    MINPRINT << "Updated parameters: "<< parameters->print() << std::endl;
    ++iter;
  }

  //Use current Hessian to update the dampening matrix
  void updateDampeningMatrix(){
    const int nparam = derivs.size();

    if(iter == 0){
      //Initialize the matrix to zero
      D = second_derivs;
      for(int i=0;i<nparam;i++)
	for(int j=0;j<nparam;j++)
	  D(i,j) = 0.;
    }
    if(mlparams.dampening_matrix == MLdampeningMatrix::Unit && iter == 0){
      for(int i=0;i<nparam;i++) D(i,i) = 1.;
    }else if(mlparams.dampening_matrix == MLdampeningMatrix::HessianDiag){
      for(int i=0;i<nparam;i++) D(i,i) = second_derivs(i,i);
    }else if(mlparams.dampening_matrix == MLdampeningMatrix::MaxHessianDiag){
      if(iter == 0) for(int i=0;i<nparam;i++) D(i,i) = second_derivs(i,i);
      else{
	for(int i=0;i<nparam;i++){
	  auto absii = fabs(second_derivs(i,i));
	  D(i,i) = absii > D(i,i) ? absii : D(i,i);
	}
      }
    }
  }


  void iterate(){ //'params' is the running set of parameters
    MINPRINT << "----------------------------------------------\n";
    MINPRINT << "calculating parameters on iteration " << iter << std::endl;
    if(iter==mlparams.max_iter)
      if(mlparams.exit_on_convergence_fail) error_exit(std::cout << prefix << "reached max iterations\n");
      else{ fail = true; return; }

    MINPRINT << "Input parameters: "<< parameters->print() << std::endl;
    MINPRINT << parameters->size() << " free parameters\n";

    if(iter>0){
      updateAlgorithmState(); //update lambda for next parameter update and check convergence if cost decreased since last iteration
      if(converged || fail) return;
    }

    //Compute cost function derivatives if necessary
    if(use_derivs_from_last_step){
      MINPRINT << "Using derivs from last step\n";
      derivs = last_derivs;
      second_derivs = last_second_derivs;
    }else{   
      MINPRINT << "Computing derivs\n";
      function.derivatives(derivs,second_derivs,*parameters);
    }
    MINPRINT << "derivatives " << derivs.print() << std::endl;
    
    updateDampeningMatrix();

    //Store cost function state as 'old' state prior to update
    *last_params = *parameters;  //input parameters to this cycle
    last_cost = cost; //cost with those inputs
    last_derivs = derivs; //derivatives of cost wrt input params
    last_second_derivs = second_derivs;
    
    MINPRINT << "Lambda = " <<lambda << std::endl;

    updateParameters();
  }
    
public:
  MarquardtLevenbergMinimizer(const CostFunction &func, const AlgorithmParameterType &_mlparams): function(func), iter(0), mlparams(_mlparams){
  }
  
  CostType fit(ParameterType &params){
    initialize(params);
    if(parameters->size() == 0){ //can't fit to zero parameters; just return cost
      converged = true;
    }else{
      while(!converged && !fail)
	iterate();
    }
    restoreEnvironment();
    params = *parameters;
    return cost;
  }
  inline bool hasConverged() const{ return converged; }
  inline int iterations() const{ return iter; }
};
#undef MINPRINT

CPSFIT_END_NAMESPACE

#endif

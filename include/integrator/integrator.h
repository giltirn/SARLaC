#ifndef CPSFIT_INTEGRATOR_H_
#define CPSFIT_INTEGRATOR_H_

//Numerical integration for real numbers

#include<config.h>
#include<utils/macros.h>
#include<cassert>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_gamma.h>
#include<functional>

CPSFIT_START_NAMESPACE

class GSLfunctionWrapperBase{
public:
  virtual double operator()(const double x) const = 0;
  virtual ~GSLfunctionWrapperBase(){}
};

//Wrapper for std::function<double(double)>
class GSLfunctionLambdaWrapper: public GSLfunctionWrapperBase{
  std::function<double(double)> f;
public:
  GSLfunctionLambdaWrapper(const std::function<double(double)> &f): f(f){}
  
  double operator()(const double x) const override{
    return f(x);
  }
};

template<typename Functional>
class GSLfunctionFunctionalWrapper: public GSLfunctionWrapperBase{
  const Functional &f;
public:
  GSLfunctionFunctionalWrapper(const Functional &f): f(f){}
  
  double operator()(const double x) const override{
    return f(x);
  }
};



static double gsl_int_func_wrap(double x, void *f){
  GSLfunctionWrapperBase const &fb = *((GSLfunctionWrapperBase const*)f);
  return fb(x);
}

struct integratorOptions{
  double epsabs; //target absolute error on integral
  double epsrel; //target relative error on integral

  //Note that the algorithm converges when
  //abserr = |result - I| <=  max( epsabs , epsrel * |I| )
  //where I is the integral
  //either can be set to 0 if only one stopping condition is preferred

  int n; //workspace size
  //For workspace size n refer to https://www.gnu.org/software/gsl/doc/html/integration.html#cquad-doubly-adaptive-integration
  

  integratorOptions(){
    epsabs = 0; //don't use absolute error
    epsrel = 1e-8; //target relative error of 1e-8
    n=100;   //default of 100 should be sufficient for most cases according to docs
             //use larger values if the integrator error is too large
  }
};



//Perform the numerical integration of f between a and b using the CQUAD algorithm
//Return a pair of doubles with
//      first: the result
//      second: the absolute error
inline std::pair<double,double> integrate(const GSLfunctionWrapperBase &f, double a, double b, const integratorOptions opts = integratorOptions() ){
    gsl_function gf;
    gf.function = &gsl_int_func_wrap;
    gf.params = (void*)&f;

    //do the integral
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(opts.n);

    double result;
    double abserr;
    assert( gsl_integration_cquad(&gf, a, b, opts.epsabs, opts.epsrel, workspace, &result, &abserr, NULL) == GSL_SUCCESS );
    gsl_integration_cquad_workspace_free(workspace);

    //Stopping conditions of cquad don't match those in the docs, there are circumstances where it can appear to succeed but still produce
    //an error larger than the condition
    if(opts.epsrel > 0. && abserr > fabs(result)*opts.epsrel ){
      error_exit(std::cout << "integrate with result " << result << " failed as abserr(" << abserr << ") > fabs(result)*opts.epsrel (" << fabs(result)*opts.epsrel << ")" );
    }
    if(opts.epsabs > 0. && abserr > opts.epsabs ){
      error_exit(std::cout << "integrate with result " << result << " failed as abserr(" << abserr << ") > opts.epsabs (" << opts.epsabs << ")" );
    }
      
    return std::pair<double,double>(result,abserr);
}

inline std::pair<double,double> integrate(const std::function<double(double)> &f, double a, double b, const integratorOptions opts = integratorOptions() ){
  return integrate(GSLfunctionLambdaWrapper(f), a, b, opts);
}

enum _Infinity { Infinity };

//Perform the semi-infinite integral between a and infinity using the QAGI algorithm
// inline std::pair<double,double> integrate(const GSLfunctionWrapperBase &f, double a, _Infinity inf, const integratorOptions opts = integratorOptions() ){
//     gsl_function gf;
//     gf.function = &gsl_int_func_wrap;
//     gf.params = (void*)&f;

//     //do the integral
//     gsl_integration_workspace * workspace = gsl_integration_workspace_alloc(opts.n);

//     double result;
//     double abserr;
//     gsl_integration_qagiu(&gf, a, opts.epsabs, opts.epsrel, opts.n, workspace, &result, &abserr);
//     gsl_integration_workspace_free(workspace);
//     return std::pair<double,double>(result,abserr);
// }

//Function remapping for infinite upper bound
//https://www.gnu.org/software/gsl/doc/html/integration.html#qagi-adaptive-integration-on-infinite-intervals
class GSLfunctionRemapUpperBoundInfinity: public GSLfunctionWrapperBase{
  GSLfunctionWrapperBase const &f;
  double a; //lower bound
public:
  GSLfunctionRemapUpperBoundInfinity(GSLfunctionWrapperBase const &f, const double a): f(f), a(a){}

  double operator()(const double t) const override{
    return f( a + (1-t)/t ) / t/t;
  }
};

inline std::pair<double,double> integrate(const GSLfunctionWrapperBase &f, double a, _Infinity inf, const integratorOptions opts = integratorOptions() ){
  GSLfunctionRemapUpperBoundInfinity fwrap(f, a);
  return integrate(fwrap, 0., 1., opts);
}





CPSFIT_END_NAMESPACE

#endif

#ifndef _CPSFIT_PVALUE_H_
#define _CPSFIT_PVALUE_H_

//An implementation of the chi^2, F and T^2 distributions, and the concept of p-value

#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <cstdlib>

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>

CPSFIT_START_NAMESPACE

struct chiSquareDistribution{
private:
  static inline double Gamma(const double x){
    gsl_sf_result r;
    int ret = gsl_sf_gamma_e(x, &r);
    if(ret != GSL_SUCCESS) error_exit(std::cout << "gsl_sf_gamma_e failed with error " << gsl_strerror(ret) << std::endl);
    return r.val;
  }
public:
  //PDF for k degrees of freedom and chisq between x and x+dx 
  static inline double PDF(const double k, const double x){
    return 1./pow(2.,k/2)/Gamma(k/2)*::pow(x,k/2-1.)*::exp(-x/2.);
  }

private:

  //Private methods
  struct pdf_params{
    double k;
  };
  static double pdf_func(double x, void *gparams){
    pdf_params *p = (pdf_params*)gparams;
    return PDF(p->k,x);
  };

  
public:
  //p( chi^2 < x; k)  explicit calculation by numerical integral
  static inline double CDF_int(const double k, const double x){
    pdf_params pdf_p;
    pdf_p.k = k;
    
    gsl_function pdf;
    pdf.function = &pdf_func;
    pdf.params = &pdf_p;

    //do the integral
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(100);
    
    double start = 0;
    double end = x;

    double epsabs = 1e-6;
    double epsrel = 1e-6;
    
    double result;
    gsl_integration_cquad(&pdf, start, end, epsabs, epsrel, workspace, &result, NULL,NULL);
    gsl_integration_cquad_workspace_free(workspace);
    return result;
  }

  //Using built-in
  static inline double CDF(const double k, const double x){
    gsl_sf_result r;

    //    These routines compute the normalized incomplete Gamma Function Q(a,x) = 1/\Gamma(a) \int_x^\infty dt t^{a-1} \exp(-t) for a > 0, x >= 0. 
    int ret = gsl_sf_gamma_inc_Q_e (k/2, x/2, &r);
    if(ret != GSL_SUCCESS) error_exit(std::cout << "gsl_sf_gamma_inc_Q_e failed with error " << gsl_strerror(ret) << std::endl);
    return 1.0-r.val;
  }


  //p = \int_x^inf dx P(k, x)     k=dof, x=chi^2
  static inline double pvalue(const double k, const double x){
    return 1.-CDF(k,x); //chi^2 is unit-normalized
  }

  inline static double mean(const double k){ return k; }
  inline static double variance(const double k){ return 2*k; }
};


//F-distribution, also known as Snedecor's F distribution or the Fisher-Snedecor distribution 
struct Fdistribution{
private:
  static inline double Beta(const double a, const double b){
    gsl_sf_result r;
    int ret = gsl_sf_beta_e(a,b, &r);
    if(ret != GSL_SUCCESS) error_exit(std::cout << "gsl_sf_beta_e failed with error " << gsl_strerror(ret) << std::endl);
    return r.val;
  }

public:
  static inline double PDF(const double x, const double d1, const double d2){
    double beta = Beta(d1/2,d2/2);
    double a = pow(d1/d2, d1/2);
    double b = pow(x, d1/2-1);
    double c = pow( 1 + d1*x/d2, -d1/2 - d2/2);

    return a*b*c/beta;

    //return sqrt( pow(d1*x, d1) * pow(d2,d2) / pow(d1*x + d2, d1+d2) ) / x / Beta(d1/2,d2/2);
  }

private:

  //Private methods
  struct pdf_params{
    double d1;
    double d2;
  };
  static double pdf_func(double x, void *gparams){
    pdf_params *p = (pdf_params*)gparams;
    return PDF(x,p->d1,p->d2);
  };

  
public:
  //p( F < x; d1,d2)  explicit calculation by numerical integral
  static inline double CDF_int(const double x, const double d1, const double d2){
    pdf_params pdf_p;
    pdf_p.d1 = d1;
    pdf_p.d2 = d2;
    
    gsl_function pdf;
    pdf.function = &pdf_func;
    pdf.params = &pdf_p;

    //do the integral
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(100);
    
    double start = 0;
    double end = x;

    double epsabs = 1e-6;
    double epsrel = 1e-6;
    
    double result;
    gsl_integration_cquad(&pdf, start, end, epsabs, epsrel, workspace, &result, NULL,NULL);
    gsl_integration_cquad_workspace_free(workspace);
    return result;
  }

  //Using built-in
  static inline double CDF(const double x, const double d1, const double d2){
    gsl_sf_result r;
    int ret = gsl_sf_beta_inc_e(d1/2, d2/2, d1*x/(d1*x + d2), &r);
    if(ret != GSL_SUCCESS) error_exit(std::cout << "gsl_sf_gamma_inc_Q_e failed with error " << gsl_strerror(ret) << std::endl);
    return r.val;
  }   
};
  

//https://en.wikipedia.org/wiki/Hotelling%27s_T-squared_distribution
//Using the sample covariance matrix computed from N samples, the quantity we usually call "chi^2" for p degrees of freedom is actually distributed according to a Hotelling's T^2 distribution T^2(p,N-1)
//distribution with 
struct TsquareDistribution{
public:
  static inline double PDF(const double x, const double p, const double n){
    return (n-p+1)/p/n * Fdistribution::PDF((n-p+1)/p/n*x, p, n-p+1);
  }

private:
  struct pdf_params{
    double p;
    double n;
  };
  static double pdf_func(double x, void *gparams){
    pdf_params *p = (pdf_params*)gparams;
    return PDF(x,p->p,p->n);
  }

public:
  //p( T2 < x; p,n)  explicit calculation by numerical integral
  static inline double CDF_int(const double x, const double p, const double n){
    pdf_params pdf_p;
    pdf_p.p = p;
    pdf_p.n = n;
    
    gsl_function pdf;
    pdf.function = &pdf_func;
    pdf.params = &pdf_p;

    //do the integral
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(100);
    
    double start = 0;
    double end = x;

    double epsabs = 1e-6;
    double epsrel = 1e-6;
    
    double result;
    gsl_integration_cquad(&pdf, start, end, epsabs, epsrel, workspace, &result, NULL,NULL);
    gsl_integration_cquad_workspace_free(workspace);
    return result;
  }

  //p(a <= T2 <= b; p,n)
  static inline double rangeInt(const double a, const double b, const double p, const double n){
    pdf_params pdf_p;
    pdf_p.p = p;
    pdf_p.n = n;
    
    gsl_function pdf;
    pdf.function = &pdf_func;
    pdf.params = &pdf_p;

    //do the integral
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(100);
    
    double start = a;
    double end = b;

    double epsabs = 1e-6;
    double epsrel = 1e-6;
    
    double result;
    gsl_integration_cquad(&pdf, start, end, epsabs, epsrel, workspace, &result, NULL,NULL);
    gsl_integration_cquad_workspace_free(workspace);
    return result;
  }

  static inline double mean(const double p, const double n){
    return n*p/(n-p-1);
  }
  static inline double variance(const double p, const double n){
    return 2*n*n*p*(n-1)/pow(n-p-1,2)/(n-p-3);
  }

  //Using built-in
  static inline double CDF(const double x, const double p, const double n){
    return Fdistribution::CDF((n-p+1)*x/p/n, p, n-p+1);
  }
   
  //p = \int_x^inf dx P(x; p,n)     p=dof, n=#samples + 1,  x="chi^2" 
  static inline double pvalue(const double x, const double p, const double n){
    return 1.-CDF(x,p,n);
  }
};

  
CPSFIT_END_NAMESPACE
#endif

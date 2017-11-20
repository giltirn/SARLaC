#ifndef _CPSFIT_PVALUE_H_
#define _CPSFIT_PVALUE_H_
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_gamma.h>
#include <utils.h>
#include <cstdlib>

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
    return 1./pow(2.,k/2)/Gamma(k/2)*pow(x,k/2-1.)*exp(-x/2.);
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


  //p = \int_x^inf chi^2(k,x)
  static inline double pvalue(const double k, const double x){
    return 1.-CDF(k,x); //chi^2 is unit-normalized
  }
    
};
  
#endif

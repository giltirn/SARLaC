#ifndef _DELTA_0_PHENO
#define _DELTA_0_PHENO

//Phenomenological curve for I=0 and I=2 pipi-scattering phase shift from Schenk, Nucl.Phys. B363 (1991) 97-116

#include<cstdio>
#include<cmath>

#include<utils/macros.h>

CPSFIT_START_NAMESPACE

class PhenoCurve{
public:
  //Calculates \delta_0 pi-pi scattering length in radians
  //s is the Mandelstam variable s=E_pipi^2  *in MeV^2*
  //I is the isospin : 0 or 2
  //mpi in MeV
  //fitcurve is A, B or C - these are the three curves in Schenk, Nucl.Phys. B363 (1991) 97-116
  //                        for I=0 the curves are upper, best and lower bounds respectively, and for I=2 they are lower, best and upper bounds respectively
  static double compute(const double &s, const int &I, const double &mpi, const char fitcurve = 'B'){
    double _4mpi2 = 4*pow(mpi,2);
    double sbrack = (s - _4mpi2)/_4mpi2;
    //printf("(s-4m_pi^2)/(4m_pi^2) = %f\n",sbrack);
    
    int i= (I == 0 ? 0 : 1);
    int fci;
    if(fitcurve == 'A') fci = 0;
    else if(fitcurve == 'B') fci = 1;
    else if(fitcurve == 'C') fci = 2;    
    else{ printf("PhenoCurve::compute invalid fitcurve\n"); exit(-1); }

    //Coefficients for I=0,2 respectively
    double a[2] = {0.2, -0.042};
    double b[2] = {0.24, -0.075};
    double sl[3][2] = {
      {pow(840,2), -pow(1000,2)},
      {pow(865,2), -pow(790,2)}, 
      {pow(890,2), -pow(680,2)}
    }; //A,B,C curves resp, units are MeV
    double c[3][2] = {
      {0.008,0.0},
      {0.0,0.0},
      {-0.015,0.0} 
    };

    double tandelta = sqrt(1. - _4mpi2/s) * ( a[i] + b[i]*sbrack + c[fci][i]*pow(sbrack,2) ) * (_4mpi2 - sl[fci][i])/(s-sl[fci][i]);
    return atan(tandelta);
  }
  static double compute_deriv(const double &s, const int &I, const double &mpi, const char fitcurve = 'B', const double &frac_shift = 1e-04){
    double ds = frac_shift * s;
    return ( compute(s+ds,I,mpi,fitcurve) - compute(s-ds,I,mpi,fitcurve) )/(2.*ds);
  }

};

CPSFIT_END_NAMESPACE
#endif

#ifndef LUSCHER_ZETA_H
#define LUSCHER_ZETA_H
#include<gsl/gsl_integration.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_gamma.h>

#include<config.h>
#include<utils/macros.h>
#include<tensors/gsl_vector.h>

//C.Kelly 2014
//Classes and methods to determine and manipulate Luscher's \zeta and \phi functions for use in the quantization condition and Lellouch-Luscher factor
//This calculaton method is based on Qi Liu's Matlab code which in turn was based on Physical review D 70 074513

//Must link against gsl and gslcblas libraries

CPSFIT_START_NAMESPACE

inline double square(const double x){ return x*x; }

//Luscher's quantization condition is
//n*\pi - \delta_0(k) = \phi(q) 
//where q = kL/2\pi  and  k is the solution of E_{\pi\pi} = 2\sqrt( m_\pi^2 + k^2 )

//Here  \tan(\phi(q)) = -\pi^{3/2}q / \zeta_{00}(1;q^2)
//This class computed \zeta_{00} and by extension, \phi, as a function of q^2, for periodic/H-parity/G-parity BCs in arbitrary spatial directions

class LuscherZeta{
  //Private variables
  GSLvector d; //H-parity or G-parity twist directions (x,y,z): there is a twist, use 1.0 or else use 0.0 
  GSLvector dnorm; //Normalised BC vector
  int N; //First integral is over integer-component vectors in Z^3. The input parameter 'N" sets the maximum magnitude of these vectors

  //Errors on the integrals
  double epsabs;
  double epsrel;

  //Private methods
  struct zeta_params{
    double q2;
    double gn2;
  };
  static double zeta_func(double t, void *zparams){
    const static double pi2 = M_PI*M_PI;
    zeta_params *p = (zeta_params*)zparams;
    return exp( t*p->q2 - pi2 * p->gn2/t ) * pow(M_PI/t, 1.5);
  };

  //If non-null, abserr is the absolute error on the result, and nevals is the number of integration steps used
  double int_zeta(const double q2, const double gn2, double *abserr = NULL, size_t *nevals = NULL) const{
    //Prepare the zeta function
    zeta_params zeta_p;
    zeta_p.q2 = q2;
    zeta_p.gn2 = gn2;

    gsl_function zeta;
    zeta.function = &zeta_func;
    zeta.params = &zeta_p;

    //do the integral
    gsl_integration_cquad_workspace * workspace = gsl_integration_cquad_workspace_alloc(100);
    
    double start = 1e-06;
    double end = 1.0;

    double result;
    gsl_integration_cquad(&zeta, start, end, epsabs, epsrel, workspace, &result, abserr, nevals);
    gsl_integration_cquad_workspace_free(workspace);
    return result;
  }

  double z3sum(const double q, const GSLvector &n) const{
    double r2 = square(n[0]+d[0]/2.0)+square(n[1]+d[1]/2.0)+square(n[2]+d[2]/2.0);
    double q2 = square(q);
    double out = 0.0;
        
    //printf("z3sum with q = %f and n = %f %f %f\n",q,n[0],n[1],n[2]);

    //sum over e^{-(r^2 - q^2)}/(r^2 - q^2)        
    out += exp(q2-r2) / (r2-q2);

    //skip integral for 0 vector
    if(n[0]==0.0 && n[1]==0.0 && n[2]==0.0)
      return out;

    //integral \int_0^1 dt e^{q^2t} (\pi/t)^{3/2} e^{-\pi^2 r^2/t} modified for APBC where appropriate (modifies r (aka n) )
    //define n_perp = n - (n \cdot dnorm) dnorm  is perpendicular to d
    double dcoeff = dot(n,dnorm);
    GSLvector np = n - dcoeff*dnorm;
    double gn2 = square(dcoeff) + np.norm2();

    double int_n = int_zeta(q*q,gn2);
    out += int_n*pow(-1, dot(n,d));

    return out;
  }
  
 public:
  LuscherZeta(): d(3), dnorm(3), N(5), epsabs(1e-06),epsrel(1e-06) {}
  LuscherZeta(const double x, const double y, const double z): d(3), dnorm(3), N(5), epsabs(1e-06),epsrel(1e-06){
    setTwists(x,y,z);
  }

  void setTwists(const double x, const double y, const double z){
    d[0] = x; d[1] = y; d[2] = z;
    double dnrm = d.norm();
    if(d[0]==0 && d[1]== 0 && d[2]==0) dnrm = 1;
    
    for(int i=0;i<3;i++){
      if(d[i] != 0.0 && d[i] != 1.0){ std::cout << "LuscherZeta::setTwists : Error, arguments must be 0 or 1\n"; exit(-1); }
      dnorm[i] = double(d[i])/dnrm; //normalize d
    }
  }

  inline void setIntegrationErrorBounds(const double eps_abs, const double eps_rel){
    epsabs = eps_abs; epsrel = eps_rel;
  }
     
  inline void setMaximumVectorMagnitude(const int iN){
    N = iN;
  }
        
  double calcZeta00(const double q) const{
    //q is a scalar modulus
    double result = 0.0;
    //outer loop over vectors in Z^3
    for(int nx = -N; nx <= N; ++nx){
      int Ny( floor(sqrt( double(N*N - nx*nx) )) );
      for(int ny = -Ny; ny <= Ny; ++ny){
	int Nz( floor(sqrt( double(N*N - nx*nx - ny*ny) ) + 0.5 ) ); //round to nearest int
	for(int nz = -Nz; nz <= Nz; ++nz){
	  GSLvector n(3);
	  n[0] = nx; n[1] = ny; n[2] = nz;
	  result += z3sum(q,n);
	}
      }
    }
    
    //constant part
    double const_part = 0.0;
    double q2 = square(q);
    bool warn = true;
    for(unsigned int l=0;l<100;l++){
      double c_aux = pow(q2,l) / gsl_sf_fact(l) / (l-0.5);
      if(l>9 && c_aux < 1e-08 * fabs(const_part)){
	warn = false;
	break;
      }
      const_part += c_aux;
    }
    if(warn) printf("Warning: reaches the maximum loop number when doing the summation\n");

    result += pow(M_PI,1.5)*const_part;
    result /= sqrt(4*M_PI);
    return result;
  }

  inline double calcPhi(const double q) const{
    return atan(-q*pow(M_PI,1.5)/calcZeta00(q));
  }

  inline double calcPhiDeriv(const double q, const double frac_shift = 1e-04) const{
    double dq = frac_shift * q;
    return ( calcPhi(q+dq) - calcPhi(q-dq) )/(2.*dq);
  }
};

//Test by comparing with numbers computed by Daiqian (based off Qi's MatLab code) 
class ZetaTesting{
public:
  ZetaTesting(){
    LuscherZeta zeta;
    
    //q  dx dy dz phi
    //0.100000	0	0	1	0.790456
    zeta.setTwists(0,0,1);
    printf("0.100000	0	0	1	Expect 0.790456 Got %.6f\n", zeta.calcPhi(0.1));
      
    //0.300000	1	1	1	1.079671
    zeta.setTwists(1,1,1);
    printf("0.300000	1	1	1	Expect 1.079671 Got %.6f\n", zeta.calcPhi(0.3));

    //0.400000	0	0	0	0.571791
    zeta.setTwists(0,0,0);
    printf("0.400000	0	0	0	Expect 0.571791  Got %.6f\n", zeta.calcPhi(0.4));

    //2.000000	1	1	1	-1.315088
    zeta.setTwists(1,1,1);
    printf("2.000000	1	1	1	Expect -1.315088 Got %.6f\n", zeta.calcPhi(2.0));
  }
};

//Compute the pipi-scattering phase shift in degrees given a pipi-energy E, a pion mass m, lattice spatial size L (assumed equal for all 3 spatial directions)
//and a zeta function/number of antiperiodic directions

inline double phaseShiftZ(const double E, const double m, const int L, const LuscherZeta &zeta){
  double k = sqrt( pow(E/2.,2.0) - m*m );
  double q = k * L/2./M_PI;
  double delta = -zeta.calcPhi(q);
  while(delta > M_PI) delta -= M_PI;
  while(delta < -M_PI) delta += M_PI;
  return delta/M_PI * 180;
}

inline double phaseShift(const double E, const double m, const int L, const std::vector<int> twists = {0,0,0}){
  LuscherZeta zeta(twists[0],twists[1],twists[2]);
  return phaseShiftZ(E,m,L,zeta);
}

CPSFIT_END_NAMESPACE

#endif


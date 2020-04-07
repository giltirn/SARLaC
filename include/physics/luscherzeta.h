#ifndef LUSCHER_ZETA_H
#define LUSCHER_ZETA_H
#include<array>
#include<cassert>
#include<gsl/gsl_integration.h>
#include<gsl/gsl_math.h>
#include<gsl/gsl_sf_gamma.h>

#include<config.h>
#include<utils/macros.h>
#include<tensors/gsl_vector.h>

//C.Kelly 2014
//Classes and methods to determine and manipulate Luscher's \zeta and \phi functions for use in the quantization condition and Lellouch-Luscher factor
//This calculaton method is based on Qi Liu's Matlab code which in turn was based on Physical review D 70 074513  (https://arxiv.org/pdf/hep-lat/0402025.pdf)

//Must link against gsl and gslcblas libraries

CPSFIT_START_NAMESPACE

//Luscher's quantization condition is
//n*\pi - \delta_0(k) = \phi(q) 
//where q = kL/2\pi  and  k is the solution of E_{\pi\pi} = 2\sqrt( m_\pi^2 + k^2 )

//Here  \tan(\phi(q)) = -\pi^{3/2}q / \zeta_{00}(1;q^2)
//This class computed \zeta_{00} and by extension, \phi, as a function of q^2, for periodic/H-parity/G-parity BCs in arbitrary spatial directions

class LuscherZeta{
  //Private variables
  GSLvector l; //H-parity or G-parity twist directions (x,y,z): there is a twist, use 1.0 or else use 0.0 

  GSLvector d; //(L/(2\pi))\vec P_CM     center of mass momentum in units of 2pi/L
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
    
    double start = 1e-12;
    double end = 1.0;

    double result;
    gsl_integration_cquad(&zeta, start, end, epsabs, epsrel, workspace, &result, abserr, nevals);
    gsl_integration_cquad_workspace_free(workspace);
    return result;
  }

  inline GSLvector gammaOpGen(const GSLvector &n, double gamma, bool inv) const{
    GSLvector dunit = d.norm() == 0. ? GSLvector(3,0.) : d/d.norm();
    GSLvector npar = dot(n,dunit)*dunit;
    GSLvector nperp = n - npar;
    return (inv ? 1./gamma : gamma) * npar + nperp;
  }
  inline GSLvector gammaOp(const GSLvector &n, double gamma) const{ return gammaOpGen(n,gamma,false); }
  inline GSLvector invGammaOp(const GSLvector &n, double gamma) const{ return gammaOpGen(n,gamma,true); }

  double z3sum(const double q, const GSLvector &n, double gamma) const{
    GSLvector r = n + d/2. + l/2.;
    r = invGammaOp(r,gamma);

    double r2 = r.norm2();
    double q2 = q*q;
    double out = 0.0;
        
    //sum over e^{-(r^2 - q^2)}/(r^2 - q^2)        
    out += exp(q2-r2) / (r2-q2);

    //skip integral for 0 vector
    if(n[0]==0.0 && n[1]==0.0 && n[2]==0.0)
      return out;
    
    //integral \gamma \int_0^1 dt e^{q^2t} (\pi/t)^{3/2} e^{-\pi^2 |\vec h|^2/t}
    //where   \vec h = \hat gamma \vec n
    GSLvector h = gammaOp(n,gamma);
    double h2 = h.norm2();

    double int_n = int_zeta(q2,h2);
    out += gamma*int_n*pow(-1, dot(n,d+l));

    return out;
  }

  void setPCMvec(const std::array<double,3> &_d){
    for(int i=0;i<3;i++) d[i] = _d[i];
  }

  //H-parity/G-parity spatial directions: 1,   otherwise 0
  void setTwists(const std::array<int,3> &_l){
    for(int i=0;i<3;i++){
      if(l[i] != 0 && l[i] != 1){ std::cout << "LuscherZeta::setTwists : Error, arguments must be 0 or 1\n"; exit(-1); }
      l[i] = _l[i];
    }
  }
  
 public:
  //\vec d = \vec P_CM / (2pi/L)
  //\vec l :  H-parity/G-parity spatial directions: 1,   otherwise 0
  LuscherZeta(const std::array<int,3> &_l, const std::array<double,3> &_d): l(3), d(3), N(5), epsabs(1e-06),epsrel(1e-06){
    setTwists(_l);
    setPCMvec(_d);
  }

  inline void setIntegrationErrorBounds(const double eps_abs, const double eps_rel){
    epsabs = eps_abs; epsrel = eps_rel;
  }
     
  inline void setMaximumVectorMagnitude(const int iN){
    N = iN;
  }
        
  //gamma: Lorentz boost computed as   E/sqrt(E^2 - |P_CM|^2)   where E is the measured pipi energy
  double calcZeta00(const double q, const double gamma) const{
    //q is a scalar modulus
    double result = 0.0;
    //outer loop over vectors in Z^3
    GSLvector n(3);
    for(int nx = -N; nx <= N; ++nx){
      int Ny( floor(sqrt( double(N*N - nx*nx) )) );
      for(int ny = -Ny; ny <= Ny; ++ny){
	int Nz( floor(sqrt( double(N*N - nx*nx - ny*ny) ) + 0.5 ) ); //round to nearest int
	for(int nz = -Nz; nz <= Nz; ++nz){
	  n[0] = nx; n[1] = ny; n[2] = nz;
	  result += z3sum(q,n,gamma);
	}
      }
    }
    
    //constant part
    double const_part = 0.0;
    double q2 = q*q;
    bool warn = true;
    for(unsigned int l=0;l<100;l++){
      double c_aux = pow(q2,l) / gsl_sf_fact(l) / (l-0.5);
      if(l>9 && c_aux < 1e-08 * fabs(const_part)){
	warn = false;
	break;
      }
      const_part += c_aux;
    }
    if(warn) printf("LuscherZeta warning: In determination reached the maximum loop number when doing the summation\n");

    result += gamma*pow(M_PI,1.5)*const_part;
    result /= sqrt(4*M_PI);
    return result;
  }

  inline double calcPhi(const double q, const double gamma = 1.) const{
    return atan(-gamma*q*pow(M_PI,1.5)/calcZeta00(q,gamma));
  }
  //Assumes gamma = 1.
  inline double calcPhiDeriv(const double q, const double frac_shift = 1e-04) const{
    double dq = frac_shift * q;
    return ( calcPhi(q+dq,1.) - calcPhi(q-dq,1.) )/(2.*dq);
  }

  const GSLvector &getd() const{ return d; }
  const GSLvector &getl() const{ return l; }
};

//Compute the pipi-scattering phase shift in degrees given a pipi-energy E, a pion mass m, lattice spatial size L (assumed equal for all 3 spatial directions)
//and a zeta function/number of antiperiodic directions
//As twists and d not dependent on energy (i.e. don't vary between jackknife samples) we can use the same Zeta instance under a thread loop
enum class DispersionRelation { Continuum, Sin, SinhSin, Fit };

inline double getSinP2(const GSLvector &p){
  double p2 = 0.;
  for(int i=0;i<3;i++){
    double k = 2*sin(p[i]/2.);
    p2 += k*k;
  }
  return p2;
}
inline double getPn(const GSLvector &p,const int n){
 double pp = 0.;
 for(int i=0;i<3;i++)
    pp += pow(p[i],n);
 return pp;
}

//Get energy from mass and momentum
double dispersionRelationGetEnergy(const DispersionRelation disp, double m, GSLvector p, const double punit = 1.){
  for(int i=0;i<3;i++) p[i] *= punit;

  if(disp == DispersionRelation::Continuum){ //E^2 = m^2 + \sum_i k_i^2 
    return sqrt( m*m + p.norm2() );
  }else if(disp == DispersionRelation::Sin){ //E^2 = m^2 + \sum_i [ 2 sin(k_i/2) ]^2
    double p2 = getSinP2(p);
    return sqrt( m*m + p2 );
  }else if(disp == DispersionRelation::SinhSin){ // 4 sinh^2(E/2) = \sum_i [ 2 sin(k_i/2) ]^2  +  m^2  https://physik.uni-graz.at/~axm/lattice2017.pdf (p16)
    double p2 = getSinP2(p);
    double h = p2 + m*m;
    h = sqrt(h/4.);
    return 2*asinh(h);
  }else if(disp == DispersionRelation::Fit){
    return sqrt( m*m + 0.995467 * getPn(p,2)  -0.144335 * getPn(p,4) );
  }else{
    assert(0);
  }
}
//Get mass from energy and momentum
double dispersionRelationGetMass(const DispersionRelation disp, double E, GSLvector p, const double punit = 1.){
  for(int i=0;i<3;i++) p[i] *= punit;

  if(disp == DispersionRelation::Continuum){ //E^2 = m^2 + \sum_i k_i^2 
    return sqrt( E*E - p.norm2() );
  }else if(disp == DispersionRelation::Sin){ //E^2 = m^2 + \sum_i [ 2 sin(k_i/2) ]^2
    double p2 = getSinP2(p);
    return sqrt( E*E - p2 );
  }else if(disp == DispersionRelation::SinhSin){ // 4 sinh^2(E/2) = \sum_i [ 2 sin(k_i/2) ]^2  +  m^2   https://physik.uni-graz.at/~axm/lattice2017.pdf (p16)
    double p2 = getSinP2(p);
    double m2 = pow(2*sinh(E/2), 2) - p2;
    return sqrt(m2);
  }else if(disp == DispersionRelation::Fit){
    return sqrt( E*E - 0.995467 * getPn(p,2)  +0.144335 * getPn(p,4) );
  }else{
    assert(0);
  }
}

double dispersionRelationGetEnergyFromEnergy(const DispersionRelation disp, GSLvector pout, double Ein, GSLvector pin, const double punit = 1.){
  for(int i=0;i<3;i++){ pin[i] *= punit; pout[i] *= punit; }

  if(disp == DispersionRelation::Continuum){ //E2^2 = E1^2 - [\sum_i k1_i^2] + \sum_i k2_i^2  
    return sqrt( Ein*Ein - pin.norm2() + pout.norm2() );
  }else if(disp == DispersionRelation::Sin){ //E^2 = E1^2 -  \sum_i [ 2 sin(k1_i/2) ]^2   + \sum_i [ 2 sin(k2_i/2) ]^2
    double pin2 = getSinP2(pin);
    double pout2 = getSinP2(pout);
    return sqrt( Ein*Ein - pin2 + pout2 );
  }else if(disp == DispersionRelation::SinhSin){ // 4 sinh^2(E2/2) =  4 sinh^2(E1/2) - \sum_i [ 2 sin(k1_i/2) ]^2 + \sum_i [ 2 sin(k2_i/2) ]^2   https://physik.uni-graz.at/~axm/lattice2017.pdf (p16)
    double pin2 = getSinP2(pin);
    double pout2 = getSinP2(pout);
    double h = pow(2*sinh(Ein/2), 2) - pin2 + pout2;
    return 2*asinh( sqrt(h/4.) );
  }else if(disp == DispersionRelation::Fit){
    return sqrt( Ein*Ein + 0.995467 * ( getPn(pout,2) - getPn(pin,2) )  -0.144335 * ( getPn(pout,4) - getPn(pin,4) ) );
  }else{
    assert(0);
  }
}

//Energy E is in the *lab frame*
//Phase shift is in degrees
inline double phaseShiftZ(const LuscherZeta &zeta, double E, const double m, const double L, const DispersionRelation dispn = DispersionRelation::Continuum, const double zero_tol =1e-10){
  const GSLvector &d = zeta.getd();

  //Get CofM energy
  double E_CM = dispersionRelationGetEnergyFromEnergy(dispn, GSLvector(3,.0), E, d, 2*M_PI/L);

  //Compute the Lorentz factor gamma = E/E_CM
  double gamma = E/E_CM;
  
  double k2 = pow(E_CM/2.,2.0) - m*m; //this should include the dispersion reln! FIXME
  if(fabs(k2) < zero_tol ) k2 = 0.;

  double k = sqrt(k2);

  double q = k * L/2./M_PI;
  double delta = -zeta.calcPhi(q,gamma);
  
  while(delta > M_PI) delta -= M_PI;
  while(delta < -M_PI) delta += M_PI;
  return delta/M_PI * 180;
}

//d= Pcm /(2pi/L) 
inline double phaseShift(const double E, const double m, const double L, const std::array<int,3> &twists = {0,0,0}, const std::array<double,3> &d = {0.,0.,0.}, const DispersionRelation dispn = DispersionRelation::Continuum){
  LuscherZeta zeta(twists, d);
  return phaseShiftZ(zeta,E,m,L,dispn);
}

CPSFIT_END_NAMESPACE

#endif


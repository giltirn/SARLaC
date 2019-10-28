#ifndef DELTA_S_1_WILSON_COEFFS_H_
#define DELTA_S_1_WILSON_COEFFS_H_

#include<config.h>
#include<utils/macros.h>
#include<physics/deltaS1_anomalous_dimension.h>

CPSFIT_START_NAMESPACE

//Equation number references are for Buchalla et al https://arxiv.org/pdf/hep-ph/9512380.pdf

inline BasicSquareMatrix<ValOrd> truncateO1e(const BasicSquareMatrix<ValOrd> &m){
  return BasicSquareMatrix<ValOrd>(m.size(), [&](const int i, const int j){ return truncateO1e(m(i,j)); });
}
inline NumericVector<ValOrd> truncateO1e(const NumericVector<ValOrd> &m){
  return NumericVector<ValOrd>(m.size(), [&](const int i){ return truncateO1e(m(i)); });
}
inline NumericVector<ValOrd> truncateO1(const NumericVector<ValOrd> &m){
  return NumericVector<ValOrd>(m.size(), [&](const int i){ return truncateO1(m(i)); });
}

struct DeltaS1WilsonCoeffs{
  typedef NumericVector<double> VectorD;
  typedef BasicSquareMatrix<double> MatrixD;
  typedef NumericVector<ValOrd> VectorO;
  typedef BasicSquareMatrix<ValOrd> MatrixO;

  static double B0(const double x){ return 1./4*(x/(1-x)+x*log(x)/pow(x-1,2)); } //VII.13
  static double C0(const double x){ return x/8*( (x-6)/(x-1)+(3*x+2)*log(x)/pow(x-1,2) ); } //VI.14
  static double D0(const double x){ return -4./9*log(x)  +  (-19*pow(x,3)+25*pow(x,2) )/36/pow(x-1,3)  +  pow(x,2)*(5*pow(x,2)-2*x-6)/18/pow(x-1,4) * log(x); } //VII.15
  static double Dt0(const double xt){ return D0(xt) - 4./9; } //VII.16
  static double E0(const double x){ return -2./3*log(x) + x*(18-11*x-pow(x,2) )/12/pow(1-x,3) + pow(x,2)*(15-16*x+4*pow(x,2) )/6/pow(1-x,4)*log(x); } //VI.15
  static double Et0(const double xt){ return E0(xt) - 2./3; } //VI.16

  //Wilson coefficients at MW for \vec v (cf VI.7)  via VII.3 - VII.12
  //\vec v(\mu) = U(mu, MW)\vec C(MW)    where U includes running and cross-threshold matching
  //xt = mt^2/MW^2  (VI.17)
  //asMW = alpha_s(MW)  (5 flavor)
  static VectorO mCatMW(const double asMW, const double xt, const double a, const double ThetaW){
    double sin2ThetaW = pow(sin(ThetaW),2);

    VectorO C(10);
    C(0) = 11./2 * _AS * asMW/(4 *M_PI);
    C(1) = 1 -11./6 * _AS * asMW/(4 * M_PI) - 35./18 * _AE * a/(4 * M_PI);
    C(2) = -asMW/(24 * M_PI) * Et0(xt) * _AS + a/(6 * M_PI)/sin2ThetaW * (2*B0(xt)+C0(xt)) * _AE;
    C(3) = asMW/(8 * M_PI) * Et0(xt) * _AS;
    C(4) = -asMW/(24 * M_PI) * Et0(xt) * _AS;
    C(5) = asMW/(8 * M_PI) * Et0(xt) * _AS;
    C(6) = a/(6 * M_PI) * (4 * C0(xt)+Dt0(xt)) * _AE;
    C(7) = 0;
    C(8) = a/(6 * M_PI) * (4 * C0(xt)+Dt0(xt)+(10*B0(xt)-4*C0(xt) )/sin2ThetaW) * _AE;
    C(9) = 0;

    return C;
  }

  //\vec v in the 4 flavor theory at mu
  //Christoph truncates to O(a^1, as^0) for the intermediate running
  static VectorO mCat4fMu(const double asMu, const double asMB, const double asMW,
			  const double xt, const double a, const double ThetaW, const int Nc, bool truncate_im = true){
    VectorO CMW = mCatMW(asMW, xt, a, ThetaW);
    MatrixO U5mbMW = DeltaS1anomalousDimension::computeUQCDQED(asMB, asMW, 2, 3, Nc, a); //5f = d=3 u=2
    MatrixO M45 = DeltaS1anomalousDimension::MatMB(asMB, a);
    MatrixO U4mumb = DeltaS1anomalousDimension::computeUQCDQED(asMu, asMB, 2, 2, Nc, a); //4f = d=2 u=2
    VectorO f;
    if(truncate_im){
      f = truncateO1(U4mumb * truncateO1e(M45 * truncateO1e(U5mbMW * CMW)));
    }else{
      f = U4mumb * M45 * U5mbMW * CMW;
    }
    return f;
  }
  //In 3f theory *at mc*
  static VectorO mCat3fMc(const VectorO &mCat4fMc_val, const double asMC, const double a, bool truncate_im = true){
    MatrixO M34 = DeltaS1anomalousDimension::MatMC(asMC, a);
    VectorO f = M34 * mCat4fMc_val;
    if(truncate_im)
      f = truncateO1(f);
    return f;
  }
  //In 3f theory at other scale
  static VectorO mCat3fMu(const VectorO &mCat3fMc_val, const MatrixO &U3f_MCtoMu, bool truncate_im = true){
    VectorO f = U3f_MCtoMu * mCat3fMc_val;
    if(truncate_im)
      f = truncateO1(f);
    return f;
  }


  //For \vec z everything is defined in terms of z_1 and z_2 at the charm mass (see VII.17)
  //z_1(mc) and z_2(mc) are the first 2 elements of the evolved \vec C in the 4-flavor theory at the charm mass
  //We have to evolve using only 2x2 evolution matrices
  static VectorO z1z2atmc(const double asMC, const double asMB, const double asMW, 
			 const double xt, const double a, const double ThetaW, const int Nc){
    VectorO CMW = mCatMW(asMW, xt, a, ThetaW).subvector(0,2);
    MatrixO U5mbMW = DeltaS1anomalousDimension::computeUQCDQED(asMB, asMW, 2, 3, Nc, a, true); //5f = d=3 u=2
    MatrixO M45 = DeltaS1anomalousDimension::MatMB(asMB, a).submatrix(0,0,2);
    MatrixO U4mcmb = DeltaS1anomalousDimension::computeUQCDQED(asMC, asMB, 2, 2, Nc, a, true); //4f = d=2 u=2

    VectorO f = U4mcmb * M45 * U5mbMW * CMW;

    return f;
  }
  //Used in 3f theory I believe
  static VectorO zatmc(const VectorO &z1z2atmc_val, const double asMC, const double a, const bool truncate_im = true){ //VII.17
    const auto &z1 = z1z2atmc_val(0);
    const auto &z2 = z1z2atmc_val(1);

    VectorO out(10);
    out(0) = z1;
    out(1) = z2;
    
    auto Fe= -4./9*(3*z1+z2);
    auto Fs= -2./3 * z2;

    out(2) = -_AS * asMC/(24*M_PI) *  Fs;
    out(3) = _AS * asMC/(8*M_PI) * Fs;
    out(4) = -_AS * asMC/(24*M_PI) * Fs;
    out(5) = _AS * asMC/(8*M_PI) * Fs;
    out(6) = _AE * a/(6*M_PI) * Fe;
    out(7) = 0;
    out(8) = _AE * a/(6*M_PI) * Fe;
    out(9) = 0;

    if(truncate_im) out = truncateO1(out);
    return out;
  }
  static VectorO zat3fmu(const VectorO &zatmc_val, const MatrixO &U3f_MCtoMu, bool truncate_im = true){
    VectorO f = U3f_MCtoMu * zatmc_val;
    if(truncate_im)
      f = truncateO1(f);
    return f;
  }


  static VectorO yatmu(const VectorO zatmu_val, const VectorO &mCat3fMu_val){
    return mCat3fMu_val - zatmu_val;
  }


};


CPSFIT_END_NAMESPACE

#endif

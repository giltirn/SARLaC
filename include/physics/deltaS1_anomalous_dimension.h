#ifndef DELTA_S_1_ANOMALOUS_DIMENSION_H_
#define DELTA_S_1_ANOMALOUS_DIMENSION_H_

#include<config.h>
#include<utils/macros.h>

#include<tensors/basic_square_matrix.h>
#include<physics/valord.h>

CPSFIT_START_NAMESPACE

void cleanup(BasicSquareMatrix<double> &m, const double tol = 1e-12){
  for(int i=0;i<m.size();i++)
    for(int j=0;j<m.size();j++) 
      if(fabs(m(i,j)) < tol) m(i,j) = 0.;
}


//2-loop anomalous dimension
struct DeltaS1anomalousDimension{
  typedef BasicSquareMatrix<double> MatrixD;
  typedef BasicSquareMatrix<ValOrd> MatrixO;

  //Nc = #colors
  //f = #flavors
  //u = #up-type quarks
  //d = #down-type quarks
  
  //\gamma_s^(0)  LO . QCD part from Eq VI.25 of https://arxiv.org/pdf/hep-ph/9512380.pdf
  //QED part from Table XIV
  static MatrixD gammas0(const double Nc, const double f, const double u, const double d){
    return MatrixD(
		  {  -6./Nc, 6., 0, 0, 0, 0, 0, 0, 0, 0,
		      6., -6./Nc, -2./3/Nc, 2./3, -2./3/Nc, 2./3, 0, 0, 0, 0,
		      0, 0, -22./3/Nc, 22./3, -4./3/Nc, 4./3, 0, 0, 0, 0,
		      0, 0, 6 - 2*f/3/Nc, -6./Nc + 2*f/3, -2*f/3/Nc, 2*f/3, 0, 0, 0, 0,
		      0, 0, 0, 0, 6./Nc, -6., 0, 0, 0, 0,
		      0, 0, -2*f/3/Nc, 2*f/3, -2*f/3/Nc, -6*(-1 + Nc*Nc)/Nc + 2*f/3, 0, 0, 0, 0,

		      0, 0, 0, 0, 0, 0, 6./Nc, -6., 0, 0,
		      0, 0, -2*(u - d/2)/3/Nc, 2*(u - d/2)/3, -2*(u - d/2)/3/Nc, 2*(u - d/2)/3, 0, -6*(-1 + Nc*Nc)/Nc, 0, 0,
		      0, 0, 2./3/Nc, -2./3, 2./3/Nc, -2./3, 0, 0, -6./Nc, 6.,
		      0, 0, -2*(u - d/2)/3/Nc, 2*(u - d/2)/3, -2*(u - d/2)/3/Nc, 2*(u - d/2)/3, 0, 0, 6., -6./Nc
		      });
  };
  //\gamma_s^(1)  NLO. QCD part from Eq. VI.26 of https://arxiv.org/pdf/hep-ph/9512380.pdf
  //QED part from Table XV
  static MatrixD gammas1(const double Nc, const double f, const double u, const double d){
    assert(Nc == 3);
    return MatrixD(
		  {-21./2-2*f/9,7./2+2./3*f,79./9,-7./3,-65./9,-7./3,0,0,0,0,
		      7./2+2*f/3,-21./2-2*f/9,-202./243,1354./81,-1192./243,904./81,0,0,0,0,
		      0,0,-5911./486+71./9*f,5983./162+f/3,-2384./243-71*f/9,1808./81-f/3,0,0,0,0,
		      0,0,379./18+56*f/243,-91./6+808*f/81,-130./9-502*f/243,-14./3+646*f/81,0,0,0,0,
		      0,0,-61*f/9,-11*f/3,71./3+61*f/9,-99+11*f/3,0,0,0,0,
		      0,0,-682*f/243,106*f/81,-225./2+1676*f/243,-1343./6+1348*f/81,0,0,0,0,

		      0,0,-61*(u-d/2)/9,-11*(u-d/2)/3,83*(u-d/2)/9,-11*(u-d/2)/3,71./3-22*f/9,-99+22*f/3,0,0,
		      0,0,-682*(u-d/2)/243,106*(u-d/2)/81,704*(u-d/2)/243,736*(u-d/2)/81,-225./2+4*f,-1343./6+68*f/9,0,0,
		      0,0,202./243+73*(u-d/2)/9,-1354./81-(u-d/2)/3,1192./243-71*(u-d/2)/9,-904./81-(u-d/2)/3,0,0,-21./2-2*f/9,7./2+2*f/3,
		      0,0,-79./9-106*(u-d/2)/243,7./3+826*(u-d/2)/81,65./9-502*(u-d/2)/243,7./3+646*(u-d/2)/81,0,0,7./2+2*f/3,-21./2-2*f/9
		      });
  }
  
  //LO EM corrections from Table XVI of https://arxiv.org/pdf/hep-ph/9512380.pdf
  static MatrixD gammae0(const double Nc, const double u, const double d){
    return MatrixD(
		  {-8./3,0,0,0,0,0,16*Nc/27,0,16*Nc/27,0,
		      0,-8./3,0,0,0,0,16./27,0,16./27,0,
		      0,0,0,0,0,0,-16./27+16*Nc*(u-d/2)/27,0,-88./27+16*Nc*(u-d/2)/27,0,
		      0,0,0,0,0,0,-16*Nc/27+16*(u-d/2)/27,0,-16*Nc/27+16*(u-d/2)/27,-8./3,
		      0,0,0,0,0,0,8./3+16*Nc*(u-d/2)/27,0,16*Nc*(u-d/2)/27,0,
		      0,0,0,0,0,0,16*(u-d/2)/27,8./3,16*(u-d/2)/27,0,
		      0,0,0,0,4./3,0,4./3+16*Nc*(u+d/4)/27,0,16*Nc*(u+d/4)/27,0,
		      0,0,0,0,0,4./3,16*(u+d/4)/27,4./3,16*(u+d/4)/27,0,
		      0,0,-4./3,0,0,0,8./27+16*Nc*(u+d/4)/27,0,-28./27+16*Nc*(u+d/4)/27,0,
		      0,0,0,-4./3,0,0,8*Nc/27+16*(u+d/4)/27,0,8*Nc/27+16*(u+d/4)/27,-4./3
		      });
  }

  //NLO EM corrections from Table XVII of https://arxiv.org/pdf/hep-ph/9512380.pdf
  static MatrixD gammase1(const double Nc, const double u, const double d){
    assert(Nc == 3);
    return MatrixD(
		   {194./9,-2./3,-88./243,88./81,-88./243,88./81,152./27,40./9,136./27,56./9,
		    25./3,-49./9,-556./729,556./243,-556./729,556./243,-484./729,-124./27,-3148./729,172./27,
		    0,0,1690./729-136*(u-d/2)/243,-1690./243+136*(u-d/2)/81,232./729-136*(u-d/2)/243,-232./243+136*(u-d/2)/81,3136./729+104*(u-d/2)/27,64./27+88*(u-d/2)/9,20272./729+184*(u-d/2)/27,-112./27+8*(u-d/2)/9,
		    0,0,-641./243-388*u/729+32*d/729,-655./81+388*u/243-32*d/243,88./243-388*u/729+32*d/729,-88./81+388*u/243-32*d/243,-152./27+3140*u/729+656*d/729,-40./9-100*u/27-16*d/27,170./27+908*u/729+1232*d/729,-14./3+148*u/27-80*d/27,
		    0,0,-136./243*(u-d/2),136./81*(u-d/2),-2-136./243*(u-d/2),6+136./81*(u-d/2),-232./9+104./27*(u-d/2),40./3+88./9*(u-d/2),184*(u-d/2)/27,8./9*(u-d/2),
		    0,0,-748*u/729+212*d/729,748*u/243-212*d/243,3-748*u/729+212*d/729,7+748*u/243-212*d/243,-2-5212*u/729+4832*d/729,182./9+188*u/27-160*d/27,-2260*u/729+2816*d/729,-140*u/27+64*d/27,
		    0,0,-136*(u+d/4)/243,136*(u+d/4)/81,-116./9-136*(u+d/4)/243,20./3+136*(u+d/4)/81,-134./9+104*(u+d/4)/27,38./3+88*(u+d/4)/9,184*(u+d/4)/27,8*(u+d/4)/9,
		    0,0,-748*u/729-106*d/729,748*u/243+106*d/243,-1-748*u/729-106*d/729,91./9+748*u/243+106*d/243,2-5212*u/729-2416*d/729,154./9+188*u/27+80*d/27,-2260*u/729-1408*d/729,-140*u/27-32*d/27,
		    0,0,7012./729-136*(u+d/4)/243,764./243+136*(u+d/4)/81,-116./729-136*(u+d/4)/243,116./243+136*(u+d/4)/81,-1568./729+104*(u+d/4)/27,-32./27+88*(u+d/4)/9,5578./729+184*(u+d/4)/27,38./27+8*(u+d/4)/9,
		       0,0,1333./243-388*u/729-16*d/729,107./81+388*u/243+16*d/243,-44./243-388*u/729-16*d/729,44./81+388*u/243+16*d/243,76./27+3140*u/729-328*d/729,20./9-100*u/27+8*d/27,140./27+908*u/729-616*d/729,-28./9+148*u/27+40*d/27});


  }

  inline static double computeBeta0(const int Nf, const int Nc = 3){
    return double(11.*Nc - 2*Nf)/3.;
  }
  inline static double computeBeta1(const int Nf, const int Nc = 3){
    double CF = double(Nc*Nc - 1)/2/Nc;
    return 34./3 * Nc*Nc - 10./3 * Nf*Nc - 2*CF*Nf;
  }

  //V is the matrix that diagonalizes \gamma_s^{(0)T}  via  V^-1 \gamma^{(0)T} V
  static MatrixD computeV(const MatrixD &g0){
    std::vector<NumericVector<std::complex<double> > > evecs;
    std::vector<std::complex<double> > evals;
    nonSymmetricMatrixEigensolve(evecs, evals, g0.transpose(), true);
    
    //Note the structure of g0 fortunately results in real eigenvectors despite being a non-symmetric matrix

    int n = g0.size();
    MatrixD out(n); //columns are eigenvectors
    for(int i=0;i<n;i++)
      for(int j=0;j<n;j++){
	assert(fabs(evecs[j](i).imag()) < 1e-10);
	double v = evecs[j](i).real();
	if(fabs(v) < 1e-10) v = 0.;

	out(i,j) = v;
      }
    return out;
  }

  static MatrixD computeg0D(const MatrixD &g0, const MatrixD &V, const MatrixD &Vinv){
    MatrixD g0D = Vinv * g0.transpose() * V;
    for(int i=0;i<g0.size();i++)
      for(int j=0;j<g0.size();j++)
	if(i!=j && fabs(g0D(i,j)) > 1e-10) assert(0);
    return g0D;
  }

  static MatrixD computeG(const MatrixD &g1, const MatrixD &V, const MatrixD &Vinv){
    MatrixD G = Vinv * g1.transpose() * V;
    return G;
  }
  static MatrixD computeH(const MatrixD &g0D, const MatrixD &G, 
			 const MatrixD &V, const MatrixD &Vinv,
			 const double beta0, const double beta1){
    const int N = g0D.size();
    MatrixD H(N);
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++){
	H(i,j) =  - G(i,j)/(2*beta0+g0D(i,i)-g0D(j,j) );
	if(i==j) H(i,j) += g0D(i,i)*beta1/2/beta0/beta0;
      }
    }
    return H;
  }

  //J is defined in Eq III.98 of https://arxiv.org/pdf/hep-ph/9512380.pdf
  static MatrixD computeJ(const MatrixD &H, const MatrixD &V, const MatrixD &Vinv){
    MatrixD J = V * H * Vinv;
    return J;
  }
  
  //U^(0)(M1, M2) as defined in Eq III.94 of https://arxiv.org/pdf/hep-ph/9512380.pdf
  //asM* = \alpha_s(M*)    M2 is a higher energy scale
  static MatrixD computeU0(const double asM1, const double asM2, 
			  const double beta0, const MatrixD &g0D,
			  const MatrixD &V, const MatrixD &Vinv){
    const int N  = g0D.size();
    MatrixD U0D(N);
    U0D.zero();
    for(int i=0;i<N;i++) U0D(i,i) = pow(asM2/asM1, g0D(i,i)/2/beta0);
    
    MatrixD U0 = V * U0D * Vinv;
    return U0;
  }

  //Compute U via eq III.93 of https://arxiv.org/pdf/hep-ph/9512380.pdf, ignoring as^2 terms
  //This is the pure QCD scale evolution
  static MatrixO computeU(const MatrixD &U0, const  MatrixD &J,
			 const double asM1, const double asM2){
    MatrixO U = U0 + _AS/(4*M_PI)*( asM1*J*U0 - asM2*U0*J ); 
    return U;
  }

  //Puts all the pieces together for pure QCD scale evolution
  static MatrixO computeU(const double asM1, const double asM2, const int u, const int d, const int Nc, const bool shift_beta0_3f = true){
    int Nf = u+d;

    double beta0 = computeBeta0(Nf, Nc);
    if(Nf == 3 && shift_beta0_3f) beta0 += 1e-08;
    double beta1 = computeBeta1(Nf, Nc);
    
    MatrixD g0 = gammas0(Nc, Nf, u,d);
    MatrixD g1 = gammas1(Nc, Nf, u,d);

    MatrixD V = computeV(g0);
    MatrixD Vinv(V); svd_inverse(Vinv,V); 

    //LO
    MatrixD g0D = computeg0D(g0, V, Vinv);
    MatrixD U0 = computeU0(asM1, asM2, beta0, g0D, V, Vinv);   

    //NLO
    MatrixD G = computeG(g1, V, Vinv);    
    MatrixD H = computeH(g0D, G, V, Vinv, beta0, beta1);
    MatrixD J = computeJ(H, V, Vinv);

    //Sum
    MatrixO U = computeU(U0, J, asM1, asM2);
    
    return U;
  }
  
  //EM corrections to evolution, defined in VII.23 of https://arxiv.org/pdf/hep-ph/9512380.pdf
  static MatrixO computeR(const MatrixD &g0D,
			 const MatrixD &ge0, const MatrixD &gse1, 
			 const MatrixD &J, const MatrixD &H, 
			 const MatrixD &V, const MatrixD &Vinv,
			 const double asM1, const double asM2, const double beta0, const double beta1){
    int N = g0D.size();
    std::vector<double> va(N);
    for(int i=0;i<N;i++) va[i] = g0D(i,i)/2/beta0;
    
    MatrixD ge0T = ge0.transpose();

    MatrixD M0 = Vinv * ge0T * V;
    MatrixD M1= Vinv * (gse1.transpose()-beta1/beta0 * ge0T  +  ge0T*J-J*ge0T ) * V; //Eq VII.28

    MatrixD K0(N);
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++)
	K0(i,j) = M0(i,j)/(va[i] - va[j] - 1) * ( pow(asM2/asM1, va[j])/asM1   - pow(asM2/asM1, va[i])/asM2 ); //VII.24

    MatrixD K11(N);
    for(int i=0;i<N;i++)
      for(int j=0;j<N;j++){
	if(fabs(va[i] - va[j])<1e-10){
	  K11(i,j) = M1(i,j) * pow(asM2/asM1, va[i]) * log(asM1/asM2);
	}else{
	  K11(i,j) = M1(i,j)/(va[i]-va[j]) * ( pow(asM2/asM1, va[j]) - pow(asM2/asM1, va[i]) );
	}
      }
    MatrixD K12 = -asM2 * K0 * H; //K0 is O(as^-1) so K12 is O(as^0)
    MatrixD K13 = asM1 * H * K0;
    
    MatrixO R = -2*M_PI/beta0 * V * (K0*_ASM1 + 1./(4*M_PI)*(K11+K12+K13) ) * Vinv; //Eq VII.23
 
    return R;
  }

  static MatrixO computeUQCDQED(const MatrixO &U,
			       const MatrixO &R,
			       const double a){
    MatrixO Ufull = U + _AE * a/(4 * M_PI) * R; //VII.22
    return Ufull;
  }

  //if cc_only, evaluate using only the current-current operators (2x2)
  //In 3f theory we get nans unless we shift beta0 by a small amount (copied from Christoph, presumably some limit issue)
  static MatrixO computeUQCDQED(const double asM1, const double asM2, const int u, const int d, const int Nc, const double a, const bool cc_only = false, const bool shift_beta0_3f = true){
    int Nf = u+d;

    double beta0 = computeBeta0(Nf, Nc);
    if(Nf == 3 && shift_beta0_3f) beta0 += 1e-08;

    double beta1 = computeBeta1(Nf, Nc);
    
    MatrixD g0 = gammas0(Nc, Nf, u,d);
    MatrixD g1 = gammas1(Nc, Nf, u,d);
    MatrixD ge0 = gammae0(Nc, u, d);
    MatrixD gse1 = gammase1(Nc, u, d);

    if(cc_only){
      g0 = g0.submatrix(0,0,2);
      g1 = g1.submatrix(0,0,2);
      ge0 = ge0.submatrix(0,0,2);
      gse1 = gse1.submatrix(0,0,2);
    }

    MatrixD V = computeV(g0);
    MatrixD Vinv(V); svd_inverse(Vinv,V); 

    //LO
    MatrixD g0D = computeg0D(g0, V, Vinv);
    MatrixD U0 = computeU0(asM1, asM2, beta0, g0D, V, Vinv);   

    //NLO
    MatrixD G = computeG(g1, V, Vinv);    
    MatrixD H = computeH(g0D, G, V, Vinv, beta0, beta1);
    MatrixD J = computeJ(H, V, Vinv);

    MatrixO R = computeR(g0D, ge0, gse1, J, H, V, Vinv, asM1, asM2, beta0, beta1);

    //Sum
    MatrixO U = computeU(U0, J, asM1, asM2);
    MatrixO Ufull = computeUQCDQED(U, R, a);

    return Ufull;
  }

  //Matching matrix between 5 and 4 flavor theory at mb  (VII.30 - VII.36)
  static MatrixO MatMB(const double as4f_mb, const double a){
    MatrixD unit(10, [&](const int i, const int j){ return i==j ? 1: 0; });
    
    double P[10] = { 0,0,-1./3,1,-1./3,1,0,0,0,0 };
    double R[10] = {0,0,0,-2,0,-2,0,1,0,1};
    
    double Pbar[10] = {0,0,0,0,0,0,1,0,1,0};
    double Rbar[10] = {0,0,6,2,6,2,-3,-1,-3,-1};

    MatrixD PR(10, [&](const int i, const int j){ return P[i]*R[j]; });
    MatrixD PbarRbar(10, [&](const int i, const int j){ return Pbar[i]*Rbar[j]; });

    MatrixO out = unit + _AS * as4f_mb/(4*M_PI) * 5./18 * PR + _AE* a/(4*M_PI) * 10./81 * PbarRbar;
    return out;
  }
  //Matching matrix between 4 and 3 flavor theory at mc  (VII.30 - VII.36)
  static MatrixO MatMC(const double as4f_mc, const double a){
    MatrixD unit(10, [&](const int i, const int j){ return i==j ? 1: 0; });
    
    double P[10] = {0,0,-1./3,1,-1./3,1,0,0,0,0};
    double R[10] = {0,0,0,1,0,1,0,1,0,1};
    
    double Pbar[10] = {0,0,0,0,0,0,1,0,1,0};
    double Rbar[10] = {0,0,3,1,3,1,3,1,3,1};

    MatrixD PR(10, [&](const int i, const int j){ return P[i]*R[j]; });
    MatrixD PbarRbar(10, [&](const int i, const int j){ return Pbar[i]*Rbar[j]; });

    MatrixO out = unit + _AS * as4f_mc/(4*M_PI) * -5./9 * PR + _AE * a/(4*M_PI) * -40./81 * PbarRbar;
    return out;
  }
    
};



CPSFIT_END_NAMESPACE
#endif

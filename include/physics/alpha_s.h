#ifndef _ALPHA_S_COMPUTE_H
#define _ALPHA_S_COMPUTE_H

#include<cassert>
#include<gsl/gsl_errno.h>
#include<gsl/gsl_min.h>

#include<config.h>
#include<utils/macros.h>
#include<parser/parser.h>

SARLAC_START_NAMESPACE

#define A_S_T_MEMBERS (double, value)(double, scale)(int, Nf)
struct a_s_t{
  GENERATE_MEMBERS(A_S_T_MEMBERS);
};
GENERATE_PARSER(a_s_t,A_S_T_MEMBERS);


#define PERTURBATIVE_INPUTS_MEMBERS \
  (int, Nc)			    \
  (double, mcmc)		    \
  (double, mbmb)		    \
  (double, mWmW)		    \
  (double, mtmt)		    \
  (double, thetaW)		    \
  (double, a_e)			    \
  (a_s_t, a_s)

// double mcmc; //charm mass at charm scale
// double mbmb; //bottom mass at bottom scale
// double mWmW; //W mass at W scale
// double mtmt; //top mass at top scale
// double thetaW; //Weinberg angle
// double a_e; //EM coupling (assumed not to run)
// a_s_t a_s; //Input alpha_s
struct PerturbativeInputs{
  GENERATE_MEMBERS(PERTURBATIVE_INPUTS_MEMBERS);

  //Set to values from Qi's thesis
  void setQisValues();

  //Set to values from ep' PRL
  void setepPRLValues();

  //Set to values from Ziyuan Bai's thesis pg 67
  void setepZiyuanValues();

  //Set to values from epsilon' 2019 paper
  void setep2019values();

  inline double xt() const{ return pow(mtmt/mWmW,2); }

  PerturbativeInputs(){
    setepPRLValues();
  }

};
GENERATE_PARSER(PerturbativeInputs, PERTURBATIVE_INPUTS_MEMBERS);

//Storage for a number of variables appropriate for perturbation theory
class PerturbativeVariables{
  double Lambda[3]; //3,4,5 flavor Lambda_QCD
  PerturbativeInputs inputs;
  bool vrb;

public:
  void initialize(const PerturbativeInputs &_inputs);

  PerturbativeVariables(const PerturbativeInputs &inputs, const bool vrb = true): vrb(vrb){
    initialize(inputs);
  }
  PerturbativeVariables(const bool vrb = true): vrb(vrb){}
  
  const PerturbativeInputs &getInputs() const{ return inputs; }
  double getLambda(const int Nf) const{
    return Lambda[Nf-3]; 
  }
};


//Use the two-loop formula of Rev.Mod.Phys. 68 (1996) 1125-1144 to compute alpha_s
//All energies in *GeV*
struct ComputeAlphaS{
  inline static double beta0(const int Nf, const int Nc = 3){
    return (11.*Nc - 2*Nf)/3.;
  }
  inline static double CF(const int Nc = 3){
    return double(Nc*Nc - 1)/2/Nc;
  }
  inline static double beta1(const int Nf, const int Nc = 3){
    return 34./3 * Nc*Nc - 10./3 * Nf*Nc - 2*CF(Nc)*Nf;
  }
  inline static double alpha_s(const double mu, const double Lambda, const int Nf, const int Nc=3){
    double b0 = beta0(Nf,Nc);
    double b1 = beta1(Nf,Nc);
    double log_mu2_over_L2 = log(mu*mu/Lambda/Lambda);
    
    double pref = 4.*M_PI/b0/log_mu2_over_L2;
    double brack = 1 - b1/b0/b0 * log(log_mu2_over_L2)/log_mu2_over_L2;
    //printf("log(mu^2/L^2)=%.8f, beta0=%.8f beta1=%.8f 0,1 loop part %.8f, 2-loop multiplier %.8f\n",log_mu2_over_L2,b0,b1,pref,brack);
    return pref*brack;
  }

  //Compute Lambda_QCD by matching alpha_s in Nf-flavor theory to an input value
  //match_at is an energy scale, match_to is a value of alpha_s
  static double computeLambda(const int Nf, const double match_at, const double match_to, const int Nc=3, const bool vrb = true);

  inline static double alpha_s(const double mu, const PerturbativeVariables &pv, const int Nf, const int Nc=3){
    return alpha_s(mu, pv.getLambda(Nf), Nf, Nc);
  }
};







SARLAC_END_NAMESPACE

#endif

#include<tensors.h>
#include<physics.h>
#include<random.h>





using namespace CPSfit;

double Lambda5sol(const double mb, const double Lambda4){
  double as4Mb = ComputeAlphaS::alpha_s(mb, Lambda4, 4, 3);
  double Lambda5 = ComputeAlphaS::computeLambda(5, mb, as4Mb, 3, false);
  return Lambda5;
}

int main(const int argc, const char** argv){
  RNG.initialize(1234);
 
  typedef NumericVector<double> VectorD;
  typedef BasicSquareMatrix<double> MatrixD;
  typedef NumericVector<ValOrd> VectorO;
  typedef BasicSquareMatrix<ValOrd> MatrixO;

  //Test compare alpha_s values with those in Christoph's notebook
  {
    std::cout << "Checking alpha_s against Christoph" << std::endl;
    int Nc = 3;
    double mb = 4.18;
    double Lambda4 = 0.3298655;
    double Lambda5 = Lambda5sol(mb, Lambda4);
    double as = ComputeAlphaS::alpha_s(91.1876, Lambda5, 5, Nc);
    std::cout << as << " expect " << 0.118405 << std::endl;

    Lambda4 = 0.325;
    as = ComputeAlphaS::alpha_s(2, Lambda4, 4, Nc);
    std::cout << as << " expect " << 0.306034 << std::endl;

    double MW = 80.385;
    Lambda4 = 0.325;
    Lambda5 = Lambda5sol(mb, Lambda4);
    as = ComputeAlphaS::alpha_s(MW, Lambda5, 5, Nc);
    std::cout << as << " expect " << 0.120402 << std::endl;
  }
    


  //Test compare to Christoph values of Wilson coeffs at MW -  exact match
  {
    std::cout << "\nChecking Wilson coeffs at MW against Christoph" << std::endl;
    double Lambda4 = 0.325;  
    double Lambda5 = Lambda5sol(4.18, Lambda4);

    double MW = 80.385;
    double a = 1./128; //alpha_EM

    double asMW = ComputeAlphaS::alpha_s(MW, Lambda5, 5, 3);

    VectorO mCatMW = DeltaS1WilsonCoeffs::mCatMW(asMW, pow(170./91, 2), a, 0.2);

    std::cout << mCatMW << std::endl;

    std::vector<std::string> christoph_vals = {
      "0.05269704770479427as",
      "1-0.0012088591597691704ae-0.017565682568264757as",
      "0.003233116759384578ae+0.0005847961951907358as",
      "-0.0017543885855722076as",
      "0.0005847961951907358as",
      "-0.0017543885855722076as",
      "0.0007363799452239127ae",
      "0",
      "-0.045181389730665246ae",
      "0"};
      
    std::cout << "Expect";
    for(int i=0;i<10;i++) std::cout << " " << christoph_vals[i];
    std::cout << std::endl;
  }
  //Reproduce Buchalla Table XI (maybe some slight procedural difference?)
  //z1, z2 at mc
  {
    std::cout << "\nChecking against Buchalla Table XI" << std::endl;
    double mc = 1.3;
    double mb = 4.4;
    double MW = 80.2;
    double mt = 170;
    double a = 1/129.;
    double thetaW = asin(0.23);
    int Nc = 3;

    double xt = pow(mt/MW,2);
    double Lambda4 = 0.325; //middle column

    double as4mc = ComputeAlphaS::alpha_s(mc, Lambda4, 4, Nc);
    double as4mb = ComputeAlphaS::alpha_s(mb, Lambda4, 4, Nc);

    double Lambda5 = ComputeAlphaS::computeLambda(5, mb, as4mb, Nc, false);
    double as5MW = ComputeAlphaS::alpha_s(MW, Lambda5, 5, Nc);

    VectorO z1z2atmc =  DeltaS1WilsonCoeffs::z1z2atmc(as4mc, as4mb, as5MW, xt, a, thetaW, Nc);
    VectorD z1z2atmc_v(2, [&](const int i){ return truncateO1(z1z2atmc(i)).value(); });

    std::cout << z1z2atmc_v << std::endl;
    std::cout << "Expect -0.412, 1.208" << std::endl;
  }

  {
    std::cout << "\nChecking \\vec z at mc against Christoph's values (note his numbers are rounded to 10^-3)" << std::endl;
    double mc = 1.275;
    double mb = 4.18;
    double MW = 80.385;
    double mt = 173.5; //170
    double a = 1/128.;
    double thetaW = asin(sqrt(0.23116));
    int Nc = 3;

    double xt = pow(mt/MW,2);
    double Lambda4 = 0.325; //middle column

    double as4mc = ComputeAlphaS::alpha_s(mc, Lambda4, 4, Nc);
    double as4mb = ComputeAlphaS::alpha_s(mb, Lambda4, 4, Nc);

    double Lambda5 = ComputeAlphaS::computeLambda(5, mb, as4mb, Nc, false);
    double as5MW = ComputeAlphaS::alpha_s(MW, Lambda5, 5, Nc);

    VectorO z1z2atmc =  DeltaS1WilsonCoeffs::z1z2atmc(as4mc, as4mb, as5MW, xt, a, thetaW, Nc);
    VectorO zatmc = DeltaS1WilsonCoeffs::zatmc(z1z2atmc, as4mc, a);

    VectorD zatmc_v(10, [&](const int i){ return zatmc(i).value(); });

    std::cout << zatmc_v << std::endl;
    std::cout << "Expect -0.421, 1.220, 0.005, -0.014, -0.005, -0.014, 0, 0, 0, 0" << std::endl;
  }
  { //Exact agreement
    std::cout << "\nChecking EM-penguin at 2 GeV against Christoph's values" << std::endl;
    //Although \vec y = \vec v - \vec z,      \vec z is zero for its last 4 entries
    
    double mu = 2;
    double mb = 4.18;
    double MW = 80.385;
    double mt = 173.5; //170
    double a = 1/128.;
    double thetaW = asin(sqrt(0.23116));
    int Nc = 3;

    double xt = pow(mt/MW,2);
    double Lambda4 = 0.325; //middle column
    
    double as4mu = ComputeAlphaS::alpha_s(mu, Lambda4, 4, Nc);
    double as4mb = ComputeAlphaS::alpha_s(mb, Lambda4, 4, Nc);

    double Lambda5 = ComputeAlphaS::computeLambda(5, mb, as4mb, Nc, false);
    double as5MW = ComputeAlphaS::alpha_s(MW, Lambda5, 5, Nc);

    VectorO mCat4fMu = DeltaS1WilsonCoeffs::mCat4fMu(as4mu, as4mb, as5MW, xt, a, thetaW, Nc, true);
    
    for(int i=6; i<10; i++) std::cout <<  truncateO1(mCat4fMu(i)).value()/a << " ";
    std::cout << "\nExpect -0.0146855, 0.0967222, -1.40952, 0.421068" << std::endl;
  }
  {
    std::cout << "\nChecking \\vec y at mc in 3f theory against Christoph's values" << std::endl;
    
    double mc = 1.275;
    double mb = 4.18;
    double MW = 80.385;
    double mt = 173.5; //170
    double a = 1/128.;
    double thetaW = asin(sqrt(0.23116));
    int Nc = 3;

    double xt = pow(mt/MW,2);
    double Lambda4 = 0.325; //middle column
    
    double as4mc = ComputeAlphaS::alpha_s(mc, Lambda4, 4, Nc);
    double as4mb = ComputeAlphaS::alpha_s(mb, Lambda4, 4, Nc);

    double Lambda5 = ComputeAlphaS::computeLambda(5, mb, as4mb, Nc, false);
    double as5MW = ComputeAlphaS::alpha_s(MW, Lambda5, 5, Nc);

    VectorO mCat4fMc = DeltaS1WilsonCoeffs::mCat4fMu(as4mc, as4mb, as5MW, xt, a, thetaW, Nc, true);
    VectorO mCat3fMc = DeltaS1WilsonCoeffs::mCat3fMc(mCat4fMc, as4mc, a, true);

    VectorO z1z2mc = DeltaS1WilsonCoeffs::z1z2atmc(as4mc, as4mb, as5MW, xt, a, thetaW, Nc);
    VectorO zmc = DeltaS1WilsonCoeffs::zatmc(z1z2mc, as4mc, a);
    VectorO ymc = DeltaS1WilsonCoeffs::yatmu(zmc, mCat3fMc);

    for(int i=0;i<6;i++)
      std::cout << ymc(i).value() << " ";
    for(int i=6;i<10;i++)
      std::cout << ymc(i).value()/a << " ";
    std::cout << std::endl;

    std::cout << "Expect 0 0 0.029 -0.058 0.005 -0.091 -0.026 0.143 -1.506 0.564" << std::endl;
  }
  {
    std::cout << "\nChecking \\vec y and \\vec z at 1 GeV in 3f theory against Christoph's values" << std::endl;
    double mu = 1.;
    double mc = 1.275;
    double mb = 4.18;
    double MW = 80.385;
    double mt = 173.5; //170
    double a = 1/128.;
    double thetaW = asin(sqrt(0.23116));
    int Nc = 3;

    double xt = pow(mt/MW,2);
    double Lambda4 = 0.325; //middle column
    
    double as4mc = ComputeAlphaS::alpha_s(mc, Lambda4, 4, Nc);
    double as4mb = ComputeAlphaS::alpha_s(mb, Lambda4, 4, Nc);

    double Lambda5 = ComputeAlphaS::computeLambda(5, mb, as4mb, Nc, false);
    double as5MW = ComputeAlphaS::alpha_s(MW, Lambda5, 5, Nc);

    double Lambda3 = ComputeAlphaS::computeLambda(3, mc, as4mc, Nc, false);
    double as3mu = ComputeAlphaS::alpha_s(mu, Lambda3, 3, Nc);

    MatrixO U3MCtoMU = DeltaS1anomalousDimension::computeUQCDQED(as3mu, as4mc, 1,2, Nc, a);

    VectorO mCat4fMc = DeltaS1WilsonCoeffs::mCat4fMu(as4mc, as4mb, as5MW, xt, a, thetaW, Nc, true);
    VectorO mCat3fMc = DeltaS1WilsonCoeffs::mCat3fMc(mCat4fMc, as4mc, a, true);
    VectorO mCat3fMu = DeltaS1WilsonCoeffs::mCat3fMu(mCat3fMc, U3MCtoMU, true);


    VectorO z1z2mc = DeltaS1WilsonCoeffs::z1z2atmc(as4mc, as4mb, as5MW, xt, a, thetaW, Nc);
    VectorO zmc = DeltaS1WilsonCoeffs::zatmc(z1z2mc, as4mc, a);
    VectorO zmu = DeltaS1WilsonCoeffs::zat3fmu(zmc, U3MCtoMU);
    
    VectorO ymu = DeltaS1WilsonCoeffs::yatmu(zmu, mCat3fMu);

    std::cout << "z" << std::endl;
    for(int i=0;i<6;i++)
      std::cout << zmu(i).value() << " ";
    for(int i=6;i<10;i++)
      std::cout << zmu(i).value()/a << " ";
    std::cout << std::endl;

    std::cout << "Expect -0.509 1.278 0.012 -0.034 0.007 -0.033 0.012 0.013 0.018 -0.008" << std::endl;

    std::cout << "\ny" << std::endl;
    for(int i=0;i<6;i++)
      std::cout << ymu(i).value() << " ";
    for(int i=6;i<10;i++)
      std::cout << ymu(i).value()/a << " ";
    std::cout << std::endl;

    std::cout << "Expect 0 0 0.328 -0.0584 -0.0004 -0.1117 -0.0273 0.1795 -1.5981 0.6999" << std::endl;
  }
 {
    std::cout << "\nChecking \\vec y and \\vec z at 2.15 GeV in 3f theory against Christoph's values from Qi's thesis" << std::endl;
    PerturbativeInputs inputs;
    inputs.setQisValues();

    PerturbativeVariables var(inputs);
    
    double mu = 2.15;
    
    double as4mc = ComputeAlphaS::alpha_s(inputs.mcmc, var.getLambda(4), 4, inputs.Nc);
    double as4mb = ComputeAlphaS::alpha_s(inputs.mbmb, var.getLambda(4), 4, inputs.Nc);
    double as5MW = ComputeAlphaS::alpha_s(inputs.mWmW, var.getLambda(5), 5, inputs.Nc);
    double as3mu = ComputeAlphaS::alpha_s(mu, var.getLambda(3), 3, inputs.Nc);

    MatrixO U3MCtoMU = DeltaS1anomalousDimension::computeUQCDQED(as3mu, as4mc, 1,2, inputs.Nc, inputs.a_e);

    VectorO mCat4fMc = DeltaS1WilsonCoeffs::mCat4fMu(as4mc, as4mb, as5MW, inputs.xt(), inputs.a_e, inputs.thetaW, inputs.Nc, true);
    VectorO mCat3fMc = DeltaS1WilsonCoeffs::mCat3fMc(mCat4fMc, as4mc, inputs.a_e, true);
    VectorO mCat3fMu = DeltaS1WilsonCoeffs::mCat3fMu(mCat3fMc, U3MCtoMU, true);


    VectorO z1z2mc = DeltaS1WilsonCoeffs::z1z2atmc(as4mc, as4mb, as5MW, inputs.xt(), inputs.a_e, inputs.thetaW, inputs.Nc);
    VectorO zmc = DeltaS1WilsonCoeffs::zatmc(z1z2mc, as4mc, inputs.a_e);
    VectorO zmu = DeltaS1WilsonCoeffs::zat3fmu(zmc, U3MCtoMU);
    
    VectorO ymu = DeltaS1WilsonCoeffs::yatmu(zmu, mCat3fMu);

    std::cout << "z" << std::endl;
    for(int i=0;i<10;i++)
      std::cout << zmu(i).value() << " ";
    std::cout << std::endl;

    std::cout << "Expect -0.29829 1.14439 -0.00243827 0.00995157 -0.00110544 0.00657457 0.0000701587 -0.0000900499 0.0000150176 0.0000656482" << std::endl;

    std::cout << "\ny" << std::endl;
    for(int i=0;i<10;i++)
      std::cout << ymu(i).value() << " ";
    std::cout << std::endl;

    std::cout << "Expect 0 0 0.024141 -0.058121 0.0102484 -0.0699711 -0.000211182 0.000776697 -0.0106787 0.0029815" << std::endl;
 }

 {
   std::cout << "\nChecking \\vec y and \\vec z at 1.531 GeV in 3f theory against 2015 ep' PRL" << std::endl;
   double mu = 1.531;

   PerturbativeInputs inputs;
   inputs.setepPRLValues();

   PerturbativeVariables var(inputs);
   
   VectorD y,z;
   DeltaS1WilsonCoeffs::compute3fWilsonCoeffs(y,z,mu,var);

   std::cout << "z" << std::endl;
   for(int i=0;i<10;i++)
     std::cout << z(i) << " ";
   std::cout << std::endl;

   std::cout << "Expect -0.3734 1.189 0.001294 -0.003691 0.002505 -0.004718 6.099e-05 -4.414e-05 3.561e-05 3.010e-05" << std::endl;

   std::cout << "\ny" << std::endl;
   for(int i=0;i<10;i++)
     std::cout << y(i) << " ";
   std::cout << std::endl;

   std::cout << "Expect 0 0 0.02720 -0.05828 0.007133 -0.08191 -0.0003401 0.0008684 -0.01046 0.003496" << std::endl;
 }
 {
   std::cout << "\nChecking \\vec y and \\vec z at 2.29 GeV in 3f theory against those given in Greg McGlynn's thesis pg 205/206 which are from Ziyuan " << std::endl;
   double mu = 2.29;

   PerturbativeInputs inputs;
   inputs.setepZiyuanValues();

   PerturbativeVariables var(inputs);
   
   VectorD y,z;
   DeltaS1WilsonCoeffs::compute3fWilsonCoeffs(y,z,mu,var);

   std::cout << "z" << std::endl;
   for(int i=0;i<10;i++)
     std::cout << z(i) << " ";
   std::cout << std::endl;

   std::cout << "Expect -0.2879390507 1.138420148 -0.002847074467 0.01169943305 -0.001628070946 0.007936074279 7.371992154e-05 -9.489150846e-05 1.513586959e-05 6.97390239e-05" << std::endl;

   std::cout << "\ny" << std::endl;
   for(int i=0;i<10;i++)
     std::cout << y(i) << " ";
   std::cout << std::endl;

   std::cout << "Expect 0 0 0.02344132124 -0.05789935776 0.01059643117 -0.06829450251 -0.0003202699843 0.0006886695743 -0.00994958219 0.002664737778" << std::endl;
 }







  std::cout << "Done\n";
  return 0;
}



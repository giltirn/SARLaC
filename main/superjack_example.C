#include <fstream>
#include <minimizer.h>
#include <cost_function.h>
#include <fitfunc.h>
#include <fitfunc_mapping.h>
#include <random.h>
#include <plot.h>
#include <distribution.h>
#include <data_series.h>
#include <parser.h>
#include <common_defs.h>
#include <sstream>
#include <superjackknife.h>

void gaussianJackknifePreciseMean(jackknifeDistributionD &d, double mean, double std_err){
  double std_dev = std_err/sqrt( d.size()-1. );
  for(int i=0;i<d.size();i++)
    d.sample(i) = gaussianRandom<double>(mean, std_dev);
  double delta = mean - d.mean();
  for(int i=0;i<d.size();i++)
    d.sample(i) += delta;    
}

inline std::string print(const superJackknifeDistribution<double> &sja){
  std::ostringstream os;
  os << "(" << sja.mean() << "," << sja.standardError() << ") = ";
  for(int i=0;i<3;i++)
    os << "[" << sja.getLayout().nSamplesEns(i) << "](" << sja.getEnsembleJackknife(i).mean() << "," << sja.getEnsembleJackknife(i).standardError() << ") ";
  return os.str();
}

int main(void){  
  RNG.initialize(1234);

  int nsamples[3] = {20,40,60};

  jackknifeDistributionD a(nsamples[0]);  gaussianJackknifePreciseMean(a,1.0,0.5);
  jackknifeDistributionD b(nsamples[1]);  gaussianJackknifePreciseMean(b,2.0,1.0);
  jackknifeDistributionD c(nsamples[2]);  gaussianJackknifePreciseMean(c,4.0,2.0);

  std::cout << "a (" << a.mean() << "," << a.standardError() << ")\n";
  std::cout << "b (" << b.mean() << "," << b.standardError() << ")\n";
  std::cout << "c (" << c.mean() << "," << c.standardError() << ")\n";
  
  superJackknifeLayout layout;
  layout.addEnsemble("A",nsamples[0]);
  layout.addEnsemble("B",nsamples[1]);
  layout.addEnsemble("C",nsamples[2]);
  
  superJackknifeDistribution<double> sja(layout, 0, a);
  superJackknifeDistribution<double> sjb(layout, 1, b);
  superJackknifeDistribution<double> sjc(layout, 2, c);

  std::cout << "sja " << print(sja) << "\nsjb " << print(sjb) << "\nsjc " << print(sjc) << '\n';

  {
    superJackknifeDistribution<double> sjd = sja + sjb;
    double err_expect = sqrt( pow(a.standardError(),2) + pow(b.standardError(),2) );
    std::cout << "sja + sjb = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    superJackknifeDistribution<double> sjd = sja - sjb;
    double err_expect = sqrt( pow(a.standardError(),2) + pow(b.standardError(),2) );
    std::cout << "sja - sjb = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    superJackknifeDistribution<double> sjd = sja * sjb;
    double err_expect = sqrt( pow(a.standardError() * b.mean()  ,2) + pow(a.mean() * b.standardError(),2) );
    std::cout << "sja * sjb = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    superJackknifeDistribution<double> sjd = sja / sjb;
    double err_expect = sqrt( pow(a.standardError() / b.mean()  ,2) + pow( a.mean() * b.standardError() / pow(b.mean(),2) ,2) );
    std::cout << "sja / sjb = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    jackknifeDistributionD t(nsamples[0]); gaussianJackknifePreciseMean(t,5.0,2.5);
    superJackknifeDistribution<double> sjt(layout, 0, t);
    
    superJackknifeDistribution<double> sjd = sja * sjt;
    double err_expect = sqrt( pow(a.standardError() * t.mean()  ,2) + pow(a.mean() * t.standardError(),2) );
    std::cout << "sjt " << print(sjt) << std::endl;
    std::cout << "sja * sjt = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    superJackknifeDistribution<double> sjd = 3.*sja / sjb;
    double err_expect = sqrt( pow(3.*a.standardError() / b.mean()  ,2) + pow( 3.*a.mean() * b.standardError() / pow(b.mean(),2) ,2) );
    std::cout << "3*sja / sjb = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    superJackknifeDistribution<double> sjd = sqrt(sjc);
    double err_expect = 0.5/sqrt(sjc.mean())*sjc.standardError();
    std::cout << "sqrt(sjc) = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  {
    superJackknifeDistribution<double> sjd = sqrt(sjc) - sjb;
    double err_expect = sqrt(
  			       pow(  0.5/sqrt(sjc.mean())*sjc.standardError(),  2)
  			     + pow(  sjb.standardError(),                       2)
  			     ); 
    std::cout << "sqrt(sjc) - sjb = " << print(sjd) << " err expect " << err_expect << std::endl;
  }
  
  std::cout << "Normal exit\n"; std::cout.flush();
  return 0;
}


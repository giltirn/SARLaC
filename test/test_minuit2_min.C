#include<config.h>

#ifndef HAVE_MINUIT2

#include<iostream>

int main(const int argc, const char **argv){
  std::cout << "Library not compiled with Minuit2" << std::endl;
  return 0;
}

#else

#include<random.h>
#include<minimizer.h>
#include<fit.h>

using namespace SARLaC;

template<typename CostFunc>
void test(const CostFunc &cost,
	  const StandardFitParams &guess){
    //Regular
    StandardFitParams params(guess);
    MarquardtLevenbergParameters<double> p; p.verbose = true; p.delta_cost_min = 1e-10;    
    MarquardtLevenbergMinimizer<CostFunc> min(cost, p);    
    double chisq = min.fit(params);
    
    std::cout << "ML : Chisq = " << chisq << "\nParams = " << params << std::endl;

    //Minuit2

    StandardFitParams params_min2(guess);
    Minuit2minimizerParams p_min2;
    p_min2.tolerance = 1e-10;
    p_min2.verbose = true;

    Minuit2minimizer<CostFunc> min_min2(cost, p_min2);
    
    double chisq_min2 = min_min2.fit(params_min2);
    
    std::cout << "Minuit2 : Chisq = " << chisq_min2 << "\nParams = " << params_min2 << std::endl;

    double reltol = 1e-05;
    
    double dchisq = fabs(chisq - chisq_min2)/chisq;
    std::cout << "Relative chi^2 difference " << dchisq << " tolerance "<< reltol << std::endl;
    //if(dchisq > reltol) assert(0);
    
    for(int i=0;i<params.size();i++){
      double dp = fabs(params(i) - params_min2(i))/params(i);
      std::cout << "Relative p(" << i << ") difference " << dp << " tolerance "<< reltol << std::endl;
      //if(dp > reltol) assert(0);
    }
}



int main(const int argc, const char** argv){
  RNG.initialize(1234);
  
  double A = 1.0;
  double E = 0.3;
  int ndata = 10;

  std::vector<double> sigma(ndata, 0.07);

  typedef correlationFunction<double, double> CorrFunc;
  CorrFunc data(ndata);
  for(int t=0;t<ndata;t++){
    data.coord(t) = t;
    data.value(t) = gaussianRandom<double>(0,0.07) * A * exp(-E*t);
  }
				       
  NumericSquareMatrix<double> corr(ndata);
  for(int i=0;i<ndata;i++)
    for(int j=0;j<ndata;j++)
      corr(i,j) = uniformRandom<double>(-0.2,0.2);
  
  corr = (corr + corr.transpose())/2.;
  for(int i=0;i<ndata;i++) corr(i,i) = 1.;

  std::cout << "Correlation matrix:\n" << corr;
  
  NumericSquareMatrix<double> inv_corr(ndata);
  svd_inverse(inv_corr,corr);

  typedef FitExp FitFunc;

  StandardFitParams guess(0.7, 0.5);
  FitFunc fitfunc;

  //Test uncorrelated fit
  {
    std::cout << "Testing uncorrelated fit" << std::endl;
    typedef UncorrelatedChisqCostFunction<FitFunc, CorrFunc> CostFunc;
    CostFunc cost(fitfunc, data, sigma);
    
    test(cost, guess);
  }


  //Test correlated fit
  {
    std::cout << "Testing correlated fit" << std::endl;
    typedef CorrelatedChisqCostFunction<FitFunc, CorrFunc> CostFunc;
    CostFunc cost(fitfunc, data, sigma, inv_corr);

    test(cost, guess);
  }

  std::cout << "Done\n";
  return 0;
}

#endif

#include<random.h>
#include<fit.h>

using namespace CPSfit;


template<typename CostFunc, typename CostFuncGSL>
void test(const CostFunc &cost, const CostFuncGSL &cost_gsl,
	  const StandardFitParams &guess){
    //Regular
    StandardFitParams params(guess);
    MarquardtLevenbergParameters<double> p; p.verbose = true; p.delta_cost_min = 1e-10;    
    MarquardtLevenbergMinimizer<CostFunc> min(cost, p);    
    double chisq = min.fit(params);
    
    std::cout << "ML : Chisq = " << chisq << "\nParams = " << params << std::endl;

    //GSL
    StandardFitParams params_gsl(guess);
    GSLtrsMinimizerParams p_gsl;
    p_gsl.stopping_conditions = { {GSLtrsStoppingCondition::StopCostDelta, 1e-10} };
    p_gsl.dampening_matrix = MLdampeningMatrix::HessianDiag;
    p_gsl.verbose = true;
    p_gsl.algorithm = GSLtrsAlgorithm::MarquardtLevenberg; 

    GSLtrsMinimizer<CostFuncGSL> min_gsl(cost_gsl, p_gsl);
    
    double chisq_gsl = min_gsl.fit(params_gsl);
    
    std::cout << "GSL ML : Chisq = " << chisq_gsl << "\nParams = " << params_gsl << std::endl;

    double reltol = 1e-05;
    
    double dchisq = fabs(chisq - chisq_gsl)/chisq;
    std::cout << "Relative chi^2 difference " << dchisq << " tolerance "<< reltol << std::endl;
    if(dchisq > reltol) assert(0);
    
    for(int i=0;i<params.size();i++){
      double dp = fabs(params(i) - params_gsl(i))/params(i);
      std::cout << "Relative p(" << i << ") difference " << dp << " tolerance "<< reltol << std::endl;
      if(dp > reltol) assert(0);
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

  std::vector<double> corr_evals;
  std::vector<NumericVector<double> > corr_evecs;
  symmetricMatrixEigensolve(corr_evecs, corr_evals, corr);

  typedef FitExp FitFunc;


  StandardFitParams guess(0.7, 0.5);
  FitFunc fitfunc;

  //Test uncorrelated fit
  {
    std::cout << "Testing uncorrelated fit" << std::endl;
    typedef UncorrelatedChisqCostFunction<FitFunc, CorrFunc> CostFunc;
    CostFunc cost(fitfunc, data, sigma);
    
    typedef UncorrelatedChisqCostFunctionTerms<FitFunc, CorrFunc> CostFuncGSL;
    CostFuncGSL cost_gsl(fitfunc, data, sigma);
    
    test(cost, cost_gsl, guess);
  }


  //Test correlated fit
  {
    std::cout << "Testing correlated fit" << std::endl;
    typedef CorrelatedChisqCostFunction<FitFunc, CorrFunc> CostFunc;
    CostFunc cost(fitfunc, data, sigma, inv_corr);
    
    typedef CorrelatedChisqCostFunctionTerms<FitFunc, CorrFunc> CostFuncGSL;
    CostFuncGSL cost_gsl(fitfunc, data, sigma, corr_evals, corr_evecs);
        
    test(cost, cost_gsl, guess);
  }

  std::cout << "Done\n";
  return 0;
}


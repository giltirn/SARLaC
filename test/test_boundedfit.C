#include<distribution.h>
#include<fit.h>
#include<data_series.h>
#include<common.h>
#include<random.h>
#include<plot.h>

using namespace SARLaC;

int main(void){
  RNG.initialize(1234);

  int Lt = 32;
  int nsample = 200;

  double snr_base = 50; 
  double snr_decay_exp = 0.15;

  double A0 = 1e6;
  double E0 = 0.35;

  double A1 = 1e6;
  double E1 = 0.6;

  rawDataCorrelationFunctionD raw(Lt);
  for(int t=0;t<Lt;t++){
    double v = A0 * exp(-E0 * t) + A1 * exp(-E1 * t);

    double snr = snr_base * exp(-snr_decay_exp * t);
    double err = v/snr * sqrt(double(nsample-1)); //Error on mean is sqrt(N-1) smaller than raw data error

    raw.coord(t) = t;
    raw.value(t).resize(nsample);
    gaussianRandom(raw.value(t), v, err);

    std::cout << t  << " " << v << " " << err << " " << raw.value(t) << std::endl;
  }

  jackknifeCorrelationFunctionD jack(Lt, [&](const int t){ return jackknifeCorrelationFunctionD::ElementType(t, jackknifeDistributionD(raw.value(t))); });

  NumericVector<jackknifeDistributionD> sigma(Lt, [&](const int t){ 
      doubleJackknifeDistributionD dj(raw.value(t));
      return jackknifeDistributionD(  sqrt(doubleJackknifeDistributionD::covariance(dj,dj)) );
    });

  std::cout << "As jackknife\n";
  for(int t=0;t<Lt;t++) std::cout << t << " " << jack.value(t) << std::endl;

  //Plot the effective mass
  jackknifeCorrelationFunctionD effmass;
  {
    FitExp ff; StandardFitParams base(1.0,1.0); 
    effmass = effectiveMass2pt(jack, ff, base, 1, Lt);
  }
  MatPlotLibScriptGenerate plot;
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
  accessor acc(effmass);
  plot.plotData(acc);
  plot.write("test_boundedfit_effmass.py", "test_boundedfit_effmass.pdf");
  
  //Try regular two-state fit
  typedef FitTwoStateExp FitFunc;
  FitFunc fitfunc;

  TwoStateFitParams guess(0.1e6, 0.1, 0.1e6, 0.2);

  {
    std::cout << "Doing regular two-state fit\n";

    typedef typename composeFitPolicy<FitFunc, standardFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;
    
    typedef typename fitter<FitPolicies>::minimizerParamsType MLparams;
    MLparams mlp;
    mlp.verbose = true;
    mlp.dampening_matrix = MLdampeningMatrix::MaxHessianDiag;

    fitter<FitPolicies> fitter(mlp);
    fitter.importFitFunc(fitfunc);
    fitter.importCostFunctionParameters(sigma);
    
    jackknifeDistribution<TwoStateFitParams> params(nsample, guess);
    jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);
    
    fitter.fit(params, chisq, chisq_per_dof, jack);
    
    std::cout << "Result of normal two-state fit\n";
    std::cout << "Params:\n" << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
    
    std::cout << "Done" << std::endl;
  }

  //Try constraining E1 to be > 0.7. We use a window fit to constrain the maximum E1 
  {
    std::cout << "Constraining 0.7 < E1 < 2\n";
    
    typedef typename composeFitPolicy<FitFunc, boundedFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;

    typedef typename fitter<FitPolicies>::minimizerParamsType MLparams;
    MLparams mlp;
    mlp.verbose = true;
    mlp.lambda_max = 1e16;
    mlp.lambda_factor = 2.;
    mlp.max_iter = 1e6;
    mlp.dampening_matrix = MLdampeningMatrix::HessianDiag;

    fitter<FitPolicies> fitter(mlp);
    fitter.importFitFunc(fitfunc);
    fitter.importCostFunctionParameters(sigma);
    fitter.setBound(3, ParameterBound::Window, 0.7,2);

    auto guess2 = guess;
    guess2.E1 = 0.8;

    jackknifeDistribution<TwoStateFitParams> params(nsample, guess2);
    jackknifeDistributionD chisq(nsample), chisq_per_dof(nsample);
    
    fitter.fit(params, chisq, chisq_per_dof, jack);

    std::cout << "Result of bound two-state fit with 0.7 < E1 < 2\n";
    std::cout << "Params:\n" << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
    
    std::cout << "Done" << std::endl;
  }

  return 0;
}

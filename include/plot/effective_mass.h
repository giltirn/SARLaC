#ifndef _EFFECTIVE_MASS_H_
#define _EFFECTIVE_MASS_H_

//Numerically invert a fit function with a zero-degree-of-freedom fit to obtain the effective mass

#include<config.h>
#include<utils/macros.h>
#include<distribution/jackknife.h>
#include<data_series/data_series.h>
#include<minimizer/minimizer.h>
#include<fit/cost_function.h>
#include<plot/effective_mass/fitfunc.h>

CPSFIT_START_NAMESPACE

//For a generic fit form and arbitrary linear combination
template<typename jackknifeTimeSeriesType, typename FitEffMassFunc>
jackknifeTimeSeriesType fitEffectiveMass(const jackknifeTimeSeriesType &edata, const FitEffMassFunc &fiteffmass){
  typedef UncorrelatedChisqCostFunction<FitEffMassFunc, dataSeries<double,double> > CostFunction;
  typedef MarquardtLevenbergMinimizer<CostFunction> MinimizerType;
  typedef typename FitEffMassFunc::ParameterType ParameterType;
  MarquardtLevenbergParameters<typename CostFunction::CostType> mlparams;
  mlparams.verbose = false;
  mlparams.exit_on_convergence_fail = false;
  std::vector<double> sigma_j(1,1.);

  const int nsample = edata.value(0).size();
  jackknifeTimeSeriesType effmass(edata.size(), nsample);
  auto orig_printer = distributionPrint<jackknifeDistribution<double> >::printer();
  distributionPrint<jackknifeDistribution<double> >::printer(new publicationDistributionPrinter<jackknifeDistribution<double> >,false);  

  std::vector<bool> erase(edata.size(),false);
  bool erase_required = false;

  for(int i=0;i<edata.size();i++){
    effmass.coord(i) = edata.coord(i);    
    
    std::vector<int> fail(omp_get_max_threads(),0); //sometimes the fit func inversion can fail on noisy data
#pragma omp parallel for
    for(int j=0;j<nsample;j++){
      dataSeries<double, double> rat_t_sample_j(1);
      rat_t_sample_j.coord(0) = edata.coord(i);
      rat_t_sample_j.value(0) = edata.value(i).sample(j);
      CostFunction costfunc(fiteffmass, rat_t_sample_j, sigma_j);
      MinimizerType fitter(costfunc, mlparams);
      ParameterType p(0.5);
      fitter.fit(p);
      if(!fitter.hasConverged()) fail[omp_get_thread_num()] = 1;
      else effmass.value(i).sample(j) = *p;
    }
    int did_fail = 0; for(int ii=0;ii<fail.size();ii++) did_fail += fail[ii];
    if(did_fail){
      std::cout << "Warning: Failed to converge on one or more samples at coord " << effmass.coord(i) << std::endl;
      erase[i] = true;
      erase_required = true;
    }else{
      jackknifeDistribution<double> y(nsample);
      jackknifeDistribution<double> yfit(nsample);
      jackknifeDistribution<double> resid(nsample);
      for(int s=0;s<nsample;s++){
	double t = edata.coord(i);
	y.sample(s) = edata.value(i).sample(s);
	singleValueContainer<double> m(effmass.value(i).sample(s));
	yfit.sample(s) = fiteffmass.value(t,m);

	resid.sample(s) = (y.sample(s) - yfit.sample(s))/y.sample(s);
      }
      std::cout << "Effmass t="<<  effmass.coord(i) << " m=" << effmass.value(i) << " y=" << y << " yfit=" << yfit << " resid=" << resid << std::endl;
    }
  }
  distributionPrint<jackknifeDistribution<double> >::printer(orig_printer);
    
  if(erase_required){
    int nkeep = 0; for(int i=0;i<edata.size();i++) if(!erase[i]) nkeep++;
    jackknifeTimeSeriesType tokeep(nkeep);
    int elem = 0;
    for(int i=0;i<edata.size();i++)
      if(!erase[i]){
	tokeep.coord(elem) = std::move(effmass.coord(i));
	tokeep.value(elem) = std::move(effmass.value(i));
	elem++;
      }
    return tokeep;
  }else return effmass;
}


//Two point effective mass for fit functions with form  A*f(m,t)
//Base should be a correctly setup parameter structure (needed so we can avoid requiring default constructors). The amplitude parameter should be set to a non-zero value. parameter_mass_index is the index of the mass parameter.
template<typename jackknifeTimeSeriesType, typename FitFunc>
jackknifeTimeSeriesType effectiveMass2pt(const jackknifeTimeSeriesType &data, const FitFunc &fitfunc, const typename FitFunc::ParameterType &base, const int parameter_mass_idx, const int Lt){
  if(data.size() == 0) return jackknifeTimeSeriesType(0);
  if(data.size() != Lt) error_exit(std::cout << "effectiveMass called with data of size " << data.size() << ". Expect 0 or Lt=" << Lt << std::endl);
  
  const int nsample = data.value(0).size();
  jackknifeTimeSeriesType ratios(Lt-1);
  for(int i=0;i<Lt-1;i++){
    double t = data.coord(i);
    assert(t == double(i));
    assert(data.coord(i+1) == t+1);

    ratios.coord(i) = t;
    ratios.value(i) = data.value(i)/data.value(i+1);
  }
  typedef Fit2ptEffectiveMass<FitFunc> FitEffMass;
  FitEffMass fiteffmass(fitfunc, base, parameter_mass_idx);
  return fitEffectiveMass<jackknifeTimeSeriesType,FitEffMass>(ratios,fiteffmass);
}

CPSFIT_END_NAMESPACE
#endif


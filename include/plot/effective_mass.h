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

SARLAC_START_NAMESPACE

//For a generic fit form and arbitrary linear combination
template<typename DistributionTimeSeriesType, typename FitEffMassFunc>
  DistributionTimeSeriesType fitEffectiveMass(const DistributionTimeSeriesType &edata, const FitEffMassFunc &fiteffmass, double guess = 0.5, bool verbose = false){
  typedef UncorrelatedChisqCostFunction<FitEffMassFunc, dataSeries<double,double> > CostFunction;
  typedef MarquardtLevenbergMinimizer<CostFunction> MinimizerType;
  typedef typename FitEffMassFunc::ParameterType ParameterType;
  typedef typename DistributionTimeSeriesType::DataType DistributionType;
  typedef iterate<DistributionType> Iter;
  
  MarquardtLevenbergParameters<typename CostFunction::CostType> mlparams;
  mlparams.verbose = false;
  mlparams.exit_on_convergence_fail = false;
  std::vector<double> sigma_j(1,1.);

  DistributionTimeSeriesType effmass(edata.size());
  auto orig_printer = distributionPrint<DistributionType>::printer();
  distributionPrint<DistributionType>::printer(new publicationDistributionPrinter<DistributionType>(2,Error),false);  

  std::vector<bool> erase(edata.size(),false);
  bool erase_required = false;

  DistributionType zero(edata.value(0)); zeroit(zero);

  int niter = Iter::size(zero);

  for(int i=0;i<edata.size();i++){
    effmass.coord(i) = edata.coord(i);    
    effmass.value(i) = edata.value(i);

    std::vector<int> fail(omp_get_max_threads(),0); //sometimes the fit func inversion can fail on noisy data
#pragma omp parallel for
    for(int j=0;j<niter;j++){
      dataSeries<double, double> rat_t_sample_j(1);
      rat_t_sample_j.coord(0) = edata.coord(i);
      rat_t_sample_j.value(0) = Iter::at(j, edata.value(i));
      CostFunction costfunc(fiteffmass, rat_t_sample_j, sigma_j);
      MinimizerType fitter(costfunc, mlparams);
      ParameterType p(guess);
      fitter.fit(p);
      if(!fitter.hasConverged()) fail[omp_get_thread_num()] = 1;
      else Iter::at(j, effmass.value(i)) = *p;
    }
    int did_fail = 0; for(int ii=0;ii<fail.size();ii++) did_fail += fail[ii];
    if(did_fail){
      if(verbose) std::cout << "Warning: Failed to converge on one or more samples at coord " << effmass.coord(i) << std::endl;
      erase[i] = true;
      erase_required = true;
    }else{
      DistributionType y(zero);
      DistributionType yfit(zero);
      DistributionType resid(zero);
      for(int s=0;s<niter;s++){
	double t = edata.coord(i);
	Iter::at(s, y) = Iter::at(s, edata.value(i));
	singleValueContainer<double> m(Iter::at(s, effmass.value(i)));
	Iter::at(s, yfit) = fiteffmass.value(t,m);

	Iter::at(s, resid) = (Iter::at(s, y) - Iter::at(s, yfit))/Iter::at(s, y);
      }
      if(verbose) std::cout << "Effmass t="<<  effmass.coord(i) << " m=" << effmass.value(i) << " y=" << y << " yfit=" << yfit << " resid=" << resid << std::endl;
    }
  }
  distributionPrint<DistributionType>::printer(orig_printer);
    
  if(erase_required){
    int nkeep = 0; for(int i=0;i<edata.size();i++) if(!erase[i]) nkeep++;
    DistributionTimeSeriesType tokeep(nkeep);
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
template<typename DistributionTimeSeriesType, typename FitFunc>
DistributionTimeSeriesType effectiveMass2pt(const DistributionTimeSeriesType &data, const FitFunc &fitfunc, const typename FitFunc::ParameterType &base, const int parameter_mass_idx, const int Lt, double guess = 0.5){
  if(data.size() == 0) return DistributionTimeSeriesType(0);
  if(data.size() != Lt) error_exit(std::cout << "effectiveMass called with data of size " << data.size() << ". Expect 0 or Lt=" << Lt << std::endl);
  
  DistributionTimeSeriesType ratios(Lt-1);
  for(int i=0;i<Lt-1;i++){
    double t = data.coord(i);
    assert(t == double(i));
    assert(data.coord(i+1) == t+1);

    ratios.coord(i) = t;
    ratios.value(i) = data.value(i)/data.value(i+1);
  }
  typedef Fit2ptEffectiveMass<FitFunc> FitEffMass;
  FitEffMass fiteffmass(fitfunc, base, parameter_mass_idx);
  return fitEffectiveMass<DistributionTimeSeriesType,FitEffMass>(ratios,fiteffmass,guess);
}

SARLAC_END_NAMESPACE
#endif


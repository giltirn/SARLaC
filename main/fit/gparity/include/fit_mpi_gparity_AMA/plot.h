#ifndef _FIT_MPI_GPARITY_AMA_PLOT_H
#define _FIT_MPI_GPARITY_AMA_PLOT_H

jackknifeTimeSeriesType effectiveMass(const jackknifeTimeSeriesType &data, const std::string &type, const int Lt, const TypeInfo &pmap){
  if(data.size() == 0) return jackknifeTimeSeriesType(0);
  if(data.size() != Lt) error_exit(std::cout << "effectiveMass called with data of size " << data.size() << ". Expect 0 or Lt=" << Lt << std::endl);
  
  const int nsample = data.value(0).size();
  dataSeries<double, jackknifeDistributionType> ratios(Lt-1);
  for(int i=0;i<Lt-1;i++){
    double t = data.coord(i);
    assert(t == double(i));
    assert(data.coord(i+1) == t+1);

    ratios.coord(i) = t;
    ratios.value(i) = data.value(i)/data.value(i+1);
  }
  typedef UncorrelatedChisqCostFunction<FitMpiEffectiveMass, dataSeries<double,double> > CostFunction;
  typedef MarquardtLevenbergMinimizer<CostFunction> MinimizerType;
  FitMpiEffectiveMass fiteffmass(2*Lt, type, pmap);
  
  MarquardtLevenbergParameters<typename CostFunction::CostType> mlparams;
  mlparams.verbose = false;
  mlparams.exit_on_convergence_fail = false;
  std::vector<double> sigma_j(1,1.);

  jackknifeTimeSeriesType effmass(ratios.size(), nsample);
  auto orig_printer = distributionPrint<jackknifeDistributionType>::printer();
  distributionPrint<jackknifeDistributionType>::printer(new publicationDistributionPrinter<jackknifeDistributionType>,false);  
  std::cout << "Effective mass for type " << type << "\n";

  std::vector<bool> erase(ratios.size(),false);
  bool erase_required = false;

  for(int i=0;i<ratios.size();i++){
    effmass.coord(i) = ratios.coord(i);    
    
    std::vector<int> fail(omp_get_max_threads(),0); //sometimes the fit func inversion can fail on noisy data
#pragma omp parallel for
    for(int j=first_sample;j<nsample;j++){
      dataSeries<double, double> rat_t_sample_j(1);
      rat_t_sample_j.coord(0) = ratios.coord(i);
      rat_t_sample_j.value(0) = ratios.value(i).sample(j);
      CostFunction costfunc(fiteffmass, rat_t_sample_j, sigma_j);
      MinimizerType fitter(costfunc, mlparams);
      FitMpiEffectiveMass::ParameterType p(0.5);
      fitter.fit(p);
      if(!fitter.hasConverged()) fail[omp_get_thread_num()] = 1;
      else effmass.value(i).sample(j) = *p;
    }
    int did_fail = 0; for(int ii=0;ii<fail.size();ii++) did_fail += fail[ii];
    if(did_fail){
      std::cout << "Failed to converge on one or more samples at coord " << effmass.coord(i) << std::endl;
      erase[i] = true;
      erase_required = true;
    }else{
      std::cout << effmass.coord(i) << " " << effmass.value(i) << std::endl;
    }
  }
  distributionPrint<jackknifeDistributionType>::printer(orig_printer);
    
  if(erase_required){
    int nkeep = 0; for(int i=0;i<ratios.size();i++) if(!erase[i]) nkeep++;
    jackknifeTimeSeriesType tokeep(nkeep);
    int elem = 0;
    for(int i=0;i<ratios.size();i++)
      if(!erase[i]){
	tokeep.coord(elem) = std::move(effmass.coord(i));
	tokeep.value(elem) = std::move(effmass.value(i));
	elem++;
      }
    return tokeep;
  }else return effmass;
}



#endif



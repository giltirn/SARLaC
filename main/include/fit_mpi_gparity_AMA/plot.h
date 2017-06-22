#ifndef _FIT_MPI_GPARITY_AMA_PLOT_H
#define _FIT_MPI_GPARITY_AMA_PLOT_H

template<typename DataSeriesType>
class dataTypeFilter{
  const DataSeriesType &series;
  std::vector<int> elems;
public:
  dataTypeFilter(const DataSeriesType &_series, const DataType type): series(_series){
    for(int i=0;i<series.size();i++)
      if(series.coord(i).type == type) elems.push_back(i);
    std::sort(elems.begin(), elems.end(),
	      [&](int a, int b) { return series.coord(a).t < series.coord(b).t; }
	      );
  }
  inline size_t size() const{ return elems.size(); }
  inline double coord(const int i) const{ return series.coord(elems[i]).t; }
  inline auto value(const int i)->decltype( series.value(elems[i]) ){ return series.value(elems[i]); }  
};


jackknifeTimeSeriesType effectiveMass(const jackknifeTimeSeriesType &data, const DataType type, const int Lt){
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
  typedef UncorrelatedChisqCostFunction<FitEffectiveMass, dataSeries<double,double> > CostFunction;
  typedef MarquardtLevenbergMinimizer<CostFunction> MinimizerType;
  FitEffectiveMass fiteffmass(2*Lt, type);
  
  MarquardtLevenbergParameters<typename CostFunction::CostType> mlparams;
  mlparams.verbose = false;
  mlparams.exit_on_convergence_fail = false;
  std::vector<double> sigma_j(1,1.);

  jackknifeTimeSeriesType effmass(ratios.size(), nsample);
  publicationPrint<> printer;
  printer << "Effective mass for type " << type << "\n";

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
      FitEffectiveMass::ParameterType p(0.5);
      fitter.fit(p);
      if(!fitter.hasConverged()) fail[omp_get_thread_num()] = 1;
      else effmass.value(i).sample(j) = *p;
    }
    int did_fail = 0; for(int i=0;i<fail.size();i++) did_fail += fail[i];
    if(did_fail){
      printer << "Failed to converge on one or more samples at coord " << effmass.coord(i) << std::endl;
      erase[i] = true;
      erase_required = true;
    }else{
      printer << effmass.coord(i) << " " << effmass.value(i) << std::endl;
    }
  }
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


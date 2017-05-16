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




template<typename DataSeriesType>
jackknifeTimeSeriesD effectiveMass(const DataSeriesType &data, const DataType type, const double Lt){
  const int nsample = data.value(0).size();
  dataTypeFilter<DataSeriesType> filter(data, type);
  dataSeries<double, jackknifeDistributionD> ratios(filter.size()-1);
  for(int i=0;i<filter.size()-1;i++){
    double t = filter.coord(i);
    assert(t == double(i));
    assert(filter.coord(i+1) == t+1);

    ratios.coord(i) = t;
    ratios.value(i) = filter.value(i)/filter.value(i+1);
  }
  typedef UncorrelatedChisqCostFunction<FitEffectiveMass, dataSeries<double,double> > CostFunction;
  typedef MarquardtLevenbergMinimizer<CostFunction> MinimizerType;
  FitEffectiveMass fiteffmass(Lt, type);
  
  MarquardtLevenbergParameters<typename CostFunction::CostType> mlparams;
  mlparams.verbose = false;
  std::vector<double> sigma_j(1,1.);

  jackknifeTimeSeriesD effmass(ratios.size(), nsample);
  publicationPrint printer;
  printer << "Effective mass for type " << toStr(type) << "\n";
  for(int i=0;i<ratios.size();i++){
    effmass.coord(i) = ratios.coord(i);
    
#pragma omp parallel for
    for(int j=0;j<nsample;j++){
      dataSeries<double, double> rat_t_sample_j(1);
      rat_t_sample_j.coord(0) = ratios.coord(i);
      rat_t_sample_j.value(0) = ratios.value(i).sample(j);
      CostFunction costfunc(fiteffmass, rat_t_sample_j, sigma_j);
      MinimizerType fitter(costfunc, mlparams);
      FitEffectiveMass::ParameterType p(0.5);
      fitter.fit(p);
      effmass.value(i).sample(j) = *p;
    }

    printer << effmass.coord(i) << " " << effmass.value(i) << std::endl;
  }
  return effmass;
}
  
#endif



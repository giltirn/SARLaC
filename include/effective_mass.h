#ifndef _EFFECTIVE_MASS_H_
#define _EFFECTIVE_MASS_H_

template<typename T>
class MLwrapper{
  T t;
public:
  ENABLE_GENERIC_ET(MLwrapper, MLwrapper<T>);
  //typedef MLwrapper<T> ET_tag;
  MLwrapper(): t(0.){}
  explicit MLwrapper(const T _t): t(_t){}
  inline T & operator()(const int i){ assert(i==0); return t; }
  inline const T &operator()(const int i) const{ assert(i==0); return t; }
  inline int size() const{ return 1; }
  inline std::string print() const{ return anyToStr<T>(t); }
  inline void resize(const int i){ if(i!=1) error_exit(std::cout << printType<MLwrapper<T> >() << " resize called with value " << i << " != 1\n"); }
  inline void zero(){ t=0.; }
  inline T& operator*(){ return t; }
  inline const T& operator*() const{ return t; }
};
template<typename T>
struct getElem<MLwrapper<T> >{
  static inline T& elem(MLwrapper<T> &v, const int i){ return *v; }
  static inline const T& elem(const MLwrapper<T> &v, const int i){ return *v; }    
  inline static int common_properties(const MLwrapper<T> &v){ return 1; }
};


//For fit functions with a single amplitude coefficient A*f(m,t) we can extract the effective mass by numerically inverting   C(t)/C(t+1) = f(m,t)/f(m,t+1)
template<typename FitFunc>
class Fit2ptEffectiveMass{
public:
  typedef typename FitFunc::ParameterType BaseParameterType;
  typedef typename FitFunc::ValueDerivativeType BaseDerivativeType;
  typedef double ValueType;
  typedef MLwrapper<double> ParameterType;
  typedef MLwrapper<double> ValueDerivativeType;
  typedef double GeneralizedCoordinate;
private:  
  FitFunc const* fitfunc;
  int params_mass_index;
  BaseParameterType base;
public:

  Fit2ptEffectiveMass(const FitFunc &ff, const BaseParameterType _base, const int pmi): fitfunc(&ff),  params_mass_index(pmi), base(_base){
    //assert(base.size() == 2);
  }

  inline double value(const double t, const ParameterType &params) const{    
    BaseParameterType p(base); p(params_mass_index) = *params;
    return fitfunc->value(t,p)/fitfunc->value(t+1,p);
  }
  
  ValueDerivativeType parameterDerivatives(const double t, const ParameterType &params) const{
    ValueDerivativeType yderivs;
    BaseParameterType p(base); p(params_mass_index) = *params;
    
    double value_t = fitfunc->value(t,p);
    BaseDerivativeType subderivs_t = fitfunc->parameterDerivatives(t,p);

    double value_tp1 = fitfunc->value(t+1,p);
    BaseDerivativeType subderivs_tp1 = fitfunc->parameterDerivatives(t+1,p);
 
    *yderivs = subderivs_t(params_mass_index)/value_tp1 - value_t/value_tp1/value_tp1 * subderivs_tp1(params_mass_index);
    return yderivs;
  }

  inline int Nparams() const{ return 1; }
};

//This is a zero degree of freedom fit to numerically invert FitEffMassFunc for the effective mass
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
  auto orig_printer = distributionPrint<jackknifeDistributionD>::printer();
  distributionPrint<jackknifeDistributionD>::printer(new publicationDistributionPrinter<jackknifeDistributionD>,false);  

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
    int did_fail = 0; for(int i=0;i<fail.size();i++) did_fail += fail[i];
    if(did_fail){
      std::cout << "Warning: Failed to converge on one or more samples at coord " << effmass.coord(i) << std::endl;
      erase[i] = true;
      erase_required = true;
    }else{
      jackknifeDistributionD y(nsample);
      jackknifeDistributionD yfit(nsample);
      jackknifeDistributionD resid(nsample);
      for(int s=0;s<nsample;s++){
	double t = edata.coord(i);
	y.sample(s) = edata.value(i).sample(s);
	MLwrapper<double> m(effmass.value(i).sample(s));
	yfit.sample(s) = fiteffmass.value(t,m);

	resid.sample(s) = (y.sample(s) - yfit.sample(s))/y.sample(s);
      }
      std::cout << "Effmass t="<<  effmass.coord(i) << " m=" << effmass.value(i) << " y=" << y << " yfit=" << yfit << " resid=" << resid << std::endl;
    }
  }
  distributionPrint<jackknifeDistributionD>::printer(orig_printer);
    
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
//Base should be a correctly setup parameter structure (needed so we can avoid requiring default constructors)
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


#endif

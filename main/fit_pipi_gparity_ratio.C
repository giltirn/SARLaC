#include<numeric_tensors.h>
#include<common_defs.h>
#include<parser.h>
#include<fit_wrapper.h>
#include<fit_wrapper_freeze.h>
#include<plot.h>
#include<effective_mass.h>

NumericVector<double> operator*(const NumericVector<double> &a, const NumericVector<double> &b){
  assert(a.size() == b.size());
  NumericVector<double> out(a.size());
  for(int i=0;i<a.size();i++) out(i) = a(i) * b(i);
  return out;
}


typedef correlationFunction<double,doubleJackknifeDistributionD> doubleJackCorrelationFunction;
typedef correlationFunction<double,jackknifeDistributionD> jackCorrelationFunction;

template<typename CorrFuncType>
typename CorrFuncType::ElementType computeRatio(const int t, const CorrFuncType &pipi, const CorrFuncType &pion_2pt){
  return typename CorrFuncType::ElementType(t,
					    (pipi.value(t+1) - pipi.value(t))/pion_2pt.value(t)/pion_2pt.value(t));
}

class FitRatio{
public:
  typedef double ValueType;
  typedef NumericVector<double> ParameterType;
  typedef NumericVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitRatio(){}
  
  ValueType value(const GeneralizedCoordinate t, const ParameterType &p) const{
    return p[0]*(1-p[1]*t);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs(2);
    yderivs[0] = 1-p[1]*t;
    yderivs[1] = -p[0]*t;
    return yderivs;
  }

  inline int Nparams() const{ return 2; }
};

class FitExactRatio{
  double T;
  double Lt;

  inline double Cpipi(const double t, const double dE, const double Epi) const{
    double Epipi = dE + 2*Epi;
    return exp(-Epipi*t) + exp(-Epipi*(T - t));
  }
  inline double dCpipi_by_ddE(const double t, const double dE, const double Epi) const{
    double Epipi = dE + 2*Epi;
    return -t*exp(-Epipi*t) - (T-t)*exp(-Epipi*(T - t));
  }
  inline double dCpipi_by_dEpi(const double t, const double dE, const double Epi) const{
    double Epipi = dE + 2*Epi;
    return -2*t*exp(-Epipi*t) - 2*(T-t)*exp(-Epipi*(T - t));
  }
  
  inline double Cpi_sq(const double t, double Epi) const{
    return exp(-2*Epi*t) + exp(-2*Epi*(Lt-t)) + 2*exp(-2*Epi*Lt);
  }
  inline double Cpi_sq_by_dEpi(const double t, double Epi) const{
    return -2*t*exp(-2*Epi*t) -2*(Lt-t)*exp(-2*Epi*(Lt-t)) - 4*Lt*exp(-2*Epi*Lt);
  }
public:
  typedef double ValueType;
  typedef NumericVector<double> ParameterType;
  typedef NumericVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitExactRatio(double _Lt, double _tsep_pipi): Lt(_Lt), T(_Lt - 2*_tsep_pipi){}

  ValueType value(const GeneralizedCoordinate t, const ParameterType &p) const{
    double Cpipi_tp1 = Cpipi(t+1,p(1),p(2));
    double Cpipi_t = Cpipi(t,p(1),p(2));
    double Cpi2 = Cpi_sq(t,p(2));

    return p(0) * ( Cpipi_tp1 - Cpipi_t )/Cpi2;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    double Cpipi_tp1 = Cpipi(t+1,p(1),p(2));
    double Cpipi_t = Cpipi(t,p(1),p(2));
    double Cpi2 = Cpi_sq(t,p(2));

    double dCpipi_tp1_by_ddE = dCpipi_by_ddE(t+1,p(1),p(2));
    double dCpipi_t_by_ddE = dCpipi_by_ddE(t,p(1),p(2));

    double dCpipi_tp1_by_dEpi = dCpipi_by_dEpi(t+1,p(1),p(2));
    double dCpipi_t_by_dEpi = dCpipi_by_dEpi(t,p(1),p(2));

    double dCpi2_by_dEpi = Cpi_sq_by_dEpi(t,p(2));
    
    ValueDerivativeType yderivs(2);
    yderivs[0] = ( Cpipi_tp1 - Cpipi_t )/Cpi2;
    yderivs[1] = p(0) * ( dCpipi_tp1_by_ddE - dCpipi_t_by_ddE )/Cpi2;
    yderivs[2] = p(0) * ( dCpipi_tp1_by_dEpi - dCpipi_t_by_dEpi )/Cpi2 - p(0)*( Cpipi_tp1 - Cpipi_t )/Cpi2/Cpi2 * dCpi2_by_dEpi;    
    return yderivs;
  }

  inline int Nparams() const{ return 3; } //A, dE, Epi
};



#define ARGS_MEMBERS \
  ( std::string, pipi_data_file ) \
  ( std::string, pion_2pt_data_file )	     \
  ( std::string, Epi_file )		     \
  ( int, Epi_idx )				     \
  ( int, Lt) \
  ( int, tsep_pipi ) \
  ( int, tmin) \
  ( int, tmax)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);



int main(const int argc, const char** argv){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);
  
  doubleJackCorrelationFunction pipi_data_dj;
  {
    HDF5reader reader(args.pipi_data_file);
    read(reader, pipi_data_dj, "data");
  }
  doubleJackCorrelationFunction pion_2pt_data_dj;
  {
    HDF5reader reader(args.pion_2pt_data_file);
    read(reader, pion_2pt_data_dj, "data");
  }

  jackknifeDistributionD Epi;
  {
    std::vector<jackknifeDistributionD> params;
    readParamsStandard(params, args.Epi_file);
    Epi = params[args.Epi_idx];
  }

  jackCorrelationFunction pipi_data_j(args.Lt, [&](const int t){ return typename jackCorrelationFunction::ElementType(t,pipi_data_dj.value(t).toJackknife()); } );
  jackCorrelationFunction pion_2pt_data_j(args.Lt, [&](const int t){ return typename jackCorrelationFunction::ElementType(t,pion_2pt_data_dj.value(t).toJackknife()); } );

  doubleJackCorrelationFunction ratio_dj(args.Lt-1, [&](const int t){ return computeRatio(t,pipi_data_dj,pion_2pt_data_dj); });
  jackCorrelationFunction ratio_j(args.Lt-1, [&](const int t){ return computeRatio(t,pipi_data_j,pion_2pt_data_j); });
  
  int nsample = pipi_data_dj.value(0).size();

  int nt_fit = args.tmax-args.tmin+1;
  
  doubleJackCorrelationFunction ratio_dj_inrange(nt_fit, [&](const int i){ return ratio_dj[args.tmin+i]; });
  jackCorrelationFunction ratio_j_inrange(nt_fit, [&](const int i){ return ratio_j[args.tmin+i]; });

  //#define FIT_LIN_APPROX
  #define FIT_EXACT
#ifdef FIT_LIN_APPROX
  typedef composeFitPolicy<FitFunc,standardFitFuncPolicy,correlatedFitPolicy>::type FitPolicies;
  typedef FitRatio FitFunc;
  FitFunc fitfunc;
#elif defined(FIT_EXACT)  
  typedef FitExactRatio FitFunc;
  typedef composeFitPolicy<FitFunc,frozenFitFuncPolicy,correlatedFitPolicy>::type FitPolicies;  
  FitFunc fitfunc(args.Lt, args.tsep_pipi);
#else
  #error "Invalid fit type"
#endif

    
  fitter<FitPolicies> fit;
  fit.importFitFunc(fitfunc);

#ifdef FIT_EXACT
  //Freeze in Epi
  typename FitPolicies::FitParameterDistribution freeze(nsample,typename FitFunc::ParameterType(3));
  for(int s=0;s<nsample;s++) freeze.sample(s)(2) = Epi.sample(s);
  fit.freeze(std::vector<int>({2}), freeze);
#endif
  
  importCostFunctionParameters<correlatedFitPolicy,FitPolicies> prepare(fit,ratio_dj_inrange);

  typename FitFunc::ParameterType guess(fitfunc.Nparams()); guess(0) = 1; guess(1) = 0;

  typename FitPolicies::FitParameterDistribution params(nsample, guess);
  jackknifeDistributionD chisq(nsample);
  jackknifeDistributionD chisq_per_dof(nsample);

  fit.fit(params, chisq, chisq_per_dof, ratio_j_inrange);

  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

  std::vector<jackknifeDistributionD> params_j(fitfunc.Nparams());
  for(int p=0;p<fitfunc.Nparams();p++){
    params_j[p] = jackknifeDistributionD(nsample, [&](const int s){ return params.sample(s)(p); });
    std::cout << "Params[" << p << "]: " << params_j[p] << std::endl;
  }
  
  jackknifeDistributionD Epipi = params_j[1] + 2.*Epi;

  std::cout << "Epi: " << Epi << std::endl;
  std::cout << "2*Epi: " << jackknifeDistributionD(2.*Epi) << std::endl;
  std::cout << "Epipi: " << Epipi << std::endl;

  {
    jackCorrelationFunction effdE(args.Lt-2);
    for(int t=0;t<args.Lt-2;t++){
      effdE.coord(t) = t;
      effdE.value(t) = (ratio_j.value(t+1) - ratio_j.value(t))/( ratio_j.value(t+1)*double(t) - ratio_j.value(t)*double(t+1) );
    }

    MatPlotLibScriptGenerate plot;
    typedef MatPlotLibScriptGenerate::kwargsType kwargsType;
    kwargsType kwargs;
    kwargs["color"] = "r";
    kwargs["alpha"] = 0.3;
    
    typedef DataSeriesAccessor<jackCorrelationFunction, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
    accessor acc(effdE);
    plot.plotData(acc);

    typedef BandRangeConstantDistributionValue<jackknifeDistributionD, DistributionPlotAccessor<jackknifeDistributionD> > band_accessor;
    band_accessor band(args.tmin, args.tmax, params_j[1]);
    plot.errorBand(band,kwargs);
        
    plot.setXaxisBounds(0, args.tmax+1);
    plot.setYaxisBounds(params_j[1].best()-10.*params_j[1].standardError(), params_j[1].best()+10.*params_j[1].standardError() );
    
    plot.write("effdE.py","effdE.pdf");    
  }


  {
    MatPlotLibScriptGenerate plot;
    typedef DataSeriesAccessor<jackCorrelationFunction, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
    accessor acc(ratio_j);
    plot.plotData(acc);

    double ymin, ymax;
    {
      typedef MatPlotLibScriptGenerate::kwargsType kwargsType;
      kwargsType kwargs;
      kwargs["color"] = "r";
      kwargs["alpha"] = 0.3;
      
      std::vector<double> x(nt_fit);
      std::vector<double> upper(nt_fit);
      std::vector<double> lower(nt_fit);
      for(int i=0;i<nt_fit;i++){
	double t = args.tmin + i;
	jackknifeDistributionD vfit(nsample, [&](const int s){ return fitfunc.value(t, params.sample(s)); });
	x[i] = t;
	double c = vfit.best();
	double e = vfit.standardError();
	if(i==0){
	  ymax = c + 10*e;
	}else if(i == nt_fit-1){
	  ymin = c - 10*e;
	}	
	upper[i] = c+e;
	lower[i] = c-e;
      }
      BandVectorAccessor band(x,upper,lower);
      plot.errorBand(band,kwargs);
    }
    plot.setYaxisBounds(ymin,ymax);
    plot.setXaxisBounds(0, args.tmax+1);
    plot.write("plot.py","plot.pdf");    
  }
  
  return 0;
}

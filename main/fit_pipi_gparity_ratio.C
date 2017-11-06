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
  
  //Params are A, m  
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

#define ARGS_MEMBERS \
  ( std::string, pipi_data_file ) \
  ( std::string, pion_2pt_data_file )	     \
  ( std::string, Epi_file )		     \
  ( int, Epi_idx )				     \
  ( int, Lt) \
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
  
  typedef FitRatio FitFunc;
  typedef composeFitPolicy<double,FitFunc,standardFitFuncPolicy,correlatedFitPolicy>::type FitPolicies;
  
  FitFunc fitfunc;
  fitter<FitPolicies> fit;
  fit.importFitFunc(fitfunc);

  importCostFunctionParameters<correlatedFitPolicy,FitPolicies> prepare(fit,ratio_dj_inrange);

  typename FitFunc::ParameterType guess(2); guess(0) = 1; guess(1) = 0;

  typename FitPolicies::jackknifeFitParameters params(nsample, guess);
  jackknifeDistributionD chisq(nsample);
  jackknifeDistributionD chisq_per_dof(nsample);

  fit.fit(params, chisq, chisq_per_dof, ratio_j_inrange);

  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

  std::vector<jackknifeDistributionD> params_j(2);
  for(int p=0;p<2;p++){
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

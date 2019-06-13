#ifndef _FIT_PIPI_GPARITY_PLOT_H_
#define _FIT_PIPI_GPARITY_PLOT_H_

#include<config.h>
#include<utils/macros.h>

#include<plot.h>
#include<containers/single_value_container.h>
#include<common.h>

#include "fit.h"

CPSFIT_START_NAMESPACE

//two-point effective energy assuming cosh form with optional constant subtraction
jackknifeCorrelationFunctionD twoPointEffectiveEnergy(const jackknifeCorrelationFunctionD &data_j,
						     const jackknifeDistributionD &fitted_Epipi,
						     const jackknifeDistributionD &fitted_constant,
						     const int Lt, const int tsep_pipi, const double Ascale, const double Cscale,
						     const bool subtract_constant){
  jackknifeCorrelationFunctionD data_j_mod(data_j);
  if(subtract_constant){
    for(int i=0;i<data_j_mod.size();i++) data_j_mod.value(i) = data_j_mod.value(i) - Cscale*fitted_constant;
  }
  FitCoshPlusConstant fitfunc(Lt, tsep_pipi, Ascale, Cscale);
  FitCoshPlusConstant::Params base(1,1,0); //amplitude irrelevant, E will be varied, set constant to 0
  std::cout << "Computing two-point effective energy" << std::endl;
  return effectiveMass2pt<jackknifeCorrelationFunctionD,FitCoshPlusConstant>(data_j,fitfunc,base,1,Lt);
}


template<typename FitFunc>
class Fit3ptPiPiEffectiveMass{
public:
  typedef typename FitFunc::ParameterType BaseParameterType;
  typedef typename FitFunc::ValueDerivativeType BaseDerivativeType;
  typedef double ValueType;
  typedef singleValueContainer<double> ParameterType;
  typedef singleValueContainer<double> ValueDerivativeType;
  typedef double GeneralizedCoordinate;
private:  
  FitFunc const* fitfunc;
  int params_mass_index;
  BaseParameterType base;
public:

  Fit3ptPiPiEffectiveMass(const FitFunc &ff, const BaseParameterType _base, const int pmi): fitfunc(&ff),  params_mass_index(pmi), base(_base){}

  //[C(t+1)-C(t)]/[C(t+2)-C(t+1)]
  
  inline double value(const double t, const ParameterType &params) const{    
    BaseParameterType p(base); p(params_mass_index) = *params;
    double value_t = fitfunc->value(t,p);
    double value_tp1 = fitfunc->value(t+1,p);
    double value_tp2 = fitfunc->value(t+2,p);

    return (value_tp1 - value_t)/(value_tp2 - value_tp1);
  }
  
  ValueDerivativeType parameterDerivatives(const double t, const ParameterType &params) const{
    ValueDerivativeType yderivs;
    BaseParameterType p(base); p(params_mass_index) = *params;
    
    double value_t = fitfunc->value(t,p);
    double deriv_t = fitfunc->parameterDerivatives(t,p)(params_mass_index);
    
    double value_tp1 = fitfunc->value(t+1,p);
    double deriv_tp1 = fitfunc->parameterDerivatives(t+1,p)(params_mass_index);
    
    double value_tp2 = fitfunc->value(t+2,p);
    double deriv_tp2 = fitfunc->parameterDerivatives(t+2,p)(params_mass_index);
    
    double value_num = value_tp1 - value_t;
    double value_den = value_tp2 - value_tp1;

    double deriv_num = deriv_tp1 - deriv_t;
    double deriv_den = deriv_tp2 - deriv_tp1;
    
    *yderivs = deriv_num/value_den - value_num/value_den/value_den*deriv_den;
    return yderivs;
  }

  inline int Nparams() const{ return 1; }
};

jackknifeCorrelationFunctionD threePointEffectiveEnergy(const jackknifeCorrelationFunctionD &data_j,
						       const jackknifeDistributionD &fitted_Epipi,
						       const int Lt, const int tsep_pipi, const double Ascale, const double Cscale){
  const int nsample = data_j.value(0).size();
  assert(data_j.size() == Lt);
  jackknifeCorrelationFunctionD ratios(Lt-2);
  for(int i=0;i<Lt-2;i++){
    double t = data_j.coord(i);
    assert(t == double(i));
    assert(data_j.coord(i+1) == t+1);
    assert(data_j.coord(i+2) == t+2);

    ratios.coord(i) = t;
    ratios.value(i) = (data_j.value(i+1) - data_j.value(i))/(data_j.value(i+2) - data_j.value(i+1));
  }
  typedef Fit3ptPiPiEffectiveMass<FitCoshPlusConstant> FitEffMass;

  FitCoshPlusConstant fitfunc(Lt, tsep_pipi, Ascale, Cscale);
  FitCoshPlusConstant::Params base(1,1,0); //amplitude and constant irrelevant, E will be varied
  std::cout << "Computing three-point effective energy" << std::endl;
  FitEffMass fiteffmass(fitfunc, base, 1);
  return fitEffectiveMass<jackknifeCorrelationFunctionD,FitEffMass>(ratios,fiteffmass);
}

  
jackknifeCorrelationFunctionD effectiveEnergy(const jackknifeCorrelationFunctionD &data_j,
					     const jackknifeDistributionD &fitted_Epipi,
					     const jackknifeDistributionD &fitted_constant,
					     const PiPiEffectiveEnergy effective_energy,
					     const int Lt, const int tsep_pipi, const double Ascale, const double Cscale){
  switch(effective_energy){
  case PiPiEffectiveEnergy::TwoPoint:
    return twoPointEffectiveEnergy(data_j,fitted_Epipi,fitted_constant, Lt, tsep_pipi, Ascale, Cscale, false);
  case PiPiEffectiveEnergy::TwoPointSubConstant:
    return twoPointEffectiveEnergy(data_j,fitted_Epipi,fitted_constant, Lt, tsep_pipi, Ascale, Cscale, true);
  case PiPiEffectiveEnergy::ThreePoint:
    return threePointEffectiveEnergy(data_j,fitted_Epipi, Lt, tsep_pipi, Ascale, Cscale);
  default:
    error_exit(std::cout << "Unknown effective energy type "<< effective_energy <<std::endl);
  }
}

void plot(const jackknifeCorrelationFunctionD &data_j,
	  const jackknifeDistributionD &fitted_Epipi,
	  const jackknifeDistributionD &fitted_constant,
	  const int fit_t_min, const int fit_t_max,
	  const PiPiEffectiveEnergy effective_energy,
	  const int Lt, const int tsep_pipi, const double Ascale, const double Cscale){
  jackknifeCorrelationFunctionD E_eff = effectiveEnergy(data_j,fitted_Epipi,fitted_constant, effective_energy, Lt, tsep_pipi, Ascale, Cscale);

  std::cout << "Effective energy:\n" << E_eff << std::endl;

  MatPlotLibScriptGenerate plotter;
  typedef MatPlotLibScriptGenerate::handleType handleType;
  MatPlotLibScriptGenerate::kwargsType plot_args;
  plot_args["color"] = "r";

  //Plot the effective energy
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD,ScalarCoordinateAccessor<double>,DistributionPlotAccessor<jackknifeDistributionD> > accessor;
  accessor effenergy_accessor(E_eff);

  handleType plotdata = plotter.plotData(effenergy_accessor,plot_args);
  
  //Plot the fitted energy across the fit range
  struct fitBand{
    double u;
    double l;
    int tmin;
    int tmax;    
    fitBand(double uu, double ll, int ttmin, int ttmax): u(uu),l(ll),tmin(ttmin),tmax(ttmax){}
    inline int size() const{ return tmax-tmin+1; }    
    inline double x(const int i) const{ return tmin+i; }
    inline double upper(const int i) const{ return u; }
    inline double lower(const int i) const{ return l; }
  };

  double se = fitted_Epipi.standardError();
  fitBand band(fitted_Epipi.best() + se, fitted_Epipi.best() - se, fit_t_min, fit_t_max);
  plot_args["alpha"] = 0.2;
  handleType fband = plotter.errorBand(band,plot_args);

  plotter.setXaxisBounds(-0.1,10.9);
  plotter.setYaxisBounds(fitted_Epipi.best()-10*se, fitted_Epipi.best()+10*se);
  plotter.setXlabel("$t$");
  plotter.setYlabel("$E_{\\pi\\pi}^{\\rm eff}$");
  
  plotter.write("effective_energy.py","effective_energy.pdf");
}

CPSFIT_END_NAMESPACE

#endif

#ifndef _CPSFIT_FITFUNC_BOUNDED_H_
#define _CPSFIT_FITFUNC_BOUNDED_H_

//This fit function allows us to constrain fit parameters to have some minimum and/or maximum value
//This is achieved by mapping the external, bounded parameter to an internal unbounded parameter under which the chi^2 is minimized

#include<config.h>
#include<utils/macros.h>
#include<containers/parameter_vector.h>
#include<fit/param_bounds.h>

CPSFIT_START_NAMESPACE

template<typename FitFunc>
class BoundedFitFunc{
public:
  typedef typename FitFunc::ValueType ValueType;
  typedef typename FitFunc::ParameterType ParameterType;
  typedef typename FitFunc::ValueDerivativeType ValueDerivativeType;
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;

private:

  const FitFunc &fitfunc;
  std::vector<boundedParameterTransform> trans;

public:
  BoundedFitFunc(const FitFunc &_fitfunc): fitfunc(_fitfunc){}

  inline void setBound(const boundedParameterTransform &tr){
    for(int t=0;t<trans.size();t++) if(tr.param == trans[t].param) error_exit(std::cout << "BoundedFitFunc::setBound cannot have multiple bounds on same parameter!\n");
    trans.push_back(tr); 
  }

  inline void setBound(int param, const ParameterBound bound, double min, double max){ //min or max will be ignored as appropriate if using Max or Min bounds
    setBound(boundedParameterTransform(param, bound, min, max));
  }

  inline ValueType value(const GeneralizedCoordinate &coord, const ParameterType &params_ub) const{
    ParameterType params_b(params_ub);
    for(int t=0;t<trans.size();t++) trans[t].mapUnboundedToBounded(params_b);
    return fitfunc.value(coord, params_b);
  }

  //Deritive wrt *unbounded* parameters
  //d/du = d/db * db/du
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &coord, const ParameterType &params_ub) const{
    ParameterType params_b(params_ub);
    for(int t=0;t<trans.size();t++) trans[t].mapUnboundedToBounded(params_b);

    ValueDerivativeType derivs_b = fitfunc.parameterDerivatives(coord, params_b);
    
    for(int t=0;t<trans.size();t++) trans[t].jacobian(derivs_b, params_ub); //convert derivatives bounded -> unbounded. Note takes in unbounded variable as 2nd param
    
    return derivs_b;
  }
  
  inline int Nparams() const{ return fitfunc.Nparams(); }
};

CPSFIT_END_NAMESPACE

#endif

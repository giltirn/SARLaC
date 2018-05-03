#ifndef _GLOBALFIT_FITFUNC_GENERAL_H_
#define _GLOBALFIT_FITFUNC_GENERAL_H_

template<typename GlobalFitPolicies>
class GlobalFit: public GlobalFitPolicies{
public:
  typedef double ValueType;
  typedef typename GlobalFitPolicies::ParameterType ParameterType;
  typedef typename GlobalFitPolicies::ParameterType ValueDerivativeType; //derivative wrt parameters
  typedef DataParams GeneralizedCoordinate;

private:

  template<typename Binding, typename FitFunc>
  ValueType value_int(const FitFunc &fitfunc, const GeneralizedCoordinate &t, const ParameterType &p) const{
    typename Binding::subsetType psub;
    ParameterType::template supersetToSubset<Binding>(psub,p,t);
    return fitfunc.value(t,psub);
  }
  template<typename Binding, typename FitFunc>
  ValueDerivativeType parameterDerivatives_int(const FitFunc &fitfunc, const GeneralizedCoordinate &t, const ParameterType &p) const{
    typename Binding::subsetType psub;
    ParameterType::template supersetToSubset<Binding>(psub,p,t);
    typename Binding::subsetType dsub = fitfunc.parameterDerivatives(t,psub);
    ValueDerivativeType out(p); out.zero();
    ParameterType::template subsetToSuperset<Binding>(out,dsub,t);
    return out;
  }

public:
  GlobalFit(const typename GlobalFitPolicies::fitFuncSetupType &setup): GlobalFitPolicies(setup){
  }

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    switch(t.type){
    case mpi2:
      return value_int<typename GlobalFitPolicies::Bindings_mpi2, decltype(this->fit_mpi2)>(this->fit_mpi2, t,p);
    case mK2:
      return value_int<typename GlobalFitPolicies::Bindings_mK2, decltype(this->fit_mK2)>(this->fit_mK2, t,p);
    case mOmega:
      return value_int<typename GlobalFitPolicies::Bindings_mOmega, decltype(this->fit_mOmega)>(this->fit_mOmega, t,p);
    default:
      error_exit(std::cout << "GlobalFit::value unknown type " << t.type << std::endl);
    }
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    switch(t.type){
    case mpi2:
      return parameterDerivatives_int<typename GlobalFitPolicies::Bindings_mpi2, decltype(this->fit_mpi2)>(this->fit_mpi2, t,p);
    case mK2:
      return parameterDerivatives_int<typename GlobalFitPolicies::Bindings_mK2, decltype(this->fit_mK2)>(this->fit_mK2, t,p);
    case mOmega:
      return parameterDerivatives_int<typename GlobalFitPolicies::Bindings_mOmega, decltype(this->fit_mOmega)>(this->fit_mOmega, t,p);
    default:
      error_exit(std::cout << "AnalyticGlobalFit::value unknown type " << t.type << std::endl);
    }
  }

  inline int Nparams() const{ return this->GlobalFitPolicies::Nparams(); }

  ParameterType guess() const{ return this->GlobalFitPolicies::guess(); }
};


struct globalFitTypes{
  typedef superJackknifeDistribution<double> DistributionType;
  typedef correlationFunction<superJackknifeDistribution<DataParams>,DistributionType> CorrelationFunctionDistribution;
};


#endif

#pragma once

class FitConstantFrozen{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return 0.1;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return ValueDerivativeType(0);
  }

  inline int Nparams() const{ return 0; }

  ParameterType guess() const{ return ParameterType(0); }
};

class FitExpP{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p[0] * ::exp(-p[1]*t);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs(2);
    yderivs[0]= ::exp(-p[1]*t);
    yderivs[1] = -p[0] * t* ::exp(-p[1]*t);
    return yderivs;
  }

  inline int Nparams() const{ return 2; }

  ParameterType guess() const{ return ParameterType({1.0,0.5}); }
};

class FitExpPlusConst{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p[0] * ::exp(-p[1]*t) + p[2];
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs(3);
    yderivs[0]= ::exp(-p[1]*t);
    yderivs[1] = -p[0] * t* ::exp(-p[1]*t);
    yderivs[2] = 1.;
    return yderivs;
  }

  inline int Nparams() const{ return 3; }

  ParameterType guess() const{ return ParameterType({1.0,0.5,0.}); }
};

#define MEMBERS (std::vector<int>, params)(std::vector<double>, values)
struct FreezeArgs{
  GENERATE_MEMBERS(MEMBERS); 
  FreezeArgs(): params(1,0),values(1,0.){  }
};
GENERATE_PARSER( FreezeArgs, MEMBERS);
#undef MEMBERS

template<typename F>
inline std::unique_ptr<genericFitFuncBase> enwrap(const F &ff, const FreezeArgs &freeze, bool do_freeze){
  auto ptr = new genericFitFuncFreezeWrapper<F>(ff, parameterVector<double>(ff.Nparams(),0.));
  if(do_freeze){
    std::vector<bool> fparams(ff.Nparams(),false);
    for(int f=0;f<freeze.params.size();f++){
      int p = freeze.params[f];
      fparams[p] = true;
      ptr->psetup(p) = freeze.values[f];
    }
    ptr->freezeParams(fparams);
  }
  return std::unique_ptr<genericFitFuncBase>(ptr);
}

inline std::unique_ptr<genericFitFuncBase> fitFuncFactory(FitFuncType type, const std::string &freeze_param_file = ""){
  FreezeArgs freeze; bool do_freeze(false);
  if(freeze_param_file.size()){
    parseOrTemplate(freeze, freeze_param_file, "freeze_params.args");
    do_freeze=true;
  }

  if(type == FitFuncType::FConstant){
    FitConstant fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FLinear){
    FitFuncLinearMultiDim<double,double,1> fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FPoly2){
    FitFuncLinearMultiDim<double,double,2> fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FPoly3){
    FitFuncLinearMultiDim<double,double,3> fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FPoly4){
    FitFuncLinearMultiDim<double,double,4> fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FPoly5){
    FitFuncLinearMultiDim<double,double,5> fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FConstantFrozen){
    FitConstantFrozen fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FExp){
    FitExpP fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else if(type == FitFuncType::FExpPlusConst){
    FitExpPlusConst fitfunc; return enwrap<decltype(fitfunc)>(fitfunc,freeze,do_freeze);
  }else{
    error_exit(std::cout << "Invalid fit function" << std::endl);
  }
}





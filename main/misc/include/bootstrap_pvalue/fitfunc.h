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

std::unique_ptr<genericFitFuncBase> fitFuncFactory(FitFuncType type){
  if(type == FitFuncType::FConstant){
    FitConstant fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitConstant>(fitfunc));
  }else if(type == FitFuncType::FLinear){
    FitFuncLinearMultiDim<double,double,1> fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitFuncLinearMultiDim<double,double,1> >(fitfunc));
  }else if(type == FitFuncType::FPoly2){
    FitFuncLinearMultiDim<double,double,2> fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitFuncLinearMultiDim<double,double,2> >(fitfunc));
  }else if(type == FitFuncType::FPoly3){
    FitFuncLinearMultiDim<double,double,3> fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitFuncLinearMultiDim<double,double,3> >(fitfunc));
  }else if(type == FitFuncType::FPoly4){
    FitFuncLinearMultiDim<double,double,4> fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitFuncLinearMultiDim<double,double,4> >(fitfunc));
  }else if(type == FitFuncType::FPoly5){
    FitFuncLinearMultiDim<double,double,5> fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitFuncLinearMultiDim<double,double,5> >(fitfunc));
  }else if(type == FitFuncType::FConstantFrozen){
    FitConstantFrozen fitfunc; 
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitConstantFrozen>(fitfunc));
  }else{
    error_exit(std::cout << "Invalid fit function" << std::endl);
  }
}

#pragma once

std::unique_ptr<genericFitFuncBase> fitFuncFactory(FitFuncType type){
  if(type == FitFuncType::FConstant){
    FitConstant fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitConstant>(fitfunc));
  }else if(type == FitFuncType::FLinear){
    FitFuncLinearMultiDim<double,double,1> fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitFuncLinearMultiDim<double,double,1> >(fitfunc));
  }else{
    error_exit(std::cout << "Invalid fit function" << std::endl);
  }
}

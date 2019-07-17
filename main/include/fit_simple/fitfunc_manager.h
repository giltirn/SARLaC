#ifndef _FIT_SIMPLE_FITFUNC_MANAGER_H_
#define _FIT_SIMPLE_FITFUNC_MANAGER_H_

#include "plot.h"

typedef parameterVector<double> parameterVectorD;

struct FitFuncManagerBase{
  int Lt;
  int t_min;
  int t_max;

  struct Options{
    bool load_guess;
    std::string guess_file;
    Options(): load_guess(false){}
  };

  Options opt;

  FitFuncManagerBase(const int Lt, const int t_min, const int t_max, const Options &opt = Options()): Lt(Lt), t_min(t_min), t_max(t_max), opt(opt){}

  virtual genericFitFuncBase const* getFitFunc() const = 0;
  virtual parameterVectorD getGuess() const = 0;
  virtual void plot(const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<parameterVectorD> &params) const = 0;
  virtual void plot(const bootstrapCorrelationFunctionD &data_j, const bootstrapDistribution<parameterVectorD> &params) const = 0;

  template<typename FF>
  parameterVectorD getGuessBase(const FF &ff) const{
   typename FF::ParameterType guess = ff.guess();
    if(opt.load_guess)
      parse(guess,opt.guess_file);
    return pconvert<parameterVectorD, typename FF::ParameterType>(guess);
  }
  
  virtual ~FitFuncManagerBase(){}
};

template<typename HyperbolicFitFunc>
struct FitFuncHyperbolicManager: public FitFuncManagerBase{
  typedef HyperbolicFitFunc FitFunc;
  simpleFitFuncWrapper<FitFunc> fitfunc;

  FitFuncHyperbolicManager(const int Lt, const int t_min, const int t_max, const FitFuncManagerBase::Options &opt = FitFuncManagerBase::Options()):
    FitFuncManagerBase(Lt, t_min, t_max, opt), fitfunc(Lt){}
    
  genericFitFuncBase const* getFitFunc() const{ return (genericFitFuncBase const*)&fitfunc; }

  parameterVectorD getGuess() const{ 
    return this->getGuessBase(fitfunc.fitfunc);
  }

  void plot(const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<parameterVectorD> &params) const{
    plotEffectiveMass(fitfunc.fitfunc,data_j,pconvert<typename FitFunc::ParameterType,parameterVectorD>(params),1, this->Lt, this->t_min, this->t_max);
  }
  void plot(const bootstrapCorrelationFunctionD &data_j, const bootstrapDistribution<parameterVectorD> &params) const{
    plotEffectiveMass(fitfunc.fitfunc,data_j,pconvert<typename FitFunc::ParameterType,parameterVectorD>(params),1, this->Lt, this->t_min, this->t_max);
  }
};

struct FitExpManager: public FitFuncManagerBase{
  typedef FitExp FitFunc;
  simpleFitFuncWrapper<FitFunc> fitfunc;

  FitExpManager(const int Lt, const int t_min, const int t_max, const FitFuncManagerBase::Options &opt = FitFuncManagerBase::Options()): FitFuncManagerBase(Lt, t_min, t_max, opt), fitfunc(FitFunc()){}

  genericFitFuncBase const* getFitFunc() const{ return (genericFitFuncBase const*)&fitfunc; }

  parameterVectorD getGuess() const{ 
    return this->getGuessBase(fitfunc.fitfunc);
  }

  void plot(const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<parameterVectorD> &params) const{
    plotEffectiveMass(fitfunc.fitfunc,data_j,pconvert<typename FitFunc::ParameterType,parameterVectorD>(params),1, this->Lt, this->t_min, this->t_max);
  }
  void plot(const bootstrapCorrelationFunctionD &data_j, const bootstrapDistribution<parameterVectorD> &params) const{
    plotEffectiveMass(fitfunc.fitfunc,data_j,pconvert<typename FitFunc::ParameterType,parameterVectorD>(params),1, this->Lt, this->t_min, this->t_max);
  }
};


struct FitConstantManager: public FitFuncManagerBase{
  typedef FitConstant FitFunc;
  simpleFitFuncWrapper<FitFunc> fitfunc;

  FitConstantManager(const int Lt, const int t_min, const int t_max, const FitFuncManagerBase::Options &opt = FitFuncManagerBase::Options()): FitFuncManagerBase(Lt, t_min, t_max, opt), fitfunc(FitFunc()){}

  genericFitFuncBase const* getFitFunc() const{ return (genericFitFuncBase const*)&fitfunc; }

  parameterVectorD getGuess() const{ 
    return this->getGuessBase(fitfunc.fitfunc);
  }

  void plot(const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<parameterVectorD> &params) const{
    plotRaw(fitfunc.fitfunc,data_j,pconvert<typename FitFunc::ParameterType,parameterVectorD>(params), this->Lt, this->t_min, this->t_max);
  }
  void plot(const bootstrapCorrelationFunctionD &data_j, const bootstrapDistribution<parameterVectorD> &params) const{
    plotRaw(fitfunc.fitfunc,data_j,pconvert<typename FitFunc::ParameterType,parameterVectorD>(params), this->Lt, this->t_min, this->t_max);
  }


};

struct FitTwoStateCoshManager: public FitFuncManagerBase{
  typedef FitTwoStateCosh FitFunc;
  simpleFitFuncWrapper<FitFunc> fitfunc;

  FitTwoStateCoshManager(const int Lt, const int t_min, const int t_max, const FitFuncManagerBase::Options &opt = FitFuncManagerBase::Options()): FitFuncManagerBase(Lt, t_min, t_max, opt), fitfunc(FitFunc(Lt)){}

  genericFitFuncBase const* getFitFunc() const{ return (genericFitFuncBase const*)&fitfunc; }

  parameterVectorD getGuess() const{ 
    return this->getGuessBase(fitfunc.fitfunc);
  }

  void plot(const jackknifeCorrelationFunctionD &data_j, const jackknifeDistribution<parameterVectorD> &params) const{
    FitCosh fcosh(this->Lt);
    plotTwoStateEffectiveMass<FitTwoStateCosh, FitCosh>(data_j, pconvert<typename FitFunc::ParameterType,parameterVectorD>(params), fcosh, fitfunc.fitfunc, 1, this->Lt, this->t_min, this->t_max); 
  }
  void plot(const bootstrapCorrelationFunctionD &data_j, const bootstrapDistribution<parameterVectorD> &params) const{
    FitCosh fcosh(this->Lt);
    plotTwoStateEffectiveMass<FitTwoStateCosh, FitCosh>(data_j, pconvert<typename FitFunc::ParameterType,parameterVectorD>(params), fcosh, fitfunc.fitfunc, 1, this->Lt, this->t_min, this->t_max); 
  }


};

std::unique_ptr< FitFuncManagerBase > getFitFuncManager(FitFuncType fitfunc, const int Lt, const int t_min, const int t_max, const FitFuncManagerBase::Options &opt = FitFuncManagerBase::Options()){
  std::unique_ptr< FitFuncManagerBase > fitfunc_manager;

  switch(fitfunc){
  case FitFuncType::FCosh:
    fitfunc_manager.reset(new FitFuncHyperbolicManager<FitCosh>(Lt,t_min,t_max,opt)); break;
  case FitFuncType::FSinh:
    fitfunc_manager.reset(new FitFuncHyperbolicManager<FitSinh>(Lt,t_min,t_max,opt)); break;
  case FitFuncType::FExp:
    fitfunc_manager.reset(new FitExpManager(Lt,t_min,t_max,opt)); break;
  case FitFuncType::FConstant:
    fitfunc_manager.reset(new FitConstantManager(Lt,t_min,t_max,opt)); break;
  case FitFuncType::FTwoStateCosh:
    fitfunc_manager.reset(new FitTwoStateCoshManager(Lt,t_min,t_max,opt)); break;
  default:
    error_exit(std::cout << "fit: Invalid fitfunc " << fitfunc << std::endl);
  }
  
  return fitfunc_manager;
}



#endif

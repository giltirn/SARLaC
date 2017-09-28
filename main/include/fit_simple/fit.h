#ifndef _FIT_SIMPLE_FIT_H_
#define _FIT_SIMPLE_FIT_H_

inline int getMassParamIdx(const FitFuncType type){
  switch(type){
  case FCosh:
  case FSinh:
  case FExp:
    return 1;
  default:
    error_exit(std::cout << "getMassParamIdx unknown fit func " << type << std::endl);
  }
}
template<typename FitFunc>
FitFunc* getFitFunc(const Args &args){ assert(0); }

template<>
FitCosh* getFitFunc<FitCosh>(const Args &args){ return new FitCosh(args.Lt); }

template<>
FitSinh* getFitFunc<FitSinh>(const Args &args){ return new FitSinh(args.Lt); }

template<>
FitExp* getFitFunc<FitExp>(const Args &args){ return new FitExp; }

template<typename FitFunc, template<typename> class CostFunctionPolicy>
void fitSpecFFcorr(const rawDataCorrelationFunctionD &data, const Args &args);

template<typename FitFunc>
inline void fitSpecFF(const rawDataCorrelationFunctionD &data, const Args &args, const CMDline &cmdline){
  return args.correlated ? fitSpecFFcorr<FitFunc,correlatedFitPolicy>(data,args,cmdline) : fitSpecFFcorr<FitFunc,uncorrelatedFitPolicy>(data,args,cmdline);
};
inline void fit(const rawDataCorrelationFunctionD &data, const Args &args, const CMDline &cmdline){
  switch(args.fitfunc){
  case FCosh:
    return fitSpecFF<FitCosh>(data,args,cmdline);
  case FSinh:
    return fitSpecFF<FitSinh>(data,args,cmdline);
  case FExp:
    return fitSpecFF<FitExp>(data,args,cmdline);    
  default:
    error_exit(std::cout << "fit: Invalid fitfunc " << args.fitfunc << std::endl);
  };
}

#endif

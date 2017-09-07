#ifndef _FIT_SIMPLE_FIT_H_
#define _FIT_SIMPLE_FIT_H_


template<typename FitPolicies, template<typename> class CostFunctionPolicy>
struct prepareFitter{};

template<typename FitPolicies>
struct prepareFitter<FitPolicies,correlatedFitPolicy>{
  NumericVector<jackknifeDistributionD> sigma;
  NumericSquareMatrix<jackknifeDistributionD> inv_corr;

  prepareFitter(fitter<FitPolicies> &fitter, const doubleJackknifeCorrelationFunctionD &data_dj){
    const int nsample = data_dj.value(0).size();
    const int nt_fit = data_dj.size();    
    NumericSquareMatrix<jackknifeDistributionD> cov(nt_fit);
    sigma.resize(nt_fit);
    for(int i=0;i<nt_fit;i++){
      cov(i,i) = doubleJackknifeDistributionD::covariance(data_dj.value(i),  data_dj.value(i));
      sigma(i) = sqrt(cov(i,i));

      for(int j=i+1;j<nt_fit;j++)
	cov(i,j) = cov(j,i) = doubleJackknifeDistributionD::covariance(data_dj.value(i),  data_dj.value(j));
    }

    NumericSquareMatrix<jackknifeDistributionD> corr(nt_fit);
    
    for(int i=0;i<nt_fit;i++){
      corr(i,i) = jackknifeDistributionD(nsample,1.);
      
      for(int j=i+1;j<nt_fit;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
    }
    
    inv_corr.resize(nt_fit, jackknifeDistributionD(nsample));
    svd_inverse(inv_corr, corr);
    std::cout << "Correlation matrix:\n" << corr << std::endl;
    std::cout << "Inverse correlation matrix:\n" << inv_corr << std::endl;  

    fitter.importCostFunctionParameters(inv_corr,sigma);
  }

};

template<typename FitPolicies>
struct prepareFitter<FitPolicies,uncorrelatedFitPolicy>{
  NumericVector<jackknifeDistributionD> sigma;

  prepareFitter(fitter<FitPolicies> &fitter, const doubleJackknifeCorrelationFunctionD &data_dj){
    int nt_fit = data_dj.size();    
    sigma.resize(nt_fit);
    for(int i=0;i<nt_fit;i++){
      sigma(i) = sqrt(doubleJackknifeDistributionD::covariance(data_dj.value(i),  data_dj.value(i)));
    }
    fitter.importCostFunctionParameters(sigma);
  }

};

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
inline void fitSpecFF(const rawDataCorrelationFunctionD &data, const Args &args){
  return args.correlated ? fitSpecFFcorr<FitFunc,correlatedFitPolicy>(data,args) : fitSpecFFcorr<FitFunc,uncorrelatedFitPolicy>(data,args);
};
inline void fit(const rawDataCorrelationFunctionD &data, const Args &args){
  switch(args.fitfunc){
  case FCosh:
    return fitSpecFF<FitCosh>(data,args);
  case FSinh:
    return fitSpecFF<FitSinh>(data,args);
  case FExp:
    return fitSpecFF<FitExp>(data,args);    
  default:
    error_exit(std::cout << "fit: Invalid fitfunc " << args.fitfunc << std::endl);
  };
}

#endif

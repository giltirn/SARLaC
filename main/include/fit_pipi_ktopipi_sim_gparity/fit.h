#ifndef _FIT_PIPI_KTOPIPI_SIM_GPARITY_FIT_H_
#define _FIT_PIPI_KTOPIPI_SIM_GPARITY_FIT_H_
#include<plot.h>
#include<gsl_eigensolve.h>

template<typename DistributionType>
class DistributionSampleAccessor{
  const DistributionType &d;
public:
  DistributionSampleAccessor(const DistributionType &_d): d(_d){}

  inline double x(const int i) const{ return i; }
  inline double y(const int i) const{ return iterate<DistributionType>::at(i,d); }
  inline double dxm(const int i) const{ return 0; }
  inline double dxp(const int i) const{ return 0; }
  inline double dym(const int i) const{ return 0; }
  inline double dyp(const int i) const{ return 0; }

  inline double upper(const int i) const{ return y(i); }
  inline double lower(const int i) const{ return y(i); }
  
  inline int size() const{ return iterate<DistributionType>::size(d); }
};

void analyzeCorrelationMatrix(const NumericSquareMatrix<jackknifeDistributionD> &corr){
  //Examine the eigenvalues and eigenvectors of the correlation matrix
  std::vector<NumericVector<jackknifeDistributionD> > evecs;
  std::vector<jackknifeDistributionD> evals;
  symmetricMatrixEigensolve(evecs,evals,corr);

  const int nsample = evals[0].size();
  const int nev = evals.size();
  jackknifeDistributionD one(nsample, 1.);
  std::vector<jackknifeDistributionD> inv_evals(nev);
  for(int i=0;i<nev;i++) inv_evals[i] = one/evals[i];
  
  std::cout << "Correlation matrix has size " << corr.size() << " and eigenvalues:\n";
  for(int i=0;i<nev;i++) std::cout << i << " " << evals[i] << std::endl;

  typedef DistributionArrayAccessor<std::vector<jackknifeDistributionD>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
  typedef DistributionSampleAccessor<jackknifeDistributionD> sample_accessor;
  {
    MatPlotLibScriptGenerate plt;
    accessor v(evals);
    plt.plotData(v);
    accessor w(inv_evals);
    plt.plotData(w);
    plt.write("corrmat_evals.py", "corrmat_evals.pdf");
  }
  {
    MatPlotLibScriptGenerate plt;
    for(int i=0;i<nev;i++){
      sample_accessor acc(evals[i]);
      plt.plotData(acc);
    }
    plt.write("corrmat_evals_samples.py", "corrmat_evals_samples.pdf");
  }
    
#ifdef HAVE_HDF5
  writeParamsStandard(evals,"corrmat_evals.hdf5");
  {
    std::vector<std::vector<jackknifeDistributionD> > tmp(nev, std::vector<jackknifeDistributionD>(nev));
    for(int i=0;i<nev;i++) for(int j=0;j<nev;j++) tmp[i][j] = evecs[i](j);
    writeParamsStandard(tmp,"corrmat_evecs.hdf5");
  }
#endif
  
}


typedef FitTwoPointThreePointSim FitFunc;

struct UncorrelatedFit;
struct CorrelatedFit;
struct PartiallyCorrelatedFit;

template<typename CorrelationStatusStruct>
struct CorrelationPolicy{};

template<>
struct CorrelationPolicy<UncorrelatedFit>{
  typedef typename composeFitPolicy<amplitudeDataCoordSim, FitFunc, frozenFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;
  std::unique_ptr<importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies> > prep;
    
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    prep.reset(new importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies>(fitter, data_combined_dj));
  }
};
template<>
struct CorrelationPolicy<CorrelatedFit>{
  typedef typename composeFitPolicy<amplitudeDataCoordSim, FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  std::unique_ptr<importCostFunctionParameters<correlatedFitPolicy,FitPolicies> > prep;
    
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    prep.reset(new importCostFunctionParameters<correlatedFitPolicy,FitPolicies>(fitter, data_combined_dj));
    analyzeCorrelationMatrix(prep->corr);
  }
};
template<>
struct CorrelationPolicy<PartiallyCorrelatedFit>{
  typedef typename composeFitPolicy<amplitudeDataCoordSim, FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  typedef FitPolicies::jackknifeMatrix jackknifeMatrix;
  typedef FitPolicies::jackknifeVector jackknifeVector;

  jackknifeMatrix corr;
  jackknifeMatrix inv_corr;
  jackknifeVector sigma;
  
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    //Only compute the correlations for data of the same type idx (pipi, Q1', Q2', Q3', ....)
    const int ndata = data_combined_dj.size();
    const int nsample = data_combined_dj.value(0).size();    

    sigma.resize(ndata);   
    
    jackknifeMatrix cov(ndata, jackknifeDistributionD(nsample,0.));
    for(int i=0;i<ndata;i++){
      cov(i,i) = doubleJackknifeDistributionD::covariance(data_combined_dj.value(i),data_combined_dj.value(i));
      sigma(i) = sqrt(cov(i,i));
      
      for(int j=i+1;j<ndata;j++){
	if(data_combined_dj.coord(i).idx == data_combined_dj.coord(j).idx)
	  cov(i,j) = cov(j,i) = doubleJackknifeDistributionD::covariance(data_combined_dj.value(i),data_combined_dj.value(j));
      }
    }

    corr = jackknifeMatrix(ndata);
    
    for(int i=0;i<ndata;i++){
      corr(i,i) = jackknifeDistributionD(nsample,1.);
      
      for(int j=i+1;j<ndata;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
    }
    
    inv_corr.resize(ndata, jackknifeDistributionD(nsample));
    svd_inverse(inv_corr, corr);

    //Test the quality of the inverse
    jackknifeMatrix test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistributionD(nsample,1.0);    
    std::cout << "|CorrMat * CorrMat^{-1} - 1|^2 = " << mod2(test) << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
    analyzeCorrelationMatrix(corr);
  }
};



template<typename CorrelationStatusStruct>
struct fitSpec{
  static void fit(const FitFunc &fitfunc,
		  const jackAmplitudeSimCorrelationFunction &data_combined_j,
		  const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,
		  const Args &args){
    typedef typename CorrelationPolicy<CorrelationStatusStruct>::FitPolicies FitPolicies;
    CorrelationPolicy<CorrelationStatusStruct> prepare;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);
    prepare.setupFitter(fitter, data_combined_dj);

    const int nsample = data_combined_j.value(0).size();
    readFrozenParams(fitter,args.freeze_file,nsample);
    
    TwoPointThreePointSimFitParams guess;  
    jackknifeDistribution<TwoPointThreePointSimFitParams> params(nsample, guess);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, data_combined_j);
    
    std::cout << "Params: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

    writeParamsStandard(params,"params_7basis.hdf5");
    
    std::vector<jackknifeDistributionD> params_10basis;
    vectorizeAndConvert10basis(params_10basis, params);

    writeParamsStandard(params_10basis,"params_10basis.hdf5");    
  }

};

void fit(const FitFunc &fitfunc,
	 const jackAmplitudeSimCorrelationFunction &data_combined_j,
	 const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,
	 const Args &args){
  switch(args.correlation_status){
  case Uncorrelated:
    return fitSpec<UncorrelatedFit>::fit(fitfunc,data_combined_j,data_combined_dj,args);
  case Correlated:
    return fitSpec<CorrelatedFit>::fit(fitfunc,data_combined_j,data_combined_dj,args);    
  case PartiallyCorrelated:
    return fitSpec<PartiallyCorrelatedFit>::fit(fitfunc,data_combined_j,data_combined_dj,args);
  default:
    error_exit(std::cout << "Unknown correlationStatus " << args.correlation_status << std::endl);    
  }
}


#endif

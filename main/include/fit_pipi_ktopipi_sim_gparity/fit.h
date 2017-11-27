#ifndef _FIT_PIPI_KTOPIPI_SIM_GPARITY_FIT_H_
#define _FIT_PIPI_KTOPIPI_SIM_GPARITY_FIT_H_
#include<plot.h>
#include<gsl_eigensolve.h>
#include<pvalue.h>

typedef FitTwoPointThreePointSim FitFunc;

struct UncorrelatedFit;
struct CorrelatedFit;
struct PartiallyCorrelatedFit;


void sort(std::vector<NumericVector<jackknifeDistributionD> > &evecs_out,
	  std::vector<jackknifeDistributionD> &evals_out,
	  const std::vector<NumericVector<jackknifeDistributionD> > &evecs_in,
	  const std::vector<jackknifeDistributionD> &evals_in){
  const int nev = evals_in.size();
  const int nsample = evals_in[0].size();

  evecs_out.resize(nev, NumericVector<jackknifeDistributionD>(nev, jackknifeDistributionD(nsample)));
  evals_out.resize(nev, jackknifeDistributionD(jackknifeDistributionD(nsample)));

  for(int s=0;s<nsample;s++){  
    std::vector< std::pair<int,double> > evs(nev);
    for(int i=0;i<nev;i++){
      evs[i].first = i;
      evs[i].second = evals_in[i].sample(s);
    }
    std::sort(evs.begin(), evs.end(),
	      [](const std::pair<int,double>& a, const std::pair<int,double>& b) {
		return a.second > b.second;
	      }
	      );
    for(int i=0;i<nev;i++){
      const int mpi = evs[i].first;
      evals_out[i].sample(s) = evals_in[mpi].sample(s);
      for(int j=0;j<nev;j++) evecs_out[i][j].sample(s) =  evecs_in[mpi][j].sample(s);
    }
  }
}

void plotSampleEvecs(const std::vector<NumericVector<jackknifeDistributionD> > &evecs, const std::string &file_stub, const int sample = 0){
  MatPlotLibScriptGenerate plot;
  
  class accessor{
    const int sample;
    const int idx;
    const std::vector<NumericVector<jackknifeDistributionD> > &evecs;
  public:
    accessor(const std::vector<NumericVector<jackknifeDistributionD> > &_evecs, const int _idx, const int _sample): idx(_idx), sample(_sample), evecs(_evecs){
    }
    
    inline double x(const int i) const{ return i; }
    inline double y(const int i) const{ return evecs[idx](i).sample(sample); }
    inline double dxm(const int i) const{ return 0; }
    inline double dxp(const int i) const{ return 0; }
    inline double dym(const int i) const{ return 0; }
    inline double dyp(const int i) const{ return 0; }
    
    inline int size() const{ return evecs.size(); }
  };

  for(int i=0;i<evecs.size();i++){
    accessor acc(evecs, i, sample);
    plot.plotData(acc);
  }

  plot.write(file_stub + ".py", file_stub + ".eps");
}

void calcEvecsEvals(std::vector<NumericVector<jackknifeDistributionD> > &evecs,
		    std::vector<jackknifeDistributionD> &evals,
		    const NumericSquareMatrix<jackknifeDistributionD> &mat,
		    const std::string &descr,
		    const bool do_scale){
  
  const int ndata = mat.size();
  const int nev = ndata;
  NumericSquareMatrix<jackknifeDistributionD> mattmp;
  NumericSquareMatrix<jackknifeDistributionD> const* matp = &mat;
  jackknifeDistributionD scale;
  if(do_scale){
    //Matrix can have very large entries. Convenient to scale such that it's easier to judge if it converged; here with the first diagonal element (ie sigma(0)^2)
    scale = mat(0,0);
    mattmp = NumericSquareMatrix<jackknifeDistributionD>(ndata, [&](const int i, const int j){ return mat(i,j)/scale; });
    matp = &mattmp;
  }

  std::vector<jackknifeDistributionD> residuals = symmetricMatrixEigensolve(evecs,evals,*matp,true);

  if(do_scale) for(int i=0;i<nev;i++) evals[i] = evals[i] * scale; //rescale the evals
  
  std::cout << "Computed eigenvalues of " << descr << " of size " << mat.size() << ":\n";
  for(int i=0;i<nev;i++) std::cout << i << " " << evals[i] << std::endl;
  
  bool fail = false;
  std::cout << "Residuals:\n";
  for(int i=0;i<residuals.size();i++){
    std::cout << i << " " << residuals[i] << std::endl;
    double min,max;
    residuals[i].range(min,max);
    if(max > 1e-10) fail = true;
  }
  if(fail) error_exit(std::cout << "Failed to compute eigenvectors/values of " << descr << ":\n" << mat << std::endl);
}


void analyzeCorrelationMatrix(const NumericSquareMatrix<jackknifeDistributionD> &corr, const NumericSquareMatrix<jackknifeDistributionD> &inv_corr, const NumericVector<jackknifeDistributionD> &sigma,
			      const jackAmplitudeSimCorrelationFunction &data, const jackknifeDistribution<TwoPointThreePointSimFitParams> &params, const FitFunc &fitfunc){
  //Examine the eigenvalues and eigenvectors of the correlation matrix
  const int nsample = sigma(0).size();
  const int nev = corr.size();

  std::vector<NumericVector<jackknifeDistributionD> > evecs;
  std::vector<jackknifeDistributionD> evals;
  calcEvecsEvals(evecs, evals, corr, "correlation matrix", false);

  //Compute the contributions of each evec towards chi^2
  const int ndata = data.size();
  NumericVector<jackknifeDistributionD> Delta(ndata, jackknifeDistributionD(nsample));
  for(int i=0;i<ndata;i++)
    for(int s=0;s<nsample;s++)
      Delta[i].sample(s) = (data.value(i).sample(s) - fitfunc.value(data.coord(i), params.sample(s)))/sigma(i).sample(s);
  
  std::vector<jackknifeDistributionD> Delta_dot_v(nev);
  for(int i=0;i<nev;i++)
    Delta_dot_v[i] = dot(Delta, evecs[i]);

  jackknifeDistributionD chisq(nsample); zeroit(chisq);
  std::vector<jackknifeDistributionD> chisq_contrib(nev);
  std::cout << "chi^2 contributions per eigenvalue:\n";
  for(int i=0;i<nev;i++){
    chisq_contrib[i] = Delta_dot_v[i] * Delta_dot_v[i] / evals[i];
    std::cout << i << " eval = " << evals[i] << " d(chi^2) = " << chisq_contrib[i] << std::endl;
    chisq = chisq + chisq_contrib[i];
  }
  std::cout << "Total chi^2 = " << chisq << std::endl;

  //Write the correlation matrix and its inverse to disk
#ifdef HAVE_HDF5
  writeParamsStandard(chisq_contrib, "corrmat_evecs_chisq_contrib.hdf5");
  writeParamsStandard(evals,"corrmat_evals.hdf5");
  {
    std::vector<std::vector<jackknifeDistributionD> > tmp(nev, std::vector<jackknifeDistributionD>(nev));
    for(int i=0;i<nev;i++) for(int j=0;j<nev;j++) tmp[i][j] = evecs[i](j);
    writeParamsStandard(tmp,"corrmat_evecs.hdf5");
  }
  {
    HDF5writer wr("corrmat.hdf5");
    write(wr,corr,"correlation_matrix");
    write(wr,inv_corr,"inverse_correlation_matrix");    
  }  
#endif

  plotSampleEvecs(evecs, "corrmat_evecs", 0);

  //Do the same for the covariance matrix  
  NumericSquareMatrix<jackknifeDistributionD> cov(ndata, [&](const int i, const int j){ return corr(i,j) * sigma(i) * sigma(j); });

  calcEvecsEvals(evecs, evals, cov, "covariance matrix", true);
  
  NumericVector<jackknifeDistributionD> Delta_abs(ndata, [&](const int i){ return Delta(i) * sigma(i); });
  
  for(int i=0;i<nev;i++)
    Delta_dot_v[i] = dot(Delta_abs, evecs[i]);

  zeroit(chisq);
  std::cout << "chi^2 contributions per eigenvalue:\n";
  for(int i=0;i<nev;i++){
    chisq_contrib[i] = Delta_dot_v[i] * Delta_dot_v[i] / evals[i];
    std::cout << i << " eval = " << evals[i] << " d(chi^2) = " << chisq_contrib[i] << std::endl;
    chisq = chisq + chisq_contrib[i];
  }
  std::cout << "Total chi^2 = " << chisq << std::endl;
  
  //Write the covariance matrix to disk
#ifdef HAVE_HDF5
  writeParamsStandard(chisq_contrib,"covmat_evecs_chisq_contrib.hdf5");
  writeParamsStandard(evals,"covmat_evals.hdf5");
  {
    std::vector<std::vector<jackknifeDistributionD> > tmp(nev, std::vector<jackknifeDistributionD>(nev));
    for(int i=0;i<nev;i++) for(int j=0;j<nev;j++) tmp[i][j] = evecs[i](j);
    writeParamsStandard(tmp,"covmat_evecs.hdf5");
  }
  {
    HDF5writer wr("covmat.hdf5");
    write(wr,cov,"covariance_matrix");
  }  
#endif

  plotSampleEvecs(evecs, "covmat_evecs", 0);
}

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
  void postFitAnalysis(const jackAmplitudeSimCorrelationFunction &data_combined_j, const jackknifeDistribution<TwoPointThreePointSimFitParams> &params, const FitFunc &fitfunc){
  }
};
template<>
struct CorrelationPolicy<CorrelatedFit>{
  typedef typename composeFitPolicy<amplitudeDataCoordSim, FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  std::unique_ptr<importCostFunctionParameters<correlatedFitPolicy,FitPolicies> > prep;
    
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    prep.reset(new importCostFunctionParameters<correlatedFitPolicy,FitPolicies>(fitter, data_combined_dj));    
  }
  void postFitAnalysis(const jackAmplitudeSimCorrelationFunction &data_combined_j, const jackknifeDistribution<TwoPointThreePointSimFitParams> &params, const FitFunc &fitfunc){
    analyzeCorrelationMatrix(prep->corr, prep->inv_corr, prep->sigma, data_combined_j, params, fitfunc);
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
  }

  void postFitAnalysis(const jackAmplitudeSimCorrelationFunction &data_combined_j, const jackknifeDistribution<TwoPointThreePointSimFitParams> &params, const FitFunc &fitfunc){
    analyzeCorrelationMatrix(corr, inv_corr, sigma, data_combined_j, params, fitfunc);
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

    double dof = chisq.sample(0)/chisq_per_dof.sample(0);
    
    std::cout << "Params: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
    std::cout << "Dof: " << dof << std::endl;
    
    writeParamsStandard(chisq, "chisq.hdf5");
    writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");

    //Compute p-value
    jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });
    std::cout << "p-value: " << pvalue << std::endl;
    writeParamsStandard(pvalue, "pvalue.hdf5");

    //Other post-fit analysis
    prepare.postFitAnalysis(data_combined_j, params, fitfunc);

    //Write params
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

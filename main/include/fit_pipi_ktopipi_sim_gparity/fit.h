#ifndef _FIT_PIPI_KTOPIPI_SIM_GPARITY_FIT_H_
#define _FIT_PIPI_KTOPIPI_SIM_GPARITY_FIT_H_
#include<plot.h>
#include<gsl_eigensolve.h>
#include<pvalue.h>

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

void plotRotatedFitFunc(const jackAmplitudeSimCorrelationFunction &data,
			const std::vector<NumericVector<jackknifeDistributionD> > &evecs, const std::vector<jackknifeDistributionD> &evals,
			const std::vector<jackknifeDistributionD> &fit_vals){
  //For *covariance matrix* write
  //\Delta'_i = (1/\sqrt(\lambda_i)) (\vec v_i^T \Delta) = ( y'_i - f'_i )
  //f'_i = (1/\sqrt(\lambda_i)) \sum_j v^i_j f(x_j)

  //Examine  (1/\sqrt(\lambda_i))v^i_j  ,  f(x_j) and their product as a function of j (time)
  const int ndata = data.size();
  const int nev = evecs.size();
  jackknifeCorrelationFunctionD vnorm(ndata);
  jackknifeCorrelationFunctionD f(ndata);
  jackknifeCorrelationFunctionD prod(ndata);

  MatPlotLibScriptGenerate plot;
  typedef DataSeriesAccessor<jackknifeCorrelationFunctionD, ScalarCoordinateAccessor<double>, DistributionPlotAccessor<jackknifeDistributionD> > accessor;
  
  for(int i=0;i<nev;i++){

    const jackknifeDistributionD &nrm = fit_vals[0];
    
    for(int j=0;j<ndata;j++){
      vnorm.coord(j) = f.coord(j) = prod.coord(j) = data.coord(j).t;
      vnorm.value(j) = evecs[i](j)/sqrt(evals[i])*nrm;
      f.value(j) = fit_vals[j]/nrm;
      prod.value(j) = vnorm.value(j) * f.value(j);
    }
    typename MatPlotLibScriptGenerate::kwargsType args;
    
    accessor acc0(vnorm);
    args["color"] = "r";
    args["marker"] = "o";
    plot.plotData(acc0,args);

    accessor acc1(f);
    args["color"] = "g";
    args["marker"] = "s";
    plot.plotData(acc1,args);
    
    accessor acc2(prod);
    args["color"] = "b";
    args["marker"] = "^";
    plot.plotData(acc2,args);
  }

  plot.invoke() << "\tdatasets = [";
  for(int i=0;i<3*nev-1;i++) plot.invoke() << "dset" << i << ", ";
  plot.invoke() << "dset" << 3*nev-1 << "]\n";
  
  plot.write("covmat_rotated_fitfunc.py","covmat_rotated_fitfunc.eps");  
}

template<typename FitFunc>
void analyzeCorrelationMatrix(const NumericSquareMatrix<jackknifeDistributionD> &corr, const NumericSquareMatrix<jackknifeDistributionD> &inv_corr, const NumericVector<jackknifeDistributionD> &sigma,
			      const jackAmplitudeSimCorrelationFunction &data, const jackknifeDistribution<typename FitFunc::Params> &params, const FitFunc &fitfunc){
  //Examine the eigenvalues and eigenvectors of the correlation matrix
  const int nsample = sigma(0).size();
  const int nev = corr.size();

  std::vector<NumericVector<jackknifeDistributionD> > evecs;
  std::vector<jackknifeDistributionD> evals;
  calcEvecsEvals(evecs, evals, corr, "correlation matrix", false);

  //Compute the contributions of each evec towards chi^2

  // \chi^2 = \vec \Delta^T M^{-1} \vec \Delta
  // \Delta_i = (y_i - f(x_i))/\sigma_i

  // M is the correlation matrix.
  // After computing evals \lambda , evecs  \vec v

  // M^{-1} = \sum_i \vec v_i (1/\lambda_i) \vec v_i^T

  // \delta\chi^2_i = (\vec \Delta^T \vec v_i) (1/\lambda_i) (\vec v_i^T \Delta)

  //Note if we define  \Delta'_i = (1/\sqrt(\lambda_i)) (\vec v_i^T \Delta) then
  // \chi^2 = \sum_{ij} \Delta'_i \delta_{ij} \Delta'_j
  //        = \vec \Delta'^T I \vec \Delta'     where I is the unit matrix.
  //Thus \Delta is rotated into an uncorrelated basis

  const int ndata = data.size();
  NumericVector<jackknifeDistributionD> Delta(ndata, jackknifeDistributionD(nsample));
  std::vector<jackknifeDistributionD> fit_vals(ndata, jackknifeDistributionD(nsample));
  for(int i=0;i<ndata;i++)
    for(int s=0;s<nsample;s++){
      fit_vals[i].sample(s) = fitfunc.value(data.coord(i), params.sample(s));
      Delta[i].sample(s) = ( data.value(i).sample(s) - fit_vals[i].sample(s) )/sigma(i).sample(s);
    }
      
  std::vector<jackknifeDistributionD> Delta_dot_v(nev);
  std::vector<jackknifeDistributionD> Delta_prime(nev);
  for(int i=0;i<nev;i++){
    Delta_dot_v[i] = dot(Delta, evecs[i]);
    Delta_prime[i] = Delta_dot_v[i]/sqrt(evals[i]);
  }
    
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
  writeParamsStandard(Delta_prime, "corrmat_delta_prime.hdf5");
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
  
  
  for(int i=0;i<nev;i++){
    Delta_dot_v[i] = dot(Delta_abs, evecs[i]);
    Delta_prime[i] = Delta_dot_v[i]/sqrt(evals[i]);
  }
  
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
  writeParamsStandard(Delta_prime, "covmat_delta_prime.hdf5");
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

  plotRotatedFitFunc(data,evecs,evals,fit_vals);

  //For *covariance matrix* write
  //\Delta'_i = (1/\sqrt(\lambda_i)) (\vec v_i^T \Delta) = ( y'_i - f'_i )
  //f'_i = (1/\sqrt(\lambda_i)) \sum_j v^i_j f(x_j)
  //Compute derivatives of f'_i with respect to parameters at minima
  const int nparams = fitfunc.Nparams();
  std::vector<std::vector<jackknifeDistributionD > > all_deriv(ndata, std::vector<jackknifeDistributionD >(nparams, jackknifeDistributionD(nsample))); //[data idx][param]
  for(int tt=0;tt<ndata;tt++){
    for(int s=0;s<nsample;s++){
      typename FitFunc::ValueDerivativeType d = fitfunc.parameterDerivatives(data.coord(tt), params.sample(s));
      for(int p=0;p<nparams;p++) all_deriv[tt][p].sample(s) = d(p);
    }
  }

  for(int i=0;i<nev;i++){
    //f'_i = (1/\sqrt(\lambda_i)) \sum_j v^i_j f(x_j)
    jackknifeDistributionD f(nsample,0.);
    for(int j=0;j<ndata;j++)
	f = f + evecs[i](j) * fit_vals[j];
    f = f / sqrt(evals[i]);
    
    //d/dp(f'_i) = (1/\sqrt(\lambda_i)) \sum_j v^i_j df(x_j)/dp
    std::cout << "Evec idx " << i << "rotated function " << f << " , derivative of rotated function wrt params and sensitivity ( df/dp * p/f ):\n";
    for(int p=0;p<nparams;p++){
      jackknifeDistributionD df_by_dp(nsample,0.);
      for(int j=0;j<ndata;j++)
	df_by_dp = df_by_dp + evecs[i](j) * all_deriv[j][p];
      df_by_dp = df_by_dp / sqrt(evals[i]);

      jackknifeDistributionD pv(nsample, [&](const int s){ return params.sample(s)(p); });
      jackknifeDistributionD sens = df_by_dp * pv / f;
      
      std::cout << p << " " << df_by_dp << " " << sens << std::endl;
    }
  }
  
}

template<typename CorrelationStatusStruct, typename FitFunc>
struct CorrelationPolicy{};

template<typename FitFunc>
struct CorrelationPolicy<UncorrelatedFit, FitFunc>{
  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, uncorrelatedFitPolicy>::type FitPolicies;
  std::unique_ptr<importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies> > prep;
    
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    prep.reset(new importCostFunctionParameters<uncorrelatedFitPolicy,FitPolicies>(fitter, data_combined_dj));
  }
  void postFitAnalysis(const jackAmplitudeSimCorrelationFunction &data_combined_j, const jackknifeDistribution<typename FitFunc::Params> &params, const FitFunc &fitfunc){
  }
};
template<typename FitFunc>
struct CorrelationPolicy<CorrelatedFit, FitFunc>{
  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  std::unique_ptr<importCostFunctionParameters<correlatedFitPolicy,FitPolicies> > prep;
    
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    prep.reset(new importCostFunctionParameters<correlatedFitPolicy,FitPolicies>(fitter, data_combined_dj));    
  }
  void postFitAnalysis(const jackAmplitudeSimCorrelationFunction &data_combined_j, const jackknifeDistribution<typename FitFunc::Params> &params, const FitFunc &fitfunc){
    analyzeCorrelationMatrix(prep->corr, prep->inv_corr, prep->sigma, data_combined_j, params, fitfunc);
  }
};
template<typename FitFunc>
struct CorrelationPolicy<PartiallyCorrelatedFit, FitFunc>{
  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;
  typedef typename FitPolicies::MatrixDistribution MatrixDistribution;
  typedef typename FitPolicies::VectorDistribution VectorDistribution;

  MatrixDistribution corr;
  MatrixDistribution inv_corr;
  VectorDistribution sigma;
  
  void setupFitter(fitter<FitPolicies> &fitter,
		   const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj){
    //Only compute the correlations for data of the same type idx (pipi, Q1', Q2', Q3', ....)
    const int ndata = data_combined_dj.size();
    const int nsample = data_combined_dj.value(0).size();    

    sigma.resize(ndata);   
    
    MatrixDistribution cov(ndata, jackknifeDistributionD(nsample,0.));
    for(int i=0;i<ndata;i++){
      cov(i,i) = doubleJackknifeDistributionD::covariance(data_combined_dj.value(i),data_combined_dj.value(i));
      sigma(i) = sqrt(cov(i,i));
      
      for(int j=i+1;j<ndata;j++){
	if(data_combined_dj.coord(i).idx == data_combined_dj.coord(j).idx)
	  cov(i,j) = cov(j,i) = doubleJackknifeDistributionD::covariance(data_combined_dj.value(i),data_combined_dj.value(j));
      }
    }

    corr = MatrixDistribution(ndata);
    
    for(int i=0;i<ndata;i++){
      corr(i,i) = jackknifeDistributionD(nsample,1.);
      
      for(int j=i+1;j<ndata;j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma(i)/sigma(j);
    }
    
    inv_corr.resize(ndata, jackknifeDistributionD(nsample));
    svd_inverse(inv_corr, corr);

    //Test the quality of the inverse
    MatrixDistribution test = corr * inv_corr;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistributionD(nsample,1.0);    
    std::cout << "|CorrMat * CorrMat^{-1} - 1|^2 = " << mod2(test) << std::endl;

    //Import
    fitter.importCostFunctionParameters(inv_corr,sigma);
  }

  void postFitAnalysis(const jackAmplitudeSimCorrelationFunction &data_combined_j, const jackknifeDistribution<typename FitFunc::Params> &params, const FitFunc &fitfunc){
    analyzeCorrelationMatrix(corr, inv_corr, sigma, data_combined_j, params, fitfunc);
  }
  
};



template<typename CorrelationStatusStruct, typename FitFunc>
struct fitSpec{
  static void fit(const FitFunc &fitfunc,
		  const jackAmplitudeSimCorrelationFunction &data_combined_j,
		  const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,
		  const Args &args){
    typedef typename CorrelationPolicy<CorrelationStatusStruct,FitFunc>::FitPolicies FitPolicies;
    CorrelationPolicy<CorrelationStatusStruct,FitFunc> prepare;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);
    prepare.setupFitter(fitter, data_combined_dj);

    const int nsample = data_combined_j.value(0).size();
    readFrozenParams(fitter,args.freeze_file,nsample);
    
    typename FitFunc::Params guess;  
    jackknifeDistribution<typename FitFunc::Params> params(nsample, guess);
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

template<typename FitFunc>
void fitFspec(const FitFunc &fitfunc,
	 const jackAmplitudeSimCorrelationFunction &data_combined_j,
	 const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,
	 const Args &args){
  switch(args.correlation_status){
  case Uncorrelated:
    return fitSpec<UncorrelatedFit,FitFunc>::fit(fitfunc,data_combined_j,data_combined_dj,args);
  case Correlated:
    return fitSpec<CorrelatedFit,FitFunc>::fit(fitfunc,data_combined_j,data_combined_dj,args);    
  case PartiallyCorrelated:
    return fitSpec<PartiallyCorrelatedFit,FitFunc>::fit(fitfunc,data_combined_j,data_combined_dj,args);
  default:
    error_exit(std::cout << "Unknown correlationStatus " << args.correlation_status << std::endl);    
  }
}

void fit(const jackAmplitudeSimCorrelationFunction &data_combined_j,
	 const doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,
	 const Args &args){
  if(args.fitfunc == OneExp){
    FitTwoPointThreePointSim fitfunc(args.Lt,args.tsep_pipi,args.Ascale,args.Cscale);
    fitFspec(fitfunc,data_combined_j,data_combined_dj,args);
  }else if(args.fitfunc == TwoExpPiPi){
    FitTwoExpTwoPointThreePointSim fitfunc(args.Lt,args.tsep_pipi,args.Ascale,args.Cscale);
    fitFspec(fitfunc,data_combined_j,data_combined_dj,args);
  }else{
    error_exit(std::cout << "Unknown fit function " << args.fitfunc << std::endl);
  }
}


#endif

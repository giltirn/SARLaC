#ifndef _CPSFIT_SIMPLE_FIT_WRAPPER_FITTER_SINGLE_H
#define _CPSFIT_SIMPLE_FIT_WRAPPER_FITTER_SINGLE_H

#include<config.h>
#include<utils/macros.h>

#include<fit/simple_fit_wrapper/fit_common.h>
#include<distribution/jackknife.h>

CPSFIT_START_NAMESPACE

//Fit a single set of means (not a distribution). You can of course use the version designed for jackknife distributions with 1 sample, but this version avoids all the overheads
class simpleSingleFitWrapper{
  INHERIT_COMMON_TYPEDEFS;

  const FitFunc &fitfunc;
  MinimizerType min_type;
  generalContainer min_params;

  NumericSquareMatrix<double> corr_mat;
  std::vector<double> sigma;
  bool have_corr_mat;
  bool corr_mat_preinverted; //is the correlation matrix provided already inverted?

  std::vector<double> freeze_values;
  std::vector<int> freeze_params;

  std::vector<boundedParameterTransform> bounded_trans;

  std::vector<Prior> priors;

  NumericSquareMatrix<double> invertCorrelationMatrix() const{
    NumericSquareMatrix<double> inv_corr_mat(corr_mat);
    double condition_number;
    svd_inverse(inv_corr_mat, corr_mat, condition_number);

    //Test the quality of the inverse
    NumericSquareMatrix<double> test = corr_mat * inv_corr_mat;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - 1.;
    double resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid << std::endl;
    return inv_corr_mat;
  }

public:  
  simpleSingleFitWrapper(const FitFunc &fitfunc, 
			 const MinimizerType min_type, 
			 const generalContainer &min_params = generalContainer()): fitfunc(fitfunc), min_type(min_type), min_params(min_params), have_corr_mat(false), corr_mat_preinverted(false){
  }
  
  //Can be changed at any time
  inline void setMinimizer(const MinimizerType _min_type, const generalContainer &_min_params= generalContainer()){
    min_type = _min_type; min_params = _min_params;
  }

  //Fix multiple parameters to input values
  inline void freeze(const std::vector<int> &_freeze_params,
	      const std::vector<double> &_freeze_values){
    freeze_params = _freeze_params;
    freeze_values = _freeze_values;
  } 
  //Fix single parameter to input value
  inline void freeze(const int idx, const double &val){
    freeze_params.push_back(idx);
    freeze_values.push_back(val);
  }
  //Add a Gaussian prior
  inline void addPrior(const double value, const double weight, const int param_idx){
    priors.push_back(Prior(value,weight,param_idx));
  }
  
  //Add bounds on parameters
  inline void setBound(int param, const ParameterBound bound, double min, double max){ //min or max will be ignored as appropriate if using Max or Min bounds
    bounded_trans.push_back( boundedParameterTransform(param, bound, min, max) );
  }
  inline void setBound(const boundedParameterTransform &t){ bounded_trans.push_back(t); }

  //Import a pre-generated covariance matrix
  void importCovarianceMatrix(const NumericSquareMatrix<double> &cov, 
			      const CostType cost_type = CostType::Correlated){
    int ndata = cov.size();
    corr_mat.resize(ndata);
    sigma.resize(ndata);

    for(int i=0;i<ndata;i++) sigma[i] = sqrt(cov(i,i));
    
    for(int i=0;i<ndata;i++){
      corr_mat(i,i) = cov(i,i)/sigma[i]/sigma[i];
      for(int j=i+1;j<ndata;j++)
	corr_mat(i,j) = corr_mat(j,i) =  cost_type == CostType::Correlated ? cov(i,j)/sigma[i]/sigma[j] : 0.;      
    }														

    have_corr_mat = true;
  }
  //Import a pre-generated correlation matrix and weights sigma   (sigma_i = sqrt(cov_ii))
  void importCorrelationMatrix(const NumericSquareMatrix<double> &corr, const std::vector<double> &sigma_in){
    corr_mat = corr;
    sigma = sigma_in;
    have_corr_mat = true;
  } 
  //Import a pre-generated *inverse* correlation matrix and weights sigma   (sigma_i = sqrt(cov_ii))
  void importInverseCorrelationMatrix(const NumericSquareMatrix<double> &inv_corr, const std::vector<double> &sigma_in){
    corr_mat = inv_corr;
    sigma = sigma_in;
    have_corr_mat = true;
    corr_mat_preinverted = true;
  }

  //Generate the covariance matrix internally from jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename T>
  void generateCovarianceMatrix(const correlationFunction<T, jackknifeDistribution<double> > &data_j, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_j.size();
    NumericSquareMatrix<double> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov(i,j) = cov(j,i) = jackknifeDistribution<double>::covariance(data_j.value(i), data_j.value(j));

    importCovarianceMatrix(cov, cost_type);
  }


  //Get the correlation matrix from the unbinned jackknife and sigma from the binned jackknife (the hybrid approach)
  template<typename T>
  void generateCovarianceMatrix(const correlationFunction<T, jackknifeDistribution<double>> &data_j_binned,
				const correlationFunction<T, jackknifeDistribution<double>> &data_j_unbinned, 
				const CostType cost_type = CostType::Correlated){
    assert(data_j_binned.size() == data_j_unbinned.size());
    assert(data_j_binned.value(0).size() <= data_j_unbinned.value(0).size() );

    int ndata = data_j_binned.size();

    NumericSquareMatrix<double> cov_unbinned(ndata);    
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov_unbinned(i,j) = cov_unbinned(j,i) = jackknifeDistribution<double>::covariance(data_j_unbinned.value(i), data_j_unbinned.value(j));
    
    std::vector<double> sigma_unbinned(ndata);
    sigma.resize(ndata);
    for(int i=0;i<ndata;i++){
      sigma[i] = sqrt( jackknifeDistribution<double>::covariance(data_j_binned.value(i), data_j_binned.value(i)) ); //sigma from binned jackknife
      sigma_unbinned[i] = sqrt(cov_unbinned(i,i));
    }
      
    //Corr from unbinned data
    corr_mat.resize(ndata);

    for(int i=0;i<ndata;i++){
      corr_mat(i,i) = cov_unbinned(i,i)/sigma_unbinned[i]/sigma_unbinned[i];
      for(int j=i+1;j<ndata;j++)
	corr_mat(i,j) = corr_mat(j,i) =  cost_type == CostType::Correlated ? cov_unbinned(i,j)/sigma_unbinned[i]/sigma_unbinned[j] : 0.;      
    }
										       
    have_corr_mat = true;
  }

  //Write the covariance matrix to a file in HDF5 format for external manipulation
  void writeCovarianceMatrixHDF5(const std::string &file) const{
    if(corr_mat_preinverted) error_exit(std::cout << "simpleFitWrapper::writeCovarianceMatrixHDF5 function inapplicable if covariance matrix inverse is precomputed\n");
#ifdef HAVE_HDF5
    if(!have_corr_mat) error_exit(std::cout << "simpleFitWrapper::writeCovarianceMatrixHDF5  No covariance/correlation matrix available. Make sure you import one before calling this method!\n");
    NumericSquareMatrix<double> cov = corr_mat;
    for(int i=0;i<cov.size();i++)
      for(int j=0;j<cov.size();j++)
	cov(i,j) = cov(i,j) * sigma[i] * sigma[j];
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
  
  inline const NumericSquareMatrix<double> & getCorrelationMatrix() const{
    assert(have_corr_mat); return corr_mat;
  }
  inline const std::vector<double> & getSigma() const{
    assert(have_corr_mat); return sigma;
  }
  
  //Note the parameter type P is translated internally into a parameterVector  (requires the usual size() and operator()(const int) methods)
  //The coordinate type is wrapped up in a generalContainer as this is only ever needed by the fit function (which knows what type it is and can retrieve it)
  //If chisq_dof_nopriors pointer is provided, the chisq computed without priors and the number of degrees of freedom without priors will be written there
  //Returns true if converged
  template<typename P, typename T>
  bool fit(P &params,
	   double &chisq,
	   double &chisq_per_dof,
	   int &dof,
	   const correlationFunction<T, double> &data,
	   std::pair<double, int>* chisq_dof_nopriors = NULL){
    if(!have_corr_mat) error_exit(std::cout << "simpleSingleFitWrapper::fit  No covariance/correlation matrix available. Make sure you import one before calling this method!\n");
    int nparam = params.size();
    int ndata = data.size();
    if(ndata == 0){
      std::cout << "Warning: Fit data container contains no data! Not performing a fit...." << std::endl;
      dof = -1;
      return false;
    }
    if(getCorrelationMatrix().size() != ndata) error_exit(std::cout << "simpleSingleFitWrapper::fit size of covariance matrix " << getCorrelationMatrix().size() << " does not match data size " << ndata << "\n");

    //Prepare inverse correlation matrix (the condition number is useful information even if we don't need the inverse explicitly)
    NumericSquareMatrix<double> inv_corr_mat = corr_mat_preinverted ? corr_mat : invertCorrelationMatrix();

    //Repack the data
    correlationFunction<generalContainer, double> data_i(ndata);
    for(int i=0;i<ndata;i++){
      data_i.coord(i) = generalContainer(data.coord(i));
      data_i.value(i) = data.value(i);
    }

    //Prepare the inner fit parameter type
    ParameterType params_i(nparam);
    for(int i=0;i<nparam;i++) params_i(i) = params(i);
    
    //Setup the sample fit function (frozen fit func needs a properly setup instance of the parameter type even if not freezing anything)
    FitFuncBounded fitfunc_bi(fitfunc);
    for(int t=0;t<bounded_trans.size();t++) fitfunc_bi.setBound(bounded_trans[t]);
      
    //Map guess parameters to internal params
    for(int t=0;t<bounded_trans.size();t++) bounded_trans[t].mapBoundedToUnbounded(params_i);

    FitFuncFrozenBounded fitfunc_i(fitfunc_bi);
    if(freeze_params.size() > 0){
      ParameterType frzp(params_i);
      for(int p=0;p<freeze_params.size();p++)
	frzp(freeze_params[p]) = freeze_values[p];
      fitfunc_i.freeze(freeze_params, frzp);
      params_i = fitfunc_i.mapParamsSupersetToSubset(params_i); //pull out the subset of parameters that are varied; here they are the same type
    }else{
      fitfunc_i.freeze({}, params_i);
    }

    //Run the fitter
    bool converged;
    switch(min_type){
    case MinimizerType::MarquardtLevenberg:
      chisq = simpleFitCommon::fitSampleML(converged, params_i, dof, data_i, inv_corr_mat, sigma, fitfunc_i, min_params, priors, chisq_dof_nopriors); break;
    case MinimizerType::GSLtrs:
      if(corr_mat_preinverted) error_exit(std::cout << "Cannot use GSLtrs with precomputed inverse covariance matrix" << std::endl);
      chisq = simpleFitCommon::fitSampleGSLtrs(converged, params_i, dof, data_i, corr_mat, sigma, fitfunc_i, min_params, priors, chisq_dof_nopriors); break;
    case MinimizerType::GSLmultimin:
      chisq = simpleFitCommon::fitSampleGSLmultimin(converged, params_i, dof, data_i, inv_corr_mat, sigma, fitfunc_i, min_params, priors, chisq_dof_nopriors); break;
    case MinimizerType::Minuit2:
      chisq = simpleFitCommon::fitSampleMinuit2(converged, params_i, dof, data_i, inv_corr_mat, sigma, fitfunc_i, min_params, priors, chisq_dof_nopriors); break;
    default:
      assert(0);
    }
      
    //If you want the fitter to fail on non-convergence, set the associated minimizer parameter which typically defaults to do so
    if(!converged) std::cout << "Warning: thread " << omp_get_thread_num() << " fit did not converge" << std::endl;

    //Re-enlarge the parameter vector with the frozen fit parameters
    if(freeze_params.size() > 0) params_i = fitfunc_i.mapParamsSubsetToSuperset(params_i);

    //Map internal parameters to external params
    for(int t=0;t<bounded_trans.size();t++) bounded_trans[t].mapUnboundedToBounded(params_i);
      
    chisq_per_dof = chisq/dof;
    for(int i=0;i<nparam;i++) params(i) = params_i(i);

    return converged;
  }  
};

CPSFIT_END_NAMESPACE

#endif

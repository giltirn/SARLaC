#ifndef _CPSFIT_SIMPLE_FIT_WRAPPER_FITTER_H
#define _CPSFIT_SIMPLE_FIT_WRAPPER_FITTER_H

#include<config.h>
#include<utils/macros.h>

#include<fit/simple_fit_wrapper/fit_common.h>
#include<distribution/jackknife.h>
#include<distribution/double_jackknife.h>
#include<distribution/block_double_jackknife.h>
#include<distribution/bootstrap.h>
#include<distribution/boot_jackknife.h>

CPSFIT_START_NAMESPACE

//BaseDistributionType is the distribution type under which the data, chisq, etc are vectorized. Usually this is jackknifeDistribution<double> 
//but it could be bootstrap for example
//The distribution type associated with the fit parameters will have a different underlying data type (eg parameterVector<double>)

template<typename BaseDistributionType>
class simpleFitWrapper{
  INHERIT_COMMON_TYPEDEFS;

  typedef typename BaseDistributionType::DataType BaseNumericType;

  const FitFunc &fitfunc;
  MinimizerType min_type;
  generalContainer min_params;

  NumericSquareMatrix<BaseDistributionType> corr_mat;
  std::vector<BaseDistributionType> sigma;
  bool have_corr_mat;

  std::vector<BaseDistributionType> freeze_values;
  std::vector<int> freeze_params;

  std::vector<boundedParameterTransform> bounded_trans;

  std::vector<Prior> priors;

  NumericSquareMatrix<BaseDistributionType> invertCorrelationMatrix() const{
    int nsample = corr_mat(0,0).size();
    NumericSquareMatrix<BaseDistributionType> inv_corr_mat(corr_mat);
    BaseDistributionType condition_number;
    svd_inverse(inv_corr_mat, corr_mat, condition_number);

    BaseDistributionType one(corr_mat(0,0)); 
    for(int i=0;i<iterate<BaseDistributionType>::size(one); i++) iterate<BaseDistributionType>::at(i, one) = 1.;

    //Test the quality of the inverse
    NumericSquareMatrix<BaseDistributionType> test = corr_mat * inv_corr_mat;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - one;
    BaseDistributionType resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number.best() << " +- " << condition_number.standardDeviation() << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid.best() << " +- " << resid.standardDeviation() << std::endl;
    return inv_corr_mat;
  }

  static inline std::vector<BaseNumericType> sample(const std::vector<BaseDistributionType> &v, const int s){
    std::vector<BaseNumericType> out(v.size()); for(int i=0;i<v.size();i++) out[i] = iterate<BaseDistributionType>::at(s, v[i]);
    return out;
  }
  static inline NumericSquareMatrix<BaseNumericType> sample(const NumericSquareMatrix<BaseDistributionType> &v, const int s){
    NumericSquareMatrix<BaseNumericType> out(v.size()); 
    for(int i=0;i<v.size();i++) 
      for(int j=0;j<v.size();j++) 
	out(i,j) = iterate<BaseDistributionType>::at(s, v(i,j)); 
    return out;
  }

  //A streambuf that filters stream output for a single thread
  class thrbuf: public std::streambuf{
  public:
    int thread;
    std::streambuf *base;

    thrbuf(const int thread, std::streambuf *base): thread(thread), base(base){}

  protected:
    std::streamsize xsputn (const char* s, std::streamsize n){
      return omp_get_thread_num() == thread ? base->sputn(s,n) : n;
    }
    int overflow(int c = std::char_traits<char>::eof()){
      return omp_get_thread_num() == thread ? base->sputc(c) : c;
    }
  };

public:  
  simpleFitWrapper(const FitFunc &fitfunc, 
		   const MinimizerType min_type, 
		   const generalContainer &min_params = generalContainer()): fitfunc(fitfunc), min_type(min_type), min_params(min_params), have_corr_mat(false){
  }
  
  //Can be changed at any time
  inline void setMinimizer(const MinimizerType _min_type, const generalContainer &_min_params= generalContainer()){
    min_type = _min_type; min_params = _min_params;
  }

  //Fix multiple parameters to input values
  inline void freeze(const std::vector<int> &_freeze_params,
	      const std::vector<BaseDistributionType> &_freeze_values){
    freeze_params = _freeze_params;
    freeze_values = _freeze_values;
  } 
  //Fix single parameter to input value
  inline void freeze(const int idx, const BaseDistributionType &val){
    freeze_params.push_back(idx);
    freeze_values.push_back(val);
  }

  //Remove all frozen parameters
  inline void resetFrozenParameters(){
    freeze_params.resize(0);
    freeze_values.resize(0);
  }

  //Add a Gaussian prior
  inline void addPrior(const double value, const double weight, const int param_idx){
    priors.push_back(Prior(value,weight,param_idx));
  }

  //Remove all priors
  inline void resetPriors(){
    priors.resize(0);
  }
  
  //Add bounds on parameters
  inline void setBound(int param, const ParameterBound bound, double min, double max){ //min or max will be ignored as appropriate if using Max or Min bounds
    bounded_trans.push_back( boundedParameterTransform(param, bound, min, max) );
  }
  inline void setBound(const boundedParameterTransform &t){ bounded_trans.push_back(t); }

  //Remove all bounds
  inline void resetBounds(){
    bounded_trans.resize(0);
  }

  //Import a pre-generated covariance matrix
  void importCovarianceMatrix(const NumericSquareMatrix<BaseDistributionType> &cov, 
			      const CostType cost_type = CostType::Correlated){
    int ndata = cov.size();
    corr_mat.resize(ndata);
    sigma.resize(ndata);

    BaseDistributionType zero(cov(0,0)); zeroit(zero);

    for(int i=0;i<ndata;i++) sigma[i] = sqrt(cov(i,i));
    
    for(int i=0;i<ndata;i++){
      corr_mat(i,i) = cov(i,i)/sigma[i]/sigma[i];
      for(int j=i+1;j<ndata;j++)
	corr_mat(i,j) = corr_mat(j,i) =  cost_type == CostType::Correlated ? cov(i,j)/sigma[i]/sigma[j] : zero;      
    }														

    have_corr_mat = true;
  }
  //Import a pre-generated correlation matrix and weights sigma   (sigma_i = sqrt(cov_ii))
  void importCorrelationMatrix(const NumericSquareMatrix<BaseDistributionType> &corr, const std::vector<BaseDistributionType> &sigma_in){
    corr_mat = corr;
    sigma = sigma_in;
    have_corr_mat = true;
  } 
  
#define JACKKNIFE_ONLY typename D=BaseDistributionType, typename std::enable_if<is_jackknife<D>::value, int>::type = 0

  //Generate the covariance matrix internally from double-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, template<typename> class V, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, doubleJackknifeDistribution<BaseNumericType,V>> &data_dj, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_dj.size();
    NumericSquareMatrix<BaseDistributionType> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov(i,j) = cov(j,i) = doubleJackknifeDistribution<BaseNumericType,V>::covariance(data_dj.value(i), data_dj.value(j));

    importCovarianceMatrix(cov, cost_type);
  }
  //Generate the covariance matrix internally from block double-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, template<typename> class V, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, blockDoubleJackknifeDistribution<BaseNumericType,V>> &data_bdj, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_bdj.size();
    NumericSquareMatrix<BaseDistributionType> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov(i,j) = cov(j,i) = blockDoubleJackknifeDistribution<BaseNumericType,V>::covariance(data_bdj.value(i), data_bdj.value(j));

    importCovarianceMatrix(cov, cost_type);
  }

  //Get the correlation matrix from the block double-jack and sigma from the regular, binned double-jackknife (the hybrid approach)
  template<typename GeneralizedCoordinate, template<typename> class V, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, doubleJackknifeDistribution<BaseNumericType,V>> &data_dj,
				const correlationFunction<GeneralizedCoordinate, blockDoubleJackknifeDistribution<BaseNumericType,V>> &data_bdj, 
				const CostType cost_type = CostType::Correlated){
    assert(data_dj.size() == data_bdj.size());
    assert(data_dj.value(0).size() == data_bdj.value(0).size()); //only inner indexing differs

    int ndata = data_dj.size();
    int nsample = data_bdj.value(0).size();

    NumericSquareMatrix<BaseDistributionType> cov_bdj(ndata);    
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov_bdj(i,j) = cov_bdj(j,i) = blockDoubleJackknifeDistribution<BaseNumericType,V>::covariance(data_bdj.value(i), data_bdj.value(j));
    
    std::vector<BaseDistributionType> sigma_bdj(ndata);
    sigma.resize(ndata);
    for(int i=0;i<ndata;i++){
      sigma[i] = sqrt( doubleJackknifeDistribution<BaseNumericType,V>::covariance(data_dj.value(i), data_dj.value(i)) ); //sigma from regular dj
      sigma_bdj[i] = sqrt(cov_bdj(i,i));
    }
      
    //Corr from bdj
    corr_mat.resize(ndata);
    BaseDistributionType zero(nsample); zeroit(zero);

    for(int i=0;i<ndata;i++){
      corr_mat(i,i) = cov_bdj(i,i)/sigma_bdj[i]/sigma_bdj[i];
      for(int j=i+1;j<ndata;j++)
	corr_mat(i,j) = corr_mat(j,i) =  cost_type == CostType::Correlated ? cov_bdj(i,j)/sigma_bdj[i]/sigma_bdj[j] : zero;      
    }
										       
    have_corr_mat = true;
  }


  //Generate the covariance matrix internally from single-jackknife data. The resulting covariance matrix is "frozen", i.e. the same for all samples
  //Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data_j, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_j.size();
    int nsample = data_j.value(0).size();
    NumericSquareMatrix<BaseDistributionType> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++){
	double cv = BaseDistributionType::covariance(data_j.value(i), data_j.value(j));
	cov(i,j) = cov(j,i) = BaseDistributionType(nsample, cv);
      }

    importCovarianceMatrix(cov, cost_type);
  }


#define BOOTSTRAP_ONLY typename D=BaseDistributionType, typename std::enable_if<is_bootstrap<D>::value, int>::type = 0

  //Generate the covariance matrix internally from boot-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, template<typename> class V, BOOTSTRAP_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, bootJackknifeDistribution<BaseNumericType,V>> &data_dj, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_dj.size();
    NumericSquareMatrix<BaseDistributionType> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov(i,j) = cov(j,i) = bootJackknifeDistribution<BaseNumericType,V>::covariance(data_dj.value(i), data_dj.value(j));

    importCovarianceMatrix(cov, cost_type);
  }

  //Generate the covariance matrix internally from bootstrap data. The resulting covariance matrix is "frozen", i.e. the same for all samples
  //Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, BOOTSTRAP_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data_j, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_j.size();
    auto binit = data_j.value(0).getInitializer();

    NumericSquareMatrix<BaseDistributionType> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++){
	double cv = BaseDistributionType::covariance(data_j.value(i), data_j.value(j));
	cov(i,j) = cov(j,i) = BaseDistributionType(cv, binit);
      }

    importCovarianceMatrix(cov, cost_type);
  }



  //Write the covariance matrix to a file in HDF5 format for external manipulation
  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    if(!have_corr_mat) error_exit(std::cout << "simpleFitWrapper::writeCovarianceMatrixHDF5  No covariance/correlation matrix available. Make sure you import one before calling this method!\n");
    NumericSquareMatrix<BaseDistributionType> cov = corr_mat;
    for(int i=0;i<cov.size();i++)
      for(int j=0;j<cov.size();j++)
	cov(i,j) = cov(i,j) * sigma[i] * sigma[j];
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
  
  inline const NumericSquareMatrix<BaseDistributionType> & getCorrelationMatrix() const{
    assert(have_corr_mat); return corr_mat;
  }
  inline const std::vector<BaseDistributionType> & getSigma() const{
    assert(have_corr_mat); return sigma;
  }
  
private:
  bool runSampleFit(ParameterType &params_s,
		    typename BaseDistributionType::DataType &chisq_s,
		    int &dof_s,
		    const correlationFunction<generalContainer, double> &data_s,
		    const int s,
		    const FitFuncFrozenBounded &fitfunc_s,
		    const NumericSquareMatrix<BaseDistributionType> &inv_corr_mat,
		    std::pair<double,int> *chisq_dof_nopriors_s_ptr) const{ 
    bool converged;
    switch(min_type){
    case MinimizerType::MarquardtLevenberg:
      chisq_s = simpleFitCommon::fitSampleML(converged, params_s, dof_s, data_s, sample(inv_corr_mat,s), sample(sigma,s), fitfunc_s, 
					     min_params, priors, chisq_dof_nopriors_s_ptr); break;
    case MinimizerType::GSLtrs:
      chisq_s = simpleFitCommon::fitSampleGSLtrs(converged, params_s, dof_s, data_s, sample(corr_mat,s), sample(sigma,s), fitfunc_s, 
						 min_params, priors, chisq_dof_nopriors_s_ptr); break;
    case MinimizerType::GSLmultimin:
      chisq_s = simpleFitCommon::fitSampleGSLmultimin(converged, params_s, dof_s, data_s, sample(inv_corr_mat,s), sample(sigma,s), fitfunc_s, 
						      min_params, priors, chisq_dof_nopriors_s_ptr); break;
    case MinimizerType::Minuit2:
      chisq_s = simpleFitCommon::fitSampleMinuit2(converged, params_s, dof_s, data_s, sample(inv_corr_mat,s), sample(sigma,s), fitfunc_s, 
						  min_params, priors, chisq_dof_nopriors_s_ptr); break;
    default:
      assert(0);
    }
    return converged;
  }

  inline FitFuncBounded setupParameterBounding(ParameterType &params_s) const{
    //Setup the sample fit function (frozen fit func needs a properly setup instance of the parameter type even if not freezing anything)
    FitFuncBounded fitfunc_bs(fitfunc);
    for(int t=0;t<bounded_trans.size();t++) fitfunc_bs.setBound(bounded_trans[t]);
      
    //Map guess parameters to internal params
    for(int t=0;t<bounded_trans.size();t++) bounded_trans[t].mapBoundedToUnbounded(params_s);

    return fitfunc_bs;
  }
  inline void convertBoundedParametersToExternalRep(ParameterType &params_s) const{
    //Map internal parameters to external params
    for(int t=0;t<bounded_trans.size();t++) bounded_trans[t].mapUnboundedToBounded(params_s);
  }  

  inline FitFuncFrozenBounded setupParameterFreezing(const FitFuncBounded &fitfunc_bs,
						     ParameterType &params_s, const int s) const{
    typedef iterate<BaseDistributionType> iter;
    FitFuncFrozenBounded fitfunc_s(fitfunc_bs);
    if(freeze_params.size() > 0){
      ParameterType frzp(params_s);
      for(int p=0;p<freeze_params.size();p++)
	frzp(freeze_params[p]) = iter::at(s, freeze_values[p]);
      fitfunc_s.freeze(freeze_params, frzp);
      params_s = fitfunc_s.mapParamsSupersetToSubset(params_s); //pull out the subset of parameters that are varied; here they are the same type
    }else{
      fitfunc_s.freeze({}, params_s);
    }
    return fitfunc_s;
  }
  inline void convertFrozenParametersToExternalRep(const FitFuncFrozenBounded &fitfunc_s, 
						   ParameterType &params_s) const{
    //Re-enlarge the parameter vector with the frozen fit parameters
    if(freeze_params.size() > 0) params_s = fitfunc_s.mapParamsSubsetToSuperset(params_s);
  }



  template<typename InputParameterType, typename GeneralizedCoordinate>
  void doSample(typename BaseDistributionType::template rebase<InputParameterType> &params,
		BaseDistributionType &chisq,
		BaseDistributionType &chisq_per_dof,
		int &dof,
		const int s,
		const correlationFunction<generalContainer, double> &data_sbase,
		const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data,
		const NumericSquareMatrix<BaseDistributionType> &inv_corr_mat,
		std::pair<BaseDistributionType, int>* chisq_dof_nopriors){
    typedef typename BaseDistributionType::template rebase<InputParameterType> ParameterDistributionType;
    typedef iterate<BaseDistributionType> iter;
    typedef iterate<ParameterDistributionType> iter_p;		

    //Get the sample data
    int ndata = data_sbase.size();
    correlationFunction<generalContainer, double> data_s = data_sbase;
    for(int i=0;i<ndata;i++) data_s.value(i) = iter::at(s, data.value(i));

    //Prepare the inner fit parameter type
    int nparam = iterate<ParameterDistributionType>::at(s, params).size();
    ParameterType params_s(nparam);
    for(int i=0;i<nparam;i++) params_s(i) = iter_p::at(s, params)(i);

    //Setup bounding including parameter transformation to internal representation
    FitFuncBounded fitfunc_bs = setupParameterBounding(params_s);

    //Setup freezing including parameter transformation to internal representation
    FitFuncFrozenBounded fitfunc_s = setupParameterFreezing(fitfunc_bs, params_s, s);

    //Prepare locations to store thread-local results
    int dof_s;
    std::pair<double,int> chisq_dof_nopriors_s;
    std::pair<double,int> *chisq_dof_nopriors_s_ptr = chisq_dof_nopriors != NULL ? &chisq_dof_nopriors_s : NULL;

    //Run the fitter
    auto & chisq_s = iter::at(s, chisq);
    bool converged = runSampleFit(params_s, chisq_s, dof_s, data_s, s, fitfunc_s, inv_corr_mat, chisq_dof_nopriors_s_ptr);

    //If you want the fitter to fail on non-convergence, set the associated minimizer parameter which typically defaults to do so
    if(!converged) std::cout << "Warning: thread " << omp_get_thread_num() << " fit did not converge on sample " << s << std::endl;

    //Re-enlarge the parameter vector with the frozen fit parameters
    convertFrozenParametersToExternalRep(fitfunc_s, params_s);

    //Map internal parameters to external params
    convertBoundedParametersToExternalRep(params_s);

    //Store outputs
    if(s==0) dof = dof_s;
      
    if(chisq_dof_nopriors){
      iter::at(s, chisq_dof_nopriors->first) = chisq_dof_nopriors_s.first;
      if(s==0) chisq_dof_nopriors->second = chisq_dof_nopriors_s.second;
    }

    iter::at(s, chisq_per_dof) = iter::at(s, chisq)/dof_s;
    for(int i=0;i<nparam;i++) iter_p::at(s, params)(i) = params_s(i);
  }
public:

  //Note the parameter type InputParameterType is translated internally into a parameterVector  (requires the usual size() and operator()(const int) methods)
  //The coordinate type is wrapped up in a generalContainer as this is only ever needed by the fit function (which knows what type it is and can retrieve it)
  //If chisq_dof_nopriors pointer is provided, the chisq computed without priors and the number of degrees of freedom without priors will be written there (distribution must have correct size
  template<typename InputParameterType, typename GeneralizedCoordinate>
  void fit(typename BaseDistributionType::template rebase<InputParameterType> &params,
	   BaseDistributionType &chisq,
	   BaseDistributionType &chisq_per_dof,
	   int &dof,
	   const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data,
	   std::pair<BaseDistributionType, int>* chisq_dof_nopriors = NULL){
    typedef typename BaseDistributionType::template rebase<InputParameterType> ParameterDistributionType;
    typedef iterate<BaseDistributionType> iter;
    typedef iterate<ParameterDistributionType> iter_p;

    if(!have_corr_mat) error_exit(std::cout << "simpleFitWrapper::fit  No covariance/correlation matrix available. Make sure you import one before calling this method!\n");

    int ndata = data.size();
    if(ndata == 0){
      std::cout << "Warning: Fit data container contains no data! Not performing a fit...." << std::endl;
      dof = -1;
      return;
    }
    int niter = iter::size(data.value(0));
    if(niter == 0){
      std::cout << "Warning: Fit data has 0 samples! Not performing a fit...." << std::endl;
      dof = -1;
      return;
    }

    assert(iter_p::size(params) == niter);
    assert(iter::size(chisq) == niter);
    assert(iter::size(chisq_per_dof) == niter);

    if(chisq_dof_nopriors) assert( iter::size(chisq_dof_nopriors->first) == niter );

    //Prepare inverse correlation matrix (the condition number is useful information even if we don't need the inverse explicitly)
    NumericSquareMatrix<BaseDistributionType> inv_corr_mat = invertCorrelationMatrix();

    correlationFunction<generalContainer, double> data_sbase(ndata);
    for(int i=0;i<ndata;i++) data_sbase.coord(i) = generalContainer(data.coord(i));

    //Minuit2 fitter annoyingly puts output on std::cout, which makes a mess with many threads; intercept and nullify cout output for all threads but 0
    std::streambuf* cout_rdbuf_orig = std::cout.rdbuf();
    thrbuf thr0_only(0,std::cout.rdbuf());
    
    if(min_type == MinimizerType::Minuit2) std::cout.rdbuf(&thr0_only);

    //Run the first iteration and use the result as a guess for the remainder, speeding up convergence
    doSample(params, chisq, chisq_per_dof, dof,
	     0, data_sbase, data, inv_corr_mat, chisq_dof_nopriors);

    //Run the rest of the samples
#pragma omp parallel for
    for(int s=1;s<niter;s++){
      iter_p::at(s,params) = iter_p::at(0,params);
      
      doSample(params, chisq, chisq_per_dof, dof,
	       s, data_sbase, data, inv_corr_mat, chisq_dof_nopriors);
    }

    if(min_type == MinimizerType::Minuit2) std::cout.rdbuf(cout_rdbuf_orig);

  }  
};

CPSFIT_END_NAMESPACE

#endif

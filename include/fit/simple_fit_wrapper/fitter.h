#ifndef _CPSFIT_SIMPLE_FIT_WRAPPER_FITTER_H
#define _CPSFIT_SIMPLE_FIT_WRAPPER_FITTER_H

#include<config.h>
#include<utils/macros.h>

#include<fit/simple_fit_wrapper/fitfunc_wrapper.h>
#include<minimizer/minimizer.h>
#include<minimizer/gsl_trs_minimizer.h>
#include<minimizer/gsl_multidim_minimizer.h>
#include<minimizer/minuit2_minimizer.h>
#include<fit/fitfunc/fitfunc_frozen.h>
#include<fit/fitfunc/fitfunc_bounded.h>
#include<tensors/numeric_square_matrix.h>
#include<distribution/jackknife.h>
#include<distribution/double_jackknife.h>
#include<distribution/block_double_jackknife.h>
#include<data_series/correlationfunction.h>
#include<fit/cost_function/correlated_chisq.h>
#include<fit/cost_function/correlated_chisq_terms.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(MinimizerType, (MarquardtLevenberg)(GSLtrs)(GSLmultimin)(Minuit2) );
GENERATE_ENUM_AND_PARSER(CostType, (Correlated)(Uncorrelated) );

class simpleFitWrapper{
  typedef genericFitFuncBase FitFunc;
  typedef BoundedFitFunc<genericFitFuncBase> FitFuncBounded;
  typedef FrozenFitFunc<FitFuncBounded> FitFuncFrozenBounded;

  const FitFunc &fitfunc;
  typedef FitFunc::ParameterType ParameterType;

  MinimizerType min_type;
  generalContainer min_params;

  NumericSquareMatrix<jackknifeDistribution<double> > corr_mat;
  std::vector<jackknifeDistribution<double>> sigma;
  bool have_corr_mat;

  std::vector<jackknifeDistribution<double> > freeze_values;
  std::vector<int> freeze_params;

  std::vector<boundedParameterTransform> bounded_trans;

  struct Prior{
    double value;
    double weight;
    int param_idx;
    Prior(const double value, const double weight, const int param_idx): value(value), weight(weight), param_idx(param_idx){}
  };
  std::vector<Prior> priors;

  template<typename CostFunc, typename FitFuncInternal>
  void addPriors(CostFunc &cost, const FitFuncInternal &fitfunc_s) const{
    for(int p=0;p<priors.size();p++){
      int subset_pidx = fitfunc_s.getParamsSubsetIndex(priors[p].param_idx);
      if(subset_pidx != -1)
	cost.addPrior(priors[p].value, priors[p].weight, subset_pidx);
    }
  }

  template<typename FitFuncInternal>
  double fitSampleML(ParameterType &params_s, int &dof,
		     const correlationFunction<generalContainer, double> &data_s,
		     const NumericSquareMatrix<double> &inv_corr_s,
		     const std::vector<double> &sigma_s,
		     const FitFuncInternal &fitfunc_s,
		     std::pair<double, int> *chisq_dof_nopriors = NULL) const{ //for frozen fit funcs, the values of the frozen parameters are sample dependent
    typedef CorrelatedChisqCostFunction<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;
    CostFunc cost(fitfunc_s, data_s, sigma_s, inv_corr_s);
    addPriors(cost, fitfunc_s);
    
    dof = cost.Ndof();
    
    typedef MarquardtLevenbergMinimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, inv_corr_s);
      chisq_dof_nopriors->first = cost_nopriors.cost(params_s);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }

    return chisq;
  }

  template<typename FitFuncInternal>
  double fitSampleGSLtrs(ParameterType &params_s, int &dof,
		     const correlationFunction<generalContainer, double> &data_s,
		     const NumericSquareMatrix<double> &corr_s,
		     const std::vector<double> &sigma_s,
		     const FitFuncInternal &fitfunc_s,
		     std::pair<double, int> *chisq_dof_nopriors = NULL) const{ //for frozen fit funcs, the values of the frozen parameters are sample dependent
    std::vector<double> corr_evals;
    std::vector<NumericVector<double> > corr_evecs;
    symmetricMatrixEigensolve(corr_evecs, corr_evals, corr_s);

    typedef CorrelatedChisqCostFunctionTerms<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;

    CostFunc cost(fitfunc_s, data_s, sigma_s, corr_evals, corr_evecs);
    addPriors(cost, fitfunc_s);
    
    dof = cost.Ndof();
    
    typedef GSLtrsMinimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, corr_evals, corr_evecs);
      NumericVector<double> cost_terms = cost_nopriors.costVector(params_s);
      chisq_dof_nopriors->first = 0.5*dot(cost_terms,cost_terms);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }
    
    return chisq;
  }

  template<typename FitFuncInternal>
  double fitSampleGSLmultimin(ParameterType &params_s, int &dof,
		     const correlationFunction<generalContainer, double> &data_s,
		     const NumericSquareMatrix<double> &inv_corr_s,
		     const std::vector<double> &sigma_s,
		     const FitFuncInternal &fitfunc_s,
		     std::pair<double, int> *chisq_dof_nopriors = NULL) const{ //for frozen fit funcs, the values of the frozen parameters are sample dependent
    typedef CorrelatedChisqCostFunction<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;
    CostFunc cost(fitfunc_s, data_s, sigma_s, inv_corr_s);
    addPriors(cost, fitfunc_s);
    
    dof = cost.Ndof();
    
    typedef GSLmultidimMinimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, inv_corr_s);
      chisq_dof_nopriors->first = cost_nopriors.cost(params_s);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }

    return chisq;
  }

  template<typename FitFuncInternal>
  double fitSampleMinuit2(ParameterType &params_s, int &dof,
			  const correlationFunction<generalContainer, double> &data_s,
			  const NumericSquareMatrix<double> &inv_corr_s,
			  const std::vector<double> &sigma_s,
			  const FitFuncInternal &fitfunc_s,
			  std::pair<double, int> *chisq_dof_nopriors = NULL) const{ //for frozen fit funcs, the values of the frozen parameters are sample dependent
#ifndef HAVE_MINUIT2
    error_exit(std::cout << "Library not compiled with Minuit2" << std::endl);
#else
    typedef CorrelatedChisqCostFunction<FitFuncInternal, correlationFunction<generalContainer, double> > CostFunc;
    CostFunc cost(fitfunc_s, data_s, sigma_s, inv_corr_s);
    addPriors(cost, fitfunc_s);
    
    dof = cost.Ndof();
    
    typedef Minuit2minimizer<CostFunc> Minimizer;
    typedef typename Minimizer::AlgorithmParameterType MParams;

    if(!min_params.is_null()) assert(min_params.is<MParams>());

    Minimizer min( cost, min_params.is_null() ? MParams() : min_params.value<MParams>());

    double chisq = min.fit(params_s);

    if(chisq_dof_nopriors){
      CostFunc cost_nopriors(fitfunc_s, data_s, sigma_s, inv_corr_s);
      chisq_dof_nopriors->first = cost_nopriors.cost(params_s);
      chisq_dof_nopriors->second = cost_nopriors.Ndof();
    }

    return chisq;
#endif
  }

  NumericSquareMatrix<jackknifeDistribution<double> > invertCorrelationMatrix() const{
    int nsample = corr_mat(0,0).size();
    NumericSquareMatrix<jackknifeDistribution<double> > inv_corr_mat(corr_mat);
    jackknifeDistribution<double> condition_number;
    svd_inverse(inv_corr_mat, corr_mat, condition_number);

    //Test the quality of the inverse
    NumericSquareMatrix<jackknifeDistribution<double> > test = corr_mat * inv_corr_mat;
    for(int i=0;i<test.size();i++) test(i,i) = test(i,i) - jackknifeDistribution<double>(nsample,1.0);    
    jackknifeDistribution<double> resid = modE(test);

    //Output the mean and standard deviation of the distributions of residual and condition number 
    std::cout << "Condition number = " << condition_number.mean() << " +- " << condition_number.standardError()/sqrt(nsample-1.) << std::endl;
    std::cout << "||CorrMat * CorrMat^{-1} - 1||_E = " << resid.mean() << " +- " << resid.standardError()/sqrt(nsample-1.) << std::endl;
    return inv_corr_mat;
  }

  static inline std::vector<double> sample(std::vector<jackknifeDistribution<double> > &v, const int s){
    std::vector<double> out(v.size()); for(int i=0;i<v.size();i++) out[i] = v[i].sample(s); 
    return out;
  }
  static inline NumericSquareMatrix<double> sample(NumericSquareMatrix<jackknifeDistribution<double> > &v, const int s){
    NumericSquareMatrix<double> out(v.size()); 
    for(int i=0;i<v.size();i++) 
      for(int j=0;j<v.size();j++) 
	out(i,j) = v(i,j).sample(s); 
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
	      const std::vector<jackknifeDistribution<double> > &_freeze_values){
    freeze_params = _freeze_params;
    freeze_values = _freeze_values;
  } 
  //Fix single parameter to input value
  inline void freeze(const int idx, const jackknifeDistribution<double> &val){
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
  void importCovarianceMatrix(const NumericSquareMatrix<jackknifeDistribution<double>> &cov, 
			      const CostType cost_type = CostType::Correlated){
    int nsample = cov(0,0).size();
    int ndata = cov.size();
    corr_mat.resize(ndata);
    sigma.resize(ndata);

    jackknifeDistribution<double> zero(nsample,0.);

    for(int i=0;i<ndata;i++) sigma[i] = sqrt(cov(i,i));
    
    for(int i=0;i<ndata;i++){
      corr_mat(i,i) = cov(i,i)/sigma[i]/sigma[i];
      for(int j=i+1;j<ndata;j++)
	corr_mat(i,j) = corr_mat(j,i) =  cost_type == CostType::Correlated ? cov(i,j)/sigma[i]/sigma[j] : zero;      
    }														

    have_corr_mat = true;
  }
  //Import a pre-generated correlation matrix and weights sigma   (sigma_i = sqrt(cov_ii))
  void importCorrelationMatrix(const NumericSquareMatrix<jackknifeDistribution<double>> &corr, const std::vector<jackknifeDistribution<double>> &sigma_in){
    corr_mat = corr;
    sigma = sigma_in;
    have_corr_mat = true;
  } 
  
  //Generate the covariance matrix internally from double-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename T>
  void generateCovarianceMatrix(const correlationFunction<T, doubleJackknifeDistribution<double>> &data_dj, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_dj.size();
    NumericSquareMatrix<jackknifeDistribution<double>> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov(i,j) = cov(j,i) = doubleJackknifeDistribution<double>::covariance(data_dj.value(i), data_dj.value(j));

    importCovarianceMatrix(cov, cost_type);
  }
  template<typename T>
  void generateCovarianceMatrix(const correlationFunction<T, blockDoubleJackknifeDistribution<double>> &data_dj, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_dj.size();
    NumericSquareMatrix<jackknifeDistribution<double>> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov(i,j) = cov(j,i) = blockDoubleJackknifeDistribution<double>::covariance(data_dj.value(i), data_dj.value(j));

    importCovarianceMatrix(cov, cost_type);
  }

  //Get the correlation matrix from the block double-jack and sigma from the regular, binned double-jackknife
  template<typename T>
  void generateCovarianceMatrix(const correlationFunction<T, doubleJackknifeDistribution<double>> &data_dj,
				const correlationFunction<T, blockDoubleJackknifeDistribution<double>> &data_bdj, 
				const CostType cost_type = CostType::Correlated){
    assert(data_dj.size() == data_bdj.size());
    assert(data_dj.value(0).size() == data_bdj.value(0).size()); //only inner indexing differs

    int ndata = data_dj.size();
    int nsample = data_bdj.value(0).size();

    NumericSquareMatrix<jackknifeDistribution<double>> cov_bdj(ndata);    
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++)
	cov_bdj(i,j) = cov_bdj(j,i) = blockDoubleJackknifeDistribution<double>::covariance(data_bdj.value(i), data_bdj.value(j));
    
    std::vector<jackknifeDistribution<double>> sigma_bdj(ndata);
    sigma.resize(ndata);
    for(int i=0;i<ndata;i++){
      sigma[i] = sqrt( doubleJackknifeDistribution<double>::covariance(data_dj.value(i), data_dj.value(i)) ); //sigma from regular dj
      sigma_bdj[i] = sqrt(cov_bdj(i,i));
    }
      
    //Corr from bdj
    corr_mat.resize(ndata);
    jackknifeDistribution<double> zero(nsample,0.);

    for(int i=0;i<ndata;i++){
      corr_mat(i,i) = cov_bdj(i,i)/sigma_bdj[i]/sigma_bdj[i];
      for(int j=i+1;j<ndata;j++)
	corr_mat(i,j) = corr_mat(j,i) =  cost_type == CostType::Correlated ? cov_bdj(i,j)/sigma_bdj[i]/sigma_bdj[j] : zero;      
    }
										       
    have_corr_mat = true;
  }


  //Generate the covariance matrix internally from single-jackknife data. The resulting covariance matrix is "frozen", i.e. the same for all samples
  //Option to use uncorrelated (diagonal) or correlated matrix
  template<typename T>
  void generateCovarianceMatrix(const correlationFunction<T, jackknifeDistribution<double>> &data_j, 
				const CostType cost_type = CostType::Correlated){
    int ndata = data_j.size();
    int nsample = data_j.value(0).size();
    NumericSquareMatrix<jackknifeDistribution<double>> cov(ndata);
    for(int i=0;i<ndata;i++)
      for(int j=i;j<ndata;j++){
	double cv = jackknifeDistribution<double>::covariance(data_j.value(i), data_j.value(j));
	cov(i,j) = cov(j,i) = jackknifeDistribution<double>(nsample, cv);
      }

    importCovarianceMatrix(cov, cost_type);
  }

  //Write the covariance matrix to a file in HDF5 format for external manipulation
  void writeCovarianceMatrixHDF5(const std::string &file) const{
#ifdef HAVE_HDF5
    if(!have_corr_mat) error_exit(std::cout << "simpleFitWrapper::writeCovarianceMatrixHDF5  No covariance/correlation matrix available. Make sure you import one before calling this method!\n");
    NumericSquareMatrix<jackknifeDistribution<double> > cov = corr_mat;
    for(int i=0;i<cov.size();i++)
      for(int j=0;j<cov.size();j++)
	cov(i,j) = cov(i,j) * sigma[i] * sigma[j];
    HDF5writer wr(file);
    write(wr, cov, "value");
#endif
  }
  
  inline const NumericSquareMatrix<jackknifeDistribution<double> > & getCorrelationMatrix() const{
    assert(have_corr_mat); return corr_mat;
  }
  inline const std::vector<jackknifeDistribution<double>> & getSigma() const{
    assert(have_corr_mat); return sigma;
  }
  
  //Note the parameter type P is translated internally into a parameterVector  (requires the usual size() and operator()(const int) methods)
  //The coordinate type is wrapped up in a generalContainer as this is only ever needed by the fit function (which knows what type it is and can retrieve it)
  //If chisq_dof_nopriors pointer is provided, the chisq computed without priors and the number of degrees of freedom without priors will be written there
  template<typename P, typename T>
  void fit(jackknifeDistribution<P> &params,
	   jackknifeDistribution<double> &chisq,
	   jackknifeDistribution<double> &chisq_per_dof,
	   int &dof,
	   const correlationFunction<T, jackknifeDistribution<double>> &data,
	   std::pair<jackknifeDistribution<double>, int>* chisq_dof_nopriors = NULL){
    if(!have_corr_mat) error_exit(std::cout << "simpleFitWrapper::fit  No covariance/correlation matrix available. Make sure you import one before calling this method!\n");

    int nsample = data.value(0).size();
    int ndata = data.size();

    if(chisq_dof_nopriors) chisq_dof_nopriors->first.resize(nsample);

    //Prepare inverse correlation matrix (the condition number is useful information even if we don't need the inverse explicitly)
    NumericSquareMatrix<jackknifeDistribution<double> > inv_corr_mat = invertCorrelationMatrix();

    correlationFunction<generalContainer, double> data_sbase(ndata);
    for(int i=0;i<ndata;i++) data_sbase.coord(i) = generalContainer(data.coord(i));

    //Minuit2 fitter annoyingly puts output on std::cout, which makes a mess with many threads; intercept and nullify cout output for all threads but 0
    std::streambuf* cout_rdbuf_orig = std::cout.rdbuf();
    thrbuf thr0_only(0,std::cout.rdbuf());
    
    if(min_type == MinimizerType::Minuit2) std::cout.rdbuf(&thr0_only);

#pragma omp parallel for
    for(int s=0;s<nsample;s++){
      //Get the sample data
      correlationFunction<generalContainer, double> data_s = data_sbase;
      for(int i=0;i<ndata;i++) data_s.value(i) = data.value(i).sample(s);

      //Prepare the inner fit parameter type
      int nparam = params.sample(s).size();
      ParameterType params_s(nparam);
      for(int i=0;i<nparam;i++) params_s(i) = params.sample(s)(i);

      //Setup the sample fit function (frozen fit func needs a properly setup instance of the parameter type even if not freezing anything)
      FitFuncBounded fitfunc_bs(fitfunc);
      for(int t=0;t<bounded_trans.size();t++) fitfunc_bs.setBound(bounded_trans[t]);
      
      //Map guess parameters to internal params
      for(int t=0;t<bounded_trans.size();t++) bounded_trans[t].mapBoundedToUnbounded(params_s);

      FitFuncFrozenBounded fitfunc_s(fitfunc_bs);
      if(freeze_params.size() > 0){
	ParameterType frzp(params_s);
	for(int p=0;p<freeze_params.size();p++)
	  frzp(freeze_params[p]) = freeze_values[p].sample(s);
	fitfunc_s.freeze(freeze_params, frzp);
	params_s = fitfunc_s.mapParamsSupersetToSubset(params_s); //pull out the subset of parameters that are varied; here they are the same type
      }else{
	fitfunc_s.freeze({}, params_s);
      }

      //Prepare locations to store thread-local results
      int dof_s;
      std::pair<double,int> chisq_dof_nopriors_s;
      std::pair<double,int> *chisq_dof_nopriors_s_ptr = chisq_dof_nopriors != NULL ? &chisq_dof_nopriors_s : NULL;

      //Run the fitter
      switch(min_type){
      case MinimizerType::MarquardtLevenberg:
	chisq.sample(s) = fitSampleML(params_s, dof_s, data_s, sample(inv_corr_mat,s), sample(sigma,s), fitfunc_s, chisq_dof_nopriors_s_ptr); break;
      case MinimizerType::GSLtrs:
	chisq.sample(s) = fitSampleGSLtrs(params_s, dof_s, data_s, sample(corr_mat,s), sample(sigma,s), fitfunc_s, chisq_dof_nopriors_s_ptr); break;
      case MinimizerType::GSLmultimin:
	chisq.sample(s) = fitSampleGSLmultimin(params_s, dof_s, data_s, sample(inv_corr_mat,s), sample(sigma,s), fitfunc_s, chisq_dof_nopriors_s_ptr); break;
      case MinimizerType::Minuit2:
	chisq.sample(s) = fitSampleMinuit2(params_s, dof_s, data_s, sample(inv_corr_mat,s), sample(sigma,s), fitfunc_s, chisq_dof_nopriors_s_ptr); break;
      default:
	assert(0);
      }
      
      //Re-enlarge the parameter vector with the frozen fit parameters
      if(freeze_params.size() > 0) params_s = fitfunc_s.mapParamsSubsetToSuperset(params_s);

      //Map internal parameters to external params
      for(int t=0;t<bounded_trans.size();t++) bounded_trans[t].mapUnboundedToBounded(params_s);
      
      //Store outputs
      if(s==0) dof = dof_s;
      
      if(chisq_dof_nopriors){
	chisq_dof_nopriors->first.sample(s) = chisq_dof_nopriors_s.first;
	if(s==0) chisq_dof_nopriors->second = chisq_dof_nopriors_s.second;
      }

      chisq_per_dof.sample(s) = chisq.sample(s)/dof_s;
      for(int i=0;i<nparam;i++) params.sample(s)(i) = params_s(i);
    }

    if(min_type == MinimizerType::Minuit2) std::cout.rdbuf(cout_rdbuf_orig);

  }  
};

CPSFIT_END_NAMESPACE

#endif

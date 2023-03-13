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

public:  
  //Constructor. Require fitfunc conforming to the wrapper in fitfunc_wrapper
  simpleFitWrapper(const FitFunc &fitfunc, 
		   const MinimizerType min_type, 
		   const generalContainer &min_params = generalContainer());
  
  //Can be changed at any time
  inline void setMinimizer(const MinimizerType _min_type, const generalContainer &_min_params= generalContainer());

  //Fix multiple parameters to input values
  inline void freeze(const std::vector<int> &_freeze_params,
		     const std::vector<BaseDistributionType> &_freeze_values);

  //Fix single parameter to input value
  inline void freeze(const int idx, const BaseDistributionType &val);

  //Remove all frozen parameters
  inline void resetFrozenParameters();

  //Add a Gaussian prior
  inline void addPrior(const double value, const double weight, const int param_idx);

  //Remove all priors
  inline void resetPriors();
  
  //Add bounds on parameters
  //min or max will be ignored as appropriate if using Max or Min bounds
  inline void setBound(int param, const ParameterBound bound, double min, double max);

  //Add bounds using a boundedParameterTransform object
  inline void setBound(const boundedParameterTransform &t);

  //Remove all bounds
  inline void resetBounds();

  //Import a pre-generated covariance matrix
  void importCovarianceMatrix(const NumericSquareMatrix<BaseDistributionType> &cov, 
			      const CostType cost_type = CostType::Correlated);

  //Import a pre-generated correlation matrix and weights sigma   (sigma_i = sqrt(cov_ii))
  void importCorrelationMatrix(const NumericSquareMatrix<BaseDistributionType> &corr, const std::vector<BaseDistributionType> &sigma_in);
  
#define JACKKNIFE_ONLY typename D=BaseDistributionType, typename std::enable_if<is_jackknife<D>::value, int>::type = 0

  //Generate the covariance matrix internally from double-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, template<typename> class V, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, doubleJackknifeDistribution<BaseNumericType,V>> &data_dj, 
				const CostType cost_type = CostType::Correlated);

  //Generate the covariance matrix internally from block double-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, template<typename> class V, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, blockDoubleJackknifeDistribution<BaseNumericType,V>> &data_bdj, 
				const CostType cost_type = CostType::Correlated);

  //Get the correlation matrix from the block double-jack and sigma from the regular, binned double-jackknife (the hybrid approach)
  template<typename GeneralizedCoordinate, template<typename> class V, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, doubleJackknifeDistribution<BaseNumericType,V>> &data_dj,
				const correlationFunction<GeneralizedCoordinate, blockDoubleJackknifeDistribution<BaseNumericType,V>> &data_bdj, 
				const CostType cost_type = CostType::Correlated);

  //Generate the covariance matrix internally from single-jackknife data. The resulting covariance matrix is "frozen", i.e. the same for all samples
  //Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, JACKKNIFE_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data_j, 
				const CostType cost_type = CostType::Correlated);

#define BOOTSTRAP_ONLY typename D=BaseDistributionType, typename std::enable_if<is_bootstrap<D>::value, int>::type = 0

  //Generate the covariance matrix internally from boot-jackknife data. Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, template<typename> class V, BOOTSTRAP_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, bootJackknifeDistribution<BaseNumericType,V>> &data_dj, 
                                const CostType cost_type = CostType::Correlated);

  //Generate the covariance matrix internally from bootstrap data. The resulting covariance matrix is "frozen", i.e. the same for all samples
  //Option to use uncorrelated (diagonal) or correlated matrix
  template<typename GeneralizedCoordinate, BOOTSTRAP_ONLY>
  void generateCovarianceMatrix(const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data_j, 
                                const CostType cost_type = CostType::Correlated);

  //Write the covariance matrix to a file in HDF5 format for external manipulation
  void writeCovarianceMatrixHDF5(const std::string &file) const;
  
  //Return the correlation matrix
  inline const NumericSquareMatrix<BaseDistributionType> & getCorrelationMatrix() const;

  //Return sigma, the weights of the covariance matrix
  inline const std::vector<BaseDistributionType> & getSigma() const;

  //Note the parameter type InputParameterType is translated internally into a parameterVector  (requires the usual size() and operator()(const int) methods)
  //The coordinate type is wrapped up in a generalContainer as this is only ever needed by the fit function (which knows what type it is and can retrieve it)
  //If chisq_dof_nopriors pointer is provided, the chisq computed without priors and the number of degrees of freedom without priors will be written there (distribution must have correct size
  template<typename InputParameterType, typename GeneralizedCoordinate>
  void fit(typename BaseDistributionType::template rebase<InputParameterType> &params,
	   BaseDistributionType &chisq,
	   BaseDistributionType &chisq_per_dof,
	   int &dof,
	   const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data,
	   std::pair<BaseDistributionType, int>* chisq_dof_nopriors = NULL);

  //A version of the above that takes the params as a vector of distributions
  template<typename GeneralizedCoordinate>
  void fit(std::vector<BaseDistributionType> &params,
	   BaseDistributionType &chisq,
	   BaseDistributionType &chisq_per_dof,
	   int &dof,
	   const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data,
	   std::pair<BaseDistributionType, int>* chisq_dof_nopriors = NULL);  
  
private:
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

  //Compute the inverse covariance matrix from the stored matrix
  NumericSquareMatrix<BaseDistributionType> invertCorrelationMatrix() const;

  //Extract a single sample from an array of distributions  
  static inline std::vector<BaseNumericType> sample(const std::vector<BaseDistributionType> &v, const int s);

  //Extract a single sample from a matrix of distributions
  static inline NumericSquareMatrix<BaseNumericType> sample(const NumericSquareMatrix<BaseDistributionType> &v, const int s);

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

  //Run the fit for a specific sample
  bool runSampleFit(ParameterType &params_s,
		    typename BaseDistributionType::DataType &chisq_s,
		    int &dof_s,
		    const correlationFunction<generalContainer, double> &data_s,
		    const int s,
		    const FitFuncFrozenBounded &fitfunc_s,
		    const NumericSquareMatrix<BaseDistributionType> &inv_corr_mat,
		    std::pair<double,int> *chisq_dof_nopriors_s_ptr) const;

  //Setup the parameter bounding
  inline FitFuncBounded setupParameterBounding(ParameterType &params_s) const;

  inline void convertBoundedParametersToExternalRep(ParameterType &params_s) const;

  //Setup the parameter freezing
  inline FitFuncFrozenBounded setupParameterFreezing(const FitFuncBounded &fitfunc_bs,
						     ParameterType &params_s, const int s) const;

  inline void convertFrozenParametersToExternalRep(const FitFuncFrozenBounded &fitfunc_s, 
						   ParameterType &params_s) const;

  //Implementation for specific sample
  template<typename InputParameterType, typename GeneralizedCoordinate>
  void doSample(typename BaseDistributionType::template rebase<InputParameterType> &params,
		BaseDistributionType &chisq,
		BaseDistributionType &chisq_per_dof,
		int &dof,
		const int s,
		const correlationFunction<generalContainer, double> &data_sbase,
		const correlationFunction<GeneralizedCoordinate, BaseDistributionType> &data,
		const NumericSquareMatrix<BaseDistributionType> &inv_corr_mat,
		std::pair<BaseDistributionType, int>* chisq_dof_nopriors);

};

#include "implementation/fitter_impl.tcc"

CPSFIT_END_NAMESPACE

#endif

#ifndef _SARLAC_BOOTSTRAP_PVALUE_H_
#define _SARLAC_BOOTSTRAP_PVALUE_H_

//An implementation of bootstrap procedure for determining a p-value in the context of a general fit
//q2 = "q^2", the generalization of \chi^2 for cases where the quantity isn't truly distributed by the \chi^2 distribution
#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>

#include<distribution/raw_data_distribution.h>
#include<distribution/bootstrap.h>
#include<random/rng.h>
#include<minimizer/minimizer.h>
#include<data_series/correlationfunction.h>
#include<fit/cost_function.h>

SARLAC_START_NAMESPACE

//Compute the bootstrap p-value given the distribution Dist
//Dist should be sorted in ascending order
double computePvalue(const double q2, const std::vector<double> &dist){
  int b_closest = -1;
  double q2diff = 1e12;
  int n = dist.size();
  for(int b=0;b<n;b++){
    double diff = fabs(q2 - dist[b]);
    if(diff < q2diff){
      q2diff = diff;
      b_closest = b;
    }else if(diff > q2diff) break; //because the values are sorted, once we get to the closest one the diff will start increasing again
  }
  return (n-b_closest-1)/double(n);
}

//Compute the boostrap p-value given the bootstrap distribution of q^2 computed from a fit to recentered data
double computePvalue(const double q2, const bootstrapDistribution<double> &dist){
  std::vector<double> q2_dist(dist.size()); for(int i=0;i<dist.size();i++) q2_dist[i] = dist.sample(i);
  std::sort(q2_dist.begin(), q2_dist.end(), [&](const double a, const double b){ return a<b; });
  return computePvalue(q2, q2_dist);
}

//Note, fit_data_cen are the data that was used to obtain the fit parameters. This does not have to be the same size as raw_data nor does its coordinate type have to be the same
//The details of how the raw data are converted into fit data is left up to the user
//fit_values should be the fit to the original data evaluated at each coordinate


//FitFunctor should be an object with method 
//double operator()(const RawDataContainer &raw_data, const correlationFunction<GeneralizedCoordinate, DataType> &corrections, const int b) const
//where raw_data are the data to be used and corrections are the shifts to the fit data that should be applied by the user prior to fitting (recentering)
//Note b is the bootstrap sample index. It doesn't have to be used but the user may find it useful in the storage of any results produced from the bootstrap ensembles

//nthread = -1 uses all available threads
//RawDataContainer is some generic structure containing the raw data

//Resampler is a Functor that acts on a raw data container and scrambles its samples according to the map provided. It should have:
//void operator()(RawDataContainer &raw_data, const std::vector<int> &map)   
//Be aware: if nsample is not a multiple of the block size it will be truncated to the nearest multiple

//If q2_boot_p is provided the q^2 distribution will be copied to this address
template<typename RawDataContainer, typename GeneralizedCoordinate, typename ValueType, typename Resampler, typename FitFunctor>
double bootstrapPvalue(const double q2,
		       const RawDataContainer &raw_data, int nsample,
		       const correlationFunction<GeneralizedCoordinate, ValueType> &fit_data_cen,
		       const correlationFunction<GeneralizedCoordinate, ValueType> &fit_values,	
		       const Resampler &resampler, FitFunctor &fitter, const std::vector<std::vector<int> > &resample_table, 
		       int nthread = -1, std::vector<double> *q2_boot_p = NULL){
  const int p_fit = fit_data_cen.size();
  assert(fit_values.size() == p_fit);

  if(nthread = -1) nthread = omp_get_max_threads();

  int nboot = resample_table.size();
  assert(resample_table[0].size() <= nsample);
  nsample = resample_table[0].size();
  
  correlationFunction<GeneralizedCoordinate, ValueType> corrections(p_fit);
  for(int i=0;i<p_fit;i++){
    assert(fit_values.coord(i) == fit_data_cen.coord(i));
    corrections.coord(i) = fit_values.coord(i);
    corrections.value(i) = fit_values.value(i) - fit_data_cen.value(i);
  }

  std::vector<RawDataContainer> thr_boot_data(nthread, raw_data);
  std::vector<double> q2_boot(nboot);

#pragma omp parallel for
  for(int b=0;b<nboot;b++){
    int me = omp_get_thread_num();
    RawDataContainer &boot_data = thr_boot_data[me];
    boot_data = raw_data;    
    resampler(boot_data, resample_table[b]);
    q2_boot[b] = fitter(boot_data, corrections, b);
  }
  std::sort(q2_boot.begin(), q2_boot.end(), [&](const double a, const double b){ return a<b; });

  if(q2_boot_p) *q2_boot_p = q2_boot;

  return computePvalue(q2, q2_boot);
}

template<typename RawDataContainer, typename GeneralizedCoordinate, typename ValueType, typename Resampler, typename FitFunctor>
double bootstrapPvalue(const double q2,
		       const RawDataContainer &raw_data, int nsample,
		       const correlationFunction<GeneralizedCoordinate, ValueType> &fit_data_cen,
		       const correlationFunction<GeneralizedCoordinate, ValueType> &fit_values,	
		       const Resampler &resampler, FitFunctor &fitter, const int nboot = 1000, 
		       const BootResampleTableType table_type = BootResampleTableType::NonOverlappingBlock,
		       const int block_size = 1, int nthread = -1, RNGstore &rng = RNG,
		       std::vector<double> *q2_boot_p = NULL){
  std::vector<std::vector<int> > otable = generateResampleTable(nsample, nboot, table_type, block_size, rng);  //[b][s]

  return bootstrapPvalue(q2, raw_data, nsample, fit_data_cen, fit_values, resampler, fitter, otable,
			 nthread, q2_boot_p);
}


//Standard fits where the raw data and fit data have the same coordinate type and the fit data are just the means of the raw data
template<typename FitFunc>
struct BasicFitFunctor{
  typedef typename FitFunc::ParameterType ParameterType;
  typedef typename FitFunc::GeneralizedCoordinate GeneralizedCoordinate;
  typedef typename FitFunc::ValueType DataType;

  typedef CorrelatedChisqCostFunction<FitFunc, correlationFunction<GeneralizedCoordinate, DataType> > CostFunc;
  typedef MarquardtLevenbergMinimizer<CostFunc> Minimizer;

  typedef typename correlationFunction<GeneralizedCoordinate, DataType>::ElementType ElementType;

  ParameterType guess;
  MarquardtLevenbergParameters<double> min_params;

  const FitFunc &fitfunc;

  BasicFitFunctor(const ParameterType &guess, const FitFunc &fitfunc, const MarquardtLevenbergParameters<double> &min_params = MarquardtLevenbergParameters<double>()):
    guess(guess), fitfunc(fitfunc), min_params(min_params){}

  void computeCov(std::vector<DataType> &sigma, NumericSquareMatrix<DataType> &inv_corr,
		  const correlationFunction<GeneralizedCoordinate, rawDataDistribution<DataType> > &data) const{
    int p = data.size();
  
    //Get covariance matrix from sample covariance
    NumericSquareMatrix<DataType> cov(p);
    sigma = std::vector<DataType>(p);
    for(int i=0;i<p;i++){
      for(int j=i; j<p; j++)
	cov(i,j) = cov(j,i) = rawDataDistribution<DataType>::covariance(data.value(i), data.value(j));
      sigma[i] = sqrt(cov(i,i));
    }
    NumericSquareMatrix<DataType> corr(p);
    for(int i=0;i<p;i++)
      for(int j=i; j<p; j++)
	corr(i,j) = corr(j,i) = cov(i,j)/sigma[i]/sigma[j];
    inv_corr = NumericSquareMatrix<DataType>(p);
  
    double cond;
    svd_inverse(inv_corr, corr, cond);
  
    if(log10(fabs(cond)) > 10){
      error_exit(std::cout << "Large condition number " << cond << " " << log10(fabs(cond)) << std::endl);
    }
  }


  double operator()(const correlationFunction<GeneralizedCoordinate, rawDataDistribution<DataType> > &raw_data, const correlationFunction<GeneralizedCoordinate, DataType> &corrections, const int b) const{
    const int p = raw_data.size(); assert(corrections.size() == p);
    correlationFunction<GeneralizedCoordinate, DataType> data_cen(raw_data.size());
    for(int i=0;i<p;i++){
      assert(raw_data.coord(i) == corrections.coord(i));
      data_cen.coord(i) = raw_data.coord(i);
      data_cen.value(i) = raw_data.value(i).mean() + corrections.value(i);
    }

    std::vector<DataType> sigma;
    NumericSquareMatrix<DataType> inv_corr;  
    computeCov(sigma, inv_corr, raw_data);
	       
    CostFunc cost(fitfunc, data_cen, sigma, inv_corr);

    Minimizer min(cost, min_params);

    ParameterType params = guess;
    return min.fit(params);
  }
};

template<typename T>
void bootstrapResampleRaw(rawDataDistribution<T> &raw, const std::vector<int> &map){
  assert(map.size() <= raw.size());
  raw = rawDataDistribution<T>(map.size(), [&](const int s){ return raw.sample(map[s]); });
}

template<typename C,typename T>
struct CorrFuncResampler{
  void operator()(correlationFunction<C, rawDataDistribution<T> > &raw, const std::vector<int> &map) const{
    for(int t=0;t<raw.size();t++) bootstrapResampleRaw(raw.value(t), map);
  }
};

template<typename GeneralizedCoordinate, typename ValueType, typename FitFunctor>
inline double bootstrapPvalue(const double q2,
			      const correlationFunction<GeneralizedCoordinate, rawDataDistribution<ValueType> > &raw_data,
			      const correlationFunction<GeneralizedCoordinate, ValueType> &fit_data_cen,
			      const correlationFunction<GeneralizedCoordinate, ValueType> &fit_values,	
			      const FitFunctor &fitter, const int nboot = 1000, const BootResampleTableType table_type = BootResampleTableType::NonOverlappingBlock,

			      const int block_size = 1, int nthread = -1, RNGstore &rng = RNG){
  return bootstrapPvalue(q2, raw_data, raw_data.value(0).size(), fit_data_cen, fit_values ,CorrFuncResampler<GeneralizedCoordinate, ValueType>(), fitter, nboot, table_type, block_size, nthread, rng);
}

//This function performs the "recentering" for correlation functions of distribution data
template<typename GeneralizedCoordinate, typename DistributionType>
void recenter(correlationFunction<GeneralizedCoordinate, DistributionType> &data, const correlationFunction<GeneralizedCoordinate, typename DistributionType::DataType> &corrections){
  assert(data.size() == corrections.size());
  for(int i=0;i<data.size();i++){
    assert(data.coord(i) == corrections.coord(i));
    for(int s=0;s<iterate<DistributionType>::size(data.value(i));s++){
      auto &v = iterate<DistributionType>::at(s, data.value(i));
      v = v + corrections.value(i);
    }
  }
}




SARLAC_END_NAMESPACE
#endif

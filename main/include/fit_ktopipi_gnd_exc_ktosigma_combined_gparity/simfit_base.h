#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_BASE_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_BASE_H

#include<config.h>
#include<utils/macros.h>
#include<pipi_common/analyze_chisq.h>
#include "simfit_common.h"

CPSFIT_START_NAMESPACE

template<typename FitFunc>
void analyzeChisqFF(const correlationFunction<SimFitCoordGen,  jackknifeDistributionD> &corr_comb_j,
		    const jackknifeDistribution<typename FitFunc::ParameterType> &params, const FitFunc &fitfunc,
		    const std::map<std::unordered_map<std::string, std::string> const*, std::string> &pmap_descr){
  struct PP{
    typedef std::map<std::unordered_map<std::string, std::string> const*, std::string> const* PtrType;
    inline static PtrType & descr(){ static PtrType p; return p; }

    inline static void print(std::ostream &os, const SimFitCoordGen &c){ os << printCoord(c, *descr()); }
    inline static std::string typeInfo(const SimFitCoordGen &c){ return descr()->find(c.param_map)->second; }   
  };
  PP::descr() = &pmap_descr;
  
  AnalyzeChisq<FitFunc,PP> chisq_analyze(corr_comb_j, fitfunc, params);
  chisq_analyze.printChisqContribs(Correlation);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Correlation);
  chisq_analyze.printChisqContribs(Covariance);
  chisq_analyze.examineProjectedDeviationContribsEvalNorm(Covariance);
}

template< template<typename, template<typename> class> class DistributionType > 
struct ResampledDataContainers{};

template< template<typename, template<typename> class> class DistributionType > 
struct SimFitDataContainers{};

template< template<typename, template<typename> class> class DistributionType > 
struct simultaneousFitBase: public simultaneousFitCommon{
  COPY_COMMON_TYPEDEFS;

  virtual void load2ptFitParams(const std::vector<PiPiOperator> &operators, const InputParamArgs &iargs, const int dist_size) = 0;

  virtual std::vector<DistributionType<Params, basic_vector> > fit(const ResampledDataContainers<DistributionType> &fit_data,
								   const std::vector<PiPiOperator> &operators,
								   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated,
								   const CovarianceMatrix covariance_matrix) = 0;

  template<typename FitFunc>
  static void runfit(std::vector<DistributionType<Params, basic_vector> > &params,
		     const SimFitDataContainers<DistributionType> &fit_data,
		     const subsetMapDescr &pmap_descr,
		     const FitFunc &fitfunc,
		     const std::vector<int> &freeze_params,
		     const DistributionType<Params, basic_vector> &freeze_vals,
		     const Params &guess,
		     const bool correlated, const CovarianceMatrix covariance_matrix){

    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, correlatedFitPolicy>::type FitPolicies;

    typedef DistributionType<double, basic_vector> DistributionTypeD;
    const int nQ = fit_data.getNq();
    params.resize(nQ);

    auto init = fit_data.getDistributionInitializer();

    std::vector<DistributionTypeD> chisq(nQ, DistributionTypeD(init)), 
      chisq_per_dof(nQ, DistributionTypeD(init)), 
      pvalue(nQ, DistributionTypeD(init));

    for(int q=0;q<nQ;q++){
      params[q] = DistributionType<Params, basic_vector>(init, guess);

      std::cout << "Performing " << q+1 << " fit" << std::endl;
      fitter<FitPolicies> fit;
      fit.importFitFunc(fitfunc);
      fit.freeze(freeze_params, freeze_vals);
      
      importCostFunctionParameters<correlatedFitPolicy, FitPolicies> import;
      fit_data.generateCovarianceMatrix(import, fit, covariance_matrix, q);

      if(!correlated) import.setUncorrelated();
          
      int ndof;
      fit.fit(params[q], chisq[q], chisq_per_dof[q], ndof, fit_data.getFitData(q));

      pvalue[q] = DistributionTypeD(init);
      for(int i=0; i<iterate<DistributionTypeD>::size(pvalue[q]); i++) 
	iterate<DistributionTypeD>::at(i, pvalue[q]) = chiSquareDistribution::pvalue(ndof, iterate<DistributionTypeD>::at(i, chisq[q]) );

      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:\n" << params[q] << std::endl;
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value(chi^2): " << pvalue[q] << std::endl;

      // std::cout << "Analysis of contributions to chi^2" << std::endl;      
      // analyzeChisqFF<typename FitPolicies::baseFitFunc>(A0_sim_j[q],params[q],fitfunc,pmap_descr);
    }  
    for(int q=0;q<nQ;q++){
      std::cout << "Q" << q+1 << std::endl;
      std::cout << "Params:\n" << params[q] << std::endl;
      std::cout << "Chisq: " << chisq[q] << std::endl;
      std::cout << "Chisq/dof: " << chisq_per_dof[q] << std::endl;
      std::cout << "p-value: " << pvalue[q] << std::endl;
    }
    writeParamsStandard(params, "params.hdf5");
    writeParamsStandard(chisq, "chisq.hdf5");
    writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
    writeParamsStandard(pvalue, "pvalue.hdf5");
  }


  virtual ~simultaneousFitBase(){}
};

CPSFIT_END_NAMESPACE

#endif

#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_BOOTSTRAP_CONTAINERS_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_BOOTSTRAP_CONTAINERS_H

#include<config.h>
#include<utils/macros.h>

#include "simfit_base.h"

CPSFIT_START_NAMESPACE


template<>
struct ResampledDataContainers<bootstrapDistribution>{
  const ResampledData<bootstrapDistributionD> &data_b;
  const ResampledData<bootJackknifeDistributionD> &data_bj;
  
  ResampledDataContainers(const ResampledData<bootstrapDistributionD> &data_b,
			  const ResampledData<bootJackknifeDistributionD> &data_bj): data_b(data_b), data_bj(data_bj){}

  const ResampledData<bootstrapDistributionD> & getFitData() const{ return data_b; }
};



template<>
struct SimFitDataContainers<bootstrapDistribution>{
  COPY_COMMON_TYPEDEFS;

  typedef correlationFunction<SimFitCoordGen, bootstrapDistributionD> SimFitCorrFuncBoot;
  typedef correlationFunction<SimFitCoordGen, bootJackknifeDistributionD> SimFitCorrFuncBJack;

  std::vector<SimFitCorrFuncBoot> A0_sim_b;
  std::vector<SimFitCorrFuncBJack> A0_sim_bj;
  
  int getNq() const{
    int nq = A0_sim_b.size();
    return nq;
  }

  template<typename FitPolicies>
  void generateCovarianceMatrix(importCostFunctionParameters<correlatedFitPolicy, FitPolicies> &import, fitter<FitPolicies> &fit, const CovarianceMatrix covariance_matrix, const int q) const{
    import.import(fit, A0_sim_bj[q]);
  }

  bootstrapInitType getDistributionInitializer() const{
    return A0_sim_b[0].value(0).getInitializer();
  }

  void generateSimData(const ResampledDataContainers<bootstrapDistribution> &fit_data,
		       const std::vector<PiPiOperator> &operators, const int tmin_k_op, const int tmin_op_snk,
		       const operatorSubsetMap &op_param_maps, const subsetMapDescr &pmap_descr, const CovarianceMatrix covariance_matrix){
    simultaneousFitCommon::generateSimData(A0_sim_b, fit_data.data_b, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    simultaneousFitCommon::generateSimData(A0_sim_bj, fit_data.data_bj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
  }
  
  const SimFitCorrFuncBoot & getFitData(const int q) const{ return A0_sim_b[q]; }

  //Func:  void [](auto &corrfunc, const int q)
  template<typename Func, typename CorrelationFunctionType>
  void applyFunctionToInternalCorrFunc(const Func &func, CorrelationFunctionType &corrfunc){
    for(int q=0;q<corrfunc.size();q++)
      func(corrfunc[q],q);
  }

  //Func:  void [](auto &corrfunc, const int q)
  template<typename Func>
  void applyFunctionToCorrFunc(const Func &func){
    applyFunctionToInternalCorrFunc(func, A0_sim_b);
    applyFunctionToInternalCorrFunc(func, A0_sim_bj);
  }

};

CPSFIT_END_NAMESPACE

#endif

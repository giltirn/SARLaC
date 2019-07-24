#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_JACKKNIFE_CONTAINERS_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_SIMFIT_JACKKNIFE_CONTAINERS_H

#include<config.h>
#include<utils/macros.h>

#include "simfit_base.h"

CPSFIT_START_NAMESPACE


template<>
struct ResampledDataContainers<jackknifeDistribution>{
  const ResampledData<jackknifeDistributionD> &data_j;
  const ResampledData<doubleJackknifeA0StorageType> &data_dj;
  const ResampledData<blockDoubleJackknifeA0StorageType> &data_bdj;
  
  ResampledDataContainers(const ResampledData<jackknifeDistributionD> &data_j,
			  const ResampledData<doubleJackknifeA0StorageType> &data_dj,
			  const ResampledData<blockDoubleJackknifeA0StorageType> &data_bdj): data_j(data_j), data_dj(data_dj), data_bdj(data_bdj){}

  const ResampledData<jackknifeDistributionD> & getFitData() const{ return data_j; }
};



template<>
struct SimFitDataContainers<jackknifeDistribution>{
  COPY_COMMON_TYPEDEFS;

  typedef correlationFunction<SimFitCoordGen, jackknifeDistributionD> SimFitCorrFuncJack;
  typedef correlationFunction<SimFitCoordGen, doubleJackknifeDistributionD> SimFitCorrFuncDJack;
  typedef correlationFunction<SimFitCoordGen, blockDoubleJackknifeDistributionD> SimFitCorrFuncBDJack;

  std::vector<SimFitCorrFuncJack> A0_sim_j;
  std::vector<SimFitCorrFuncDJack> A0_sim_dj;
  std::vector<SimFitCorrFuncBDJack> A0_sim_bdj;
  
  int getNq() const{
    int nq = A0_sim_j.size();
    return nq;
  }

  template<typename FitPolicies>
  void generateCovarianceMatrix(importCostFunctionParameters<correlatedFitPolicy, FitPolicies> &import, fitter<FitPolicies> &fit, const CovarianceMatrix covariance_matrix, const int q) const{
    bool do_dj = covariance_matrix == CovarianceMatrix::Regular;
    bool do_bdj = covariance_matrix == CovarianceMatrix::Block;

    if(do_dj) import.import(fit, A0_sim_dj[q]);
    if(do_bdj) import.import(fit, A0_sim_bdj[q]);
  }

  int getDistributionInitializer() const{
    return A0_sim_j[0].value(0).size();
  }

  void generateSimData(const ResampledDataContainers<jackknifeDistribution> &fit_data,
		       const std::vector<PiPiOperator> &operators, const int tmin_k_op, const int tmin_op_snk,
		       const operatorSubsetMap &op_param_maps, const subsetMapDescr &pmap_descr, const CovarianceMatrix covariance_matrix){
    bool do_dj = covariance_matrix == CovarianceMatrix::Regular;
    bool do_bdj = covariance_matrix == CovarianceMatrix::Block;
    simultaneousFitCommon::generateSimData(A0_sim_j, fit_data.data_j, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    if(do_dj) simultaneousFitCommon::generateSimData(A0_sim_dj, fit_data.data_dj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
    if(do_bdj) simultaneousFitCommon::generateSimData(A0_sim_bdj, fit_data.data_bdj, operators, tmin_k_op, tmin_op_snk, op_param_maps, pmap_descr);
  }
  
  const SimFitCorrFuncJack & getFitData(const int q) const{ return A0_sim_j[q]; }

  //Func:  void [](auto &corrfunc, const int q)
  template<typename Func, typename CorrelationFunctionType>
  void applyFunctionToInternalCorrFunc(const Func &func, CorrelationFunctionType &corrfunc){
    for(int q=0;q<corrfunc.size();q++)
      func(corrfunc[q],q);
  }

  //Func:  void [](auto &corrfunc, const int q)
  template<typename Func>
  void applyFunctionToCorrFunc(const Func &func){
    applyFunctionToInternalCorrFunc(func, A0_sim_j);
    applyFunctionToInternalCorrFunc(func, A0_sim_dj);
    applyFunctionToInternalCorrFunc(func, A0_sim_bdj);
  }

};

CPSFIT_END_NAMESPACE

#endif

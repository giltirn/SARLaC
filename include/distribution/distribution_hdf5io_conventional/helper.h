#ifndef _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_HELPER_H_
#define _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_HELPER_H_

//These helper classes implement certain actions for mapping between Distribution<compound-type> and vector<Distribution<POD-type>
///They should be specialized if the actions need to be customized

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>
#include<distribution/distribution_hdf5io_conventional/utils.h>

CPSFIT_START_NAMESPACE

//For a generic vector-type distribution
template<typename DistributionOfPODtype, typename DistributionOfStructType>
struct standardIOhelper{
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    const int nsample = from.size();
    to.resize(from.size());
    for(int s=0;s<nsample;s++) to.sample(s) = from.sample(s)(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample){
    to.resize(nsample);
    typedef typename DistributionOfStructType::DataType StructType;
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<nsample;s++) to.sample(s) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    const int nsample = from.size();
    for(int s=0;s<nsample;s++){
      to.sample(s)(param_idx) = from.sample(s);
    }
  }
};

CPSFIT_END_NAMESPACE
#endif

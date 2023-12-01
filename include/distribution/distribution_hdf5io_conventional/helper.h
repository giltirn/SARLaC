#ifndef _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_HELPER_H_
#define _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_HELPER_H_

//These helper classes implement certain actions for mapping between Distribution<compound-type> and vector<Distribution<POD-type>
///They should be specialized if the actions need to be customized

#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>
#include<distribution/distribution_hdf5io_conventional/utils.h>

SARLAC_START_NAMESPACE

template<typename D> struct iterate;

//For a generic vector-type distribution
template<typename DistributionOfPODtype, typename DistributionOfStructType>
struct standardIOhelper{
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    typedef iterate<DistributionOfPODtype> iter_to;
    typedef iterate<DistributionOfStructType> iter_from;

    to.resize(from.getInitializer());
    for(int s=0;s<iter_from::size(from);s++) iter_to::at(s, to) = iter_from::at(s, from)(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const typename DistributionOfStructType::initType &init){
    to.resize(init);
    typedef typename DistributionOfStructType::DataType StructType;
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);

    typedef iterate<DistributionOfStructType> iter_to;
    for(int s=0;s<iter_to::size(to);s++) iter_to::at(s, to) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    typedef iterate<DistributionOfPODtype> iter_from;
    typedef iterate<DistributionOfStructType> iter_to;

    for(int s=0;s<iter_from::size(from);s++)
      iter_to::at(s, to)(param_idx) = iter_from::at(s, from);    
  }
};

SARLAC_END_NAMESPACE
#endif

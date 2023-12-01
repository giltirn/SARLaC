#ifndef _DOUBLE_JACKKNIFE_DIST_HDF5IO_H_
#define _DOUBLE_JACKKNIFE_DIST_HDF5IO_H_


#include<config.h>
#include<distribution/distribution_hdf5io_conventional/helper.h>
#include<distribution/distribution_hdf5io_conventional/io_format.h>
#include<distribution/double_jackknife/class.h>

SARLAC_START_NAMESPACE

template<typename PODtype, typename StructType, template<typename> class V1, template<typename> class V2>
struct standardIOhelper<doubleJackknifeDistribution<PODtype,V1>, doubleJackknifeDistribution<StructType,V2> >{
  typedef doubleJackknifeDistribution<PODtype,V1> DistributionOfPODtype;
  typedef doubleJackknifeDistribution<StructType,V2> DistributionOfStructType;
  
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    const int nsample = from.size();
    to.resize(from.size());
    for(int s=0;s<nsample;s++)
      for(int t=0;t<nsample-1;t++)      
	to.sample(s).sample(t) = from.sample(s).sample(t)(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample){
    to.resize(nsample);
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<nsample;s++)
      for(int t=0;t<nsample-1;t++)
	to.sample(s).sample(t) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    const int nsample = from.size();
    for(int s=0;s<nsample;s++){
      for(int t=0;t<nsample-1;t++){
	to.sample(s).sample(t)(param_idx) = from.sample(s).sample(t);
      }
    }
  }
};


template<typename T, template<typename> class V>
struct getDataType<doubleJackknifeDistribution<T,V>,1>{ typedef T type; };

template<typename T, template<typename> class V>
struct standardIOformat<doubleJackknifeDistribution<T,V> >{
  inline static void write(HDF5writer &writer, const std::vector<doubleJackknifeDistribution<T,V> > &value, const std::string &tag){ SARLaC::write(writer,value,tag); } //no option to flatten for double jack
  inline static void read(HDF5reader &reader, std::vector<doubleJackknifeDistribution<T,V> > &value, const std::string &tag){ SARLaC::read(reader,value,tag); }
  inline static void write(HDF5writer &writer, const std::vector<std::vector<doubleJackknifeDistribution<T,V> > > &value, const std::string &tag){ SARLaC::write(writer,value,tag); }
  inline static void read(HDF5reader &reader, std::vector<std::vector<doubleJackknifeDistribution<T,V> > > &value, const std::string &tag){ SARLaC::read(reader,value,tag); }  
};

SARLAC_END_NAMESPACE
#endif

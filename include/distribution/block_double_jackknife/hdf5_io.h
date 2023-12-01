#ifndef _BLOCK_DOUBLE_JACKKNIFE_DIST_HDF5IO_H_
#define _BLOCK_DOUBLE_JACKKNIFE_DIST_HDF5IO_H_


#include<config.h>
#include<distribution/distribution_hdf5io_conventional/helper.h>
#include<distribution/distribution_hdf5io_conventional/io_format.h>
#include<distribution/block_double_jackknife/class.h>

SARLAC_START_NAMESPACE

template<typename T, template<typename> class V>
void write(HDF5writer &writer, const blockDoubleJackknifeDistribution<T,V> &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.sampleVector(),"data");
  write(writer,value.nSamplesUnbinned(),"nsamples_unbinned");
  write(writer,value.binSize(),"bin_size");
  writer.leave();
}
template<typename T, template<typename> class V>
void read(HDF5reader &reader, blockDoubleJackknifeDistribution<T,V> &value, const std::string &tag){
  reader.enter(tag);
  int nsample, bin_size;
  read(reader,nsample,"nsamples_unbinned");
  read(reader,bin_size,"bin_size");
  value.resize(nsample, bin_size);
  read(reader,value.sampleVector(),"data");
  reader.leave();
}




template<typename PODtype, typename StructType, template<typename> class V1, template<typename> class V2>
struct standardIOhelper<blockDoubleJackknifeDistribution<PODtype,V1>, blockDoubleJackknifeDistribution<StructType,V2> >{
  typedef blockDoubleJackknifeDistribution<PODtype,V1> DistributionOfPODtype;
  typedef blockDoubleJackknifeDistribution<StructType,V2> DistributionOfStructType;
  
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    to.resize(from.nSamplesUnbinned(), from.binSize());
    for(int s=0;s<from.size();s++)
      for(int t=0;t<from.sample(s).size();t++)      
	to.sample(s).sample(t) = from.sample(s).sample(t)(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample, const int bin_size){
    to.resize(nsample, bin_size);
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<to.size();s++)
      for(int t=0;t<to.sample(s).size();t++)
	to.sample(s).sample(t) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    for(int s=0;s<from.size();s++){
      for(int t=0;t<from.sample(s).size();t++){
	to.sample(s).sample(t)(param_idx) = from.sample(s).sample(t);
      }
    }
  }
};


template<typename T, template<typename> class V>
struct getDataType<blockDoubleJackknifeDistribution<T,V>,1>{ typedef T type; };

template<typename T, template<typename> class V>
struct standardIOformat<blockDoubleJackknifeDistribution<T,V> >{
  inline static void write(HDF5writer &writer, const std::vector<blockDoubleJackknifeDistribution<T,V> > &value, const std::string &tag){ SARLaC::write(writer,value,tag); } //no option to flatten for double jack
  inline static void read(HDF5reader &reader, std::vector<blockDoubleJackknifeDistribution<T,V> > &value, const std::string &tag){ SARLaC::read(reader,value,tag); }
  inline static void write(HDF5writer &writer, const std::vector<std::vector<blockDoubleJackknifeDistribution<T,V> > > &value, const std::string &tag){ SARLaC::write(writer,value,tag); }
  inline static void read(HDF5reader &reader, std::vector<std::vector<blockDoubleJackknifeDistribution<T,V> > > &value, const std::string &tag){ SARLaC::read(reader,value,tag); }  
};



//NOTE: Should think about adding ability to write flattened double-jackknife
// template<typename T, template<typename> class V> //override for distributions with extra information
// struct extraFlattenIO<blockDoubleJackknifeDistribution<T,V> >{
//   static inline void write(HDF5writer &writer, const std::vector<blockDoubleJackknifeDistribution<T> > &value){
//     std::vector<int> nsample(value.size());
//     std::vector<int> bin_size(value.size());
//     for(int i=0;i<value.size();i++){ nsample[i] = value[i].nSamplesUnbinned(); bin_size[i] = value[i].binSize(); }
//     SARLaC::write(writer,nsample,"nsample");
//     SARLaC::write(writer,bin_size,"bin_size");
//   }    
//   static inline void read(HDF5reader &reader, std::vector<superJackknifeDistribution<T> > &value){
//     std::vector<int> nsample;
//     std::vector<int> bin_size;
//     SARLaC::read(reader,nsample,"nsample");
//     SARLaC::read(reader,bin_size,"bin_size");
//     assert(value.size() == nsample.size());
//     assert(value.size() == bin_size.size());
//     for(int i=0;i<value.size();i++) value[i].resize(nsample[i], bin_size[i]);
//   }
//   static inline void write(HDF5writer &writer, const std::vector<std::vector<superJackknifeDistribution<T> > > &value){
//     int sz = 0; for(int i=0;i<value.size();i++) sz += value[i].size();
//     std::vector<int> nsample(sz);
//     std::vector<int> bin_size(sz);
//     int k=0;
//     for(int i=0;i<value.size();i++){
//       for(int j=0;j<value[i].size();j++){ 
// 	nsample[k] = value[i][j].nSamplesUnbinned(); 
// 	bin_size[k] = value[i][j].binSize(); 
// 	++k;
//       }
//     }

//     SARLaC::write(writer,nsample,"nsample");
//     SARLaC::write(writer,bin_size,"bin_size");
//   }
//   static inline void read(HDF5reader &reader, std::vector<std::vector<superJackknifeDistribution<T> > > &value){
//     std::vector<int> nsample;
//     std::vector<int> bin_size;
//     SARLaC::read(reader,nsample,"nsample");
//     SARLaC::read(reader,bin_size,"bin_size");

//     int k=0;
//     for(int i=0;i<value.size();i++){
//       for(int j=0;j<value[i].size();j++){ 
// 	value[i][j].resize(nsample[k],bin_size[k]);
// 	++k;
//       }
//     }
//   }
// };


SARLAC_END_NAMESPACE
#endif

#ifndef _BOOT_JACKKNIFE_DIST_HDF5IO_H_
#define _BOOT_JACKKNIFE_DIST_HDF5IO_H_


#include<config.h>
#include<distribution/distribution_hdf5io_conventional/helper.h>
#include<distribution/distribution_hdf5io_conventional/io_format.h>
#include<distribution/boot_jackknife/class.h>

CPSFIT_START_NAMESPACE

#ifdef HAVE_HDF5

template<typename T, template<typename> class V>
void write(HDF5writer &writer, const bootJackknifeDistribution<T,V> &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.origEnsJackknife(),"base_jack");
  write(writer,value.sampleVector(),"data");
  write(writer,value.confidence(),"confidence");
  writer.leave();
}
template<typename T, template<typename> class V>
void read(HDF5reader &reader, bootJackknifeDistribution<T,V> &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.origEnsJackknife(),"base_jack");
  read(reader,value.sampleVector(),"data");
  read(reader,value.confidence(),"confidence");
  reader.leave();
}

#endif

template<typename PODtype, typename StructType, template<typename> class V1, template<typename> class V2>
struct standardIOhelper<bootJackknifeDistribution<PODtype,V1>, bootJackknifeDistribution<StructType,V2> >{
  typedef bootJackknifeDistribution<PODtype,V1> DistributionOfPODtype;
  typedef bootJackknifeDistribution<StructType,V2> DistributionOfStructType;
  
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    to.resize(from.getInitializer());

    for(int s=0;s<from.size();s++)
      for(int t=0;t<from.sample(s).size();t++)      
	to.sample(s).sample(t) = from.sample(s).sample(t)(param_idx);
    
    for(int t=0;t<from.origEnsJackknife().size();t++)      
      to.origEnsJackknife().sample(t) = from.origEnsJackknife().sample(t)(param_idx);

  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const bootJackknifeInitType &init){
    to.resize(init);
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<to.size();s++)
      for(int t=0;t<to.sample(s).size();t++)
	to.sample(s).sample(t) = base;

    for(int t=0;t<to.origEnsJackknife().size();t++) 
      to.origEnsJackknife().sample(t) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    for(int s=0;s<to.size();s++){
      for(int t=0;t<to.sample(s).size();t++){
	to.sample(s).sample(t)(param_idx) = from.sample(s).sample(t);
      }
    }
    for(int t=0;t<to.origEnsJackknife().size();t++) 
      to.origEnsJackknife().sample(t)(param_idx) = from.origEnsJackknife().sample(t);

  }
};


template<typename T, template<typename> class V>
struct getDataType<bootJackknifeDistribution<T,V>,1>{ typedef T type; };

template<typename T, template<typename> class V>
struct standardIOformat<bootJackknifeDistribution<T,V> >{
  inline static void write(HDF5writer &writer, const std::vector<bootJackknifeDistribution<T,V> > &value, const std::string &tag){ CPSfit::write(writer,value,tag); } //no option to flatten for boot jack
  inline static void read(HDF5reader &reader, std::vector<bootJackknifeDistribution<T,V> > &value, const std::string &tag){ CPSfit::read(reader,value,tag); }
  inline static void write(HDF5writer &writer, const std::vector<std::vector<bootJackknifeDistribution<T,V> > > &value, const std::string &tag){ CPSfit::write(writer,value,tag); }
  inline static void read(HDF5reader &reader, std::vector<std::vector<bootJackknifeDistribution<T,V> > > &value, const std::string &tag){ CPSfit::read(reader,value,tag); }  
};

CPSFIT_END_NAMESPACE
#endif

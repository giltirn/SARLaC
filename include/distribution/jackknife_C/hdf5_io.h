#ifndef _JACKKNIFEC_HDF5_IO_H_
#define _JACKKNIFEC_HDF5_IO_H_

#include<config.h>

#include<distribution/jackknife_C/class.h>
#include<distribution/distribution_hdf5io_basic.h>
#include<distribution/distribution_hdf5io_conventional/helper.h>

SARLAC_START_NAMESPACE

#ifdef HAVE_HDF5

//Specialization for jackknifeCdistribution
template<typename T, template<typename> class V>
void write(HDF5writer &writer, const jackknifeCdistribution<T,V> &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.cen,"cen");
  write(writer,value._data,"data");
  writer.leave();
}
template<typename T, template<typename> class V>
void read(HDF5reader &reader, jackknifeCdistribution<T,V> &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.cen,"cen");
  read(reader,value._data,"data");  
  reader.leave();
}


template<typename T, template<typename> class V> //override for distributions with extra information
struct extraFlattenIO<jackknifeCdistribution<T,V> >{
  static inline void write(HDF5writer &writer, const std::vector<jackknifeCdistribution<T,V> > &value){
    std::vector<T> cen(value.size());
    for(int i=0;i<cen.size();i++) cen[i] = value[i].best();
    SARLaC::write(writer,cen,"central");
  }    
  static inline void read(HDF5reader &reader, std::vector<jackknifeCdistribution<T,V> > &value){
    std::vector<T> cen;
    SARLaC::read(reader,cen,"central");
    for(int i=0;i<cen.size();i++) value[i].best() = cen[i];
  }
  static inline void write(HDF5writer &writer, const std::vector<std::vector<jackknifeCdistribution<T,V> > > &value){
    int sz = 0; for(int i=0;i<value.size();i++) sz += value[i].size();
    
    std::vector<T> cen(sz);
    int off=0;
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++)
	cen[off++] = value[i][j].best();

    SARLaC::write(writer,cen,"central");
  }
  static inline void read(HDF5reader &reader, std::vector<std::vector<jackknifeCdistribution<T,V> > > &value){
    std::vector<T> cen;
    SARLaC::read(reader,cen,"central");

    int off=0;
    for(int i=0;i<value.size();i++) //has already been resized
      for(int j=0;j<value[i].size();j++)
	value[i][j].best() = cen[off++];
  } 
};

#endif //HAVE_HDF5

//Specialize the helper class for conversion to and from vector<Distribution<POD-type> >
template<typename PODtype, typename StructType, template<typename> class V1, template<typename> class V2>
struct standardIOhelper<jackknifeCdistribution<PODtype,V1>, jackknifeCdistribution<StructType,V2> >{
  typedef jackknifeCdistribution<PODtype,V1> DistributionOfPODtype;
  typedef jackknifeCdistribution<StructType,V2> DistributionOfStructType;

  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    const int nsample = from.size();
    to.resize(from.size());
    for(int s=0;s<nsample;s++) to.sample(s) = from.sample(s)(param_idx);
    to.best() = from.best()(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample){
    to.resize(nsample);
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<nsample;s++) to.sample(s) = base;
    to.best() = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    const int nsample = from.size();
    for(int s=0;s<nsample;s++){
      to.sample(s)(param_idx) = from.sample(s);
    }
    to.best()(param_idx) = from.best();
  }
};



SARLAC_END_NAMESPACE



#endif

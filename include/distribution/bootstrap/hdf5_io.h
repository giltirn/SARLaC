#ifndef _BOOTSTRAP_HDF5_IO_H_
#define _BOOTSTRAP_HDF5_IO_H_

#include<config.h>

#include<distribution/bootstrap/class.h>
#include<distribution/distribution_hdf5io_basic.h>
#include<distribution/distribution_hdf5io_conventional/helper.h>

SARLAC_START_NAMESPACE

#ifdef HAVE_HDF5

//Specialization for bootstrapDistribution
template<typename T, template<typename> class V>
void write(HDF5writer &writer, const bootstrapDistribution<T,V> &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value._confidence,"confidence");
  write(writer,value.avg,"avg");
  write(writer,value._data,"data");
  writer.leave();
}
template<typename T, template<typename> class V>
void read(HDF5reader &reader, bootstrapDistribution<T,V> &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value._confidence,"confidence");
  read(reader,value.avg,"avg");
  read(reader,value._data,"data");  
  reader.leave();
}


template<typename T, template<typename> class V> //override for distributions with extra information
struct extraFlattenIO<bootstrapDistribution<T,V> >{
  static inline void write(HDF5writer &writer, const std::vector<bootstrapDistribution<T,V> > &value){
    std::vector<int> conf(value.size());
    for(int i=0;i<conf.size();i++) conf[i] = value[i].confidence();
    SARLaC::write(writer,conf,"confidence");

    std::vector<T> avg(value.size());
    for(int i=0;i<avg.size();i++) avg[i] = value[i].best();
    SARLaC::write(writer,avg,"avg");
  }    
  static inline void read(HDF5reader &reader, std::vector<bootstrapDistribution<T,V> > &value){
    std::vector<int> conf;
    SARLaC::read(reader,conf,"confidence");
    for(int i=0;i<conf.size();i++)
      value[i].confidence() = conf[i];

    std::vector<T> avg;
    SARLaC::read(reader,avg,"avg");
    for(int i=0;i<avg.size();i++) value[i].best() = avg[i];
  }
  static inline void write(HDF5writer &writer, const std::vector<std::vector<bootstrapDistribution<T,V> > > &value){
    int sz = 0; for(int i=0;i<value.size();i++) sz += value[i].size();

    std::vector<int> conf(sz);
    std::vector<T> avg(sz);
    int off=0;
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++){
	avg[off] = value[i][j].best();
	conf[off] = value[i][j].confidence();
	++off;
      }

    SARLaC::write(writer,conf,"confidence");
    SARLaC::write(writer,avg,"avg");
  }
  static inline void read(HDF5reader &reader, std::vector<std::vector<bootstrapDistribution<T,V> > > &value){
    std::vector<int> conf;
    std::vector<T> avg;
    SARLaC::read(reader,conf,"confidence");
    SARLaC::read(reader,avg,"avg");

    int off=0;
    for(int i=0;i<value.size();i++) //has already been resized
      for(int j=0;j<value[i].size();j++){
	value[i][j].best() = avg[off];
	value[i][j].confidence() = conf[off];
	++off;
      }
  } 
};

#endif //HAVE_HDF5

SARLAC_END_NAMESPACE



#endif

#ifndef _JACKKNIFE_DIST_XMLIO_H_
#define _JACKKNIFE_DIST_XMLIO_H_

//For XML IO compatible with UKfit "boot_print" results

#include<config.h>
#include<distribution/jackknife/class.h>
#include<serialize/xml_serialize.h>

SARLAC_START_NAMESPACE

template<typename DistributionType>
struct UKvalenceDistributionContainer{
  int Nentries;
  std::vector<DistributionType> list;
};
template<typename DistributionType>
void read(XMLreader &reader, UKvalenceDistributionContainer<DistributionType> &v, const std::string &tag){
  reader.enter(tag);
  read(reader,v.Nentries,"Nentries");
  read(reader,v.list,"list");
  reader.leave();
}

template<typename T, template<typename> class V>
void read(XMLreader &reader, jackknifeDistribution<T,V> &v, const std::string &tag){
  reader.enter(tag);
#define GETIT(type,name) type name; read(reader,name,#name)

  GETIT(std::string, SampleType);
  assert(SampleType == "Jackknife");
  
  GETIT(int, Nmeas);
  GETIT(double, avg);
  GETIT(std::vector<T>, values);
  
  if(values.size() != Nmeas) error_exit(std::cout << "read(XMLreader &, jackknifeDistribution<T> &, const std::string &) file states Nmeas=" << Nmeas << " but read " << values.size() << " samples!\n");
  v.resize(Nmeas);
  for(int i=0;i<Nmeas;i++) v.sample(i) = values[i];
  
#undef GETIT
  reader.leave();
}

SARLAC_END_NAMESPACE
#endif

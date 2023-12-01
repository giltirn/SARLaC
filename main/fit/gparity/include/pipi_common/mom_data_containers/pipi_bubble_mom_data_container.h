#pragma once
#include "mom_data_container_common.h"
#include "../base_data_containers/pipi_bubble_data_container.h"

SARLAC_START_NAMESPACE

template<typename _ContainerType,typename Extra=empty_t>
class bubbleDataAllMomentaBase: public Extra{
  friend Extra;
public:
  typedef _ContainerType ContainerType;
  typedef std::pair<SourceOrSink, threeMomentum> keyType;
  typedef typename std::map<keyType, ContainerType>::iterator iterator;  
  typedef typename std::map<keyType, ContainerType>::const_iterator const_iterator;  
  typedef typename ContainerType::DistributionType DistributionType;
private:
  std::map<keyType, ContainerType> B;
  int Lt;
  int tsep_pipi;
  
  ContainerType & get(const keyType &key, bool lock){
    auto it = B.find(key);
    if(it == B.end()){
      if(lock) error_exit(std::cout << "bubbleDataAllMomenta::get Could not find requested source/sink and momentum:" << (int)key.first << " " << key.second << std::endl );

      it = B.insert(std::make_pair(key, ContainerType())).first;
      it->second.setup(key.first,Lt,tsep_pipi);
    }
    return it->second;
  }

public:
  bubbleDataAllMomentaBase(const int _Lt, const int tsep_pipi): Lt(_Lt), tsep_pipi(tsep_pipi){}
  bubbleDataAllMomentaBase(){}

  void setup(const int _Lt, const int _tsep_pipi){
    Lt = _Lt; tsep_pipi = _tsep_pipi;
  }
  
  const ContainerType &operator()(const SourceOrSink src_snk, const threeMomentum &mom) const{
    bubbleDataAllMomentaBase<ContainerType,Extra> *t = const_cast<bubbleDataAllMomentaBase<ContainerType,Extra> *>(this);
    return const_cast<const ContainerType &>( t->get(keyType(src_snk,mom),true) );
  }
  ContainerType &operator()(const SourceOrSink src_snk, const threeMomentum &mom){
    return this->get(keyType(src_snk,mom),false);
  }
  const ContainerType &operator()(const keyType &key) const{
    bubbleDataAllMomentaBase<ContainerType,Extra> *t = const_cast<bubbleDataAllMomentaBase<ContainerType,Extra> *>(this);
    return const_cast<const ContainerType &>( t->get(key,true) );
  }
  ContainerType &operator()(const keyType &key){
    return this->get(key,false);
  }

  inline int getLt() const{ return Lt; }
  inline int getTsepPiPi() const{ return tsep_pipi; }

  inline const_iterator begin() const{ return B.begin(); }
  inline const_iterator end() const{ return B.end(); }

  inline iterator begin(){ return B.begin(); }
  inline iterator end(){ return B.end(); }

  inline size_t getNmomenta() const{ return B.size(); }

  GENERATE_HDF5_SERIALIZE_METHOD((B)(Lt)(tsep_pipi));
};
#ifdef HAVE_HDF5
template<typename C,typename E>
inline void write(HDF5writer &writer, const bubbleDataAllMomentaBase<C,E> &d, const std::string &tag){ d.write(writer,tag); }
template<typename C,typename E>
inline void read(HDF5reader &reader, bubbleDataAllMomentaBase<C,E> &d, const std::string &tag){ d.read(reader,tag); }
#endif

template<typename BubbleDataContainerType>
struct bubbleDataAllMomentaExtra{
  typedef bubbleDataAllMomentaBase<BubbleDataContainerType, bubbleDataAllMomentaExtra<BubbleDataContainerType> > full_type;
  inline full_type & upcast(){ return *static_cast< full_type* >(this); }

  inline void bin(const int bin_size){
    full_type & me = upcast();
    for(auto it = me.begin(); it != me.end(); it++)
	it->second.bin(bin_size);  
  }
};


typedef bubbleDataAllMomentaBase<bubbleData, bubbleDataAllMomentaExtra<bubbleData> > bubbleDataAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataZ, bubbleDataAllMomentaExtra<bubbleDataZ> > bubbleDataAllMomentaZ;
typedef bubbleDataAllMomentaBase<bubbleDataJack> bubbleDataJackAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataDoubleJack> bubbleDataDoubleJackAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataBlockDoubleJack> bubbleDataBlockDoubleJackAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataBoot> bubbleDataBootAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataBootJack> bubbleDataBootJackAllMomenta;

template<typename DistributionType>
struct _bubbleDataAllMomentaTypeSelector{};
template<>
struct _bubbleDataAllMomentaTypeSelector<rawDataDistributionD>{ typedef bubbleDataAllMomenta type; };
template<>
struct _bubbleDataAllMomentaTypeSelector<rawDataDistribution<std::complex<double> > >{ typedef bubbleDataAllMomentaZ type; };
template<>
struct _bubbleDataAllMomentaTypeSelector<jackknifeDistributionD>{ typedef bubbleDataJackAllMomenta type; };
template<>
struct _bubbleDataAllMomentaTypeSelector<doubleJackknifeDistributionD>{ typedef bubbleDataDoubleJackAllMomenta type; };
template<>
struct _bubbleDataAllMomentaTypeSelector<blockDoubleJackknifeDistributionD>{ typedef bubbleDataBlockDoubleJackAllMomenta type; };
template<>
struct _bubbleDataAllMomentaTypeSelector<bootstrapDistributionD>{ typedef bubbleDataBootAllMomenta type; };
template<>
struct _bubbleDataAllMomentaTypeSelector<bootJackknifeDistributionD>{ typedef bubbleDataBootJackAllMomenta type; };


template<typename DistributionType>
using bubbleDataAllMomentaSelect = typename _bubbleDataAllMomentaTypeSelector<DistributionType>::type;

inline bubbleDataAllMomenta reIm(const bubbleDataAllMomentaZ &in, const int reim){
  bubbleDataAllMomenta out(in.getLt(), in.getTsepPiPi());
  for(auto it = in.begin(); it != in.end(); ++it)
    out(it->first) = reIm(it->second, reim);
  return out;
}

SARLAC_END_NAMESPACE

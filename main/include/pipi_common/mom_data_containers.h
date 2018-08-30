#ifndef PIPI_MOM_DATA_CONTAINERS_H
#define PIPI_MOM_DATA_CONTAINERS_H

#include <boost/serialization/map.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/timer/timer.hpp>

#include<config.h>
#include<utils/macros.h>

#include "data_containers.h"
#include "threemomentum.h"

CPSFIT_START_NAMESPACE

template<typename _ContainerType, typename Extra = empty_t>
class figureDataAllMomentaBase: public Extra{
  friend Extra;  
public:
  typedef _ContainerType ContainerType;
  typedef typename std::map<sinkSourceMomenta, ContainerType>::const_iterator const_iterator;  
  typedef typename std::map<sinkSourceMomenta, ContainerType>::iterator iterator;  
private:
  typedef std::map<sinkSourceMomenta, ContainerType> MapType;
  MapType C;
  MapType D;
  MapType R;
  MapType V;
  int Lt;
  int Nsample;

  MapType* getMap(const char fig){
    switch(fig){
    case 'C':
      return &C;
    case 'D':
      return &D;
    case 'R':
      return &R;
    case 'V':
      return &V;
    default:
      error_exit(std::cout << "figureDataAllMomenta::getMap invalid figure " << fig << std::endl);
    }
  }
  MapType const* getMap(const char fig) const{
    return const_cast<MapType const*>( const_cast<figureDataAllMomentaBase<_ContainerType,Extra>* >(this)->getMap(fig) );
  }
  ContainerType & get(const char fig, const sinkSourceMomenta &mom, bool lock){
    MapType* mp = getMap(fig);
    typename MapType::iterator it = mp->find(mom);
    if(it == mp->end()){
      if(lock) error_exit(std::cout << "figureDataAllMomenta::get Could not find requested momentum\n");

      it = mp->insert(std::make_pair(mom, ContainerType())).first;
      it->second.setup(Lt, Nsample);
    }
    return it->second;
  }

public:
  figureDataAllMomentaBase(const int _Lt, const int _Nsample): Lt(_Lt), Nsample(_Nsample){}
  figureDataAllMomentaBase(){}

  void setup(const int _Lt, const int _Nsample){
    Lt = _Lt; Nsample = _Nsample;
  }
  
  const ContainerType &operator()(const char fig, const sinkSourceMomenta &mom) const{
    figureDataAllMomentaBase<ContainerType,Extra> *t = const_cast<figureDataAllMomentaBase<ContainerType,Extra> *>(this);
    return const_cast<const ContainerType &>( t->get(fig,mom,true) );
  }
  ContainerType &operator()(const char fig, const sinkSourceMomenta &mom){
    return this->get(fig,mom,false);
  }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return Nsample; }

  const_iterator begin(const char fig) const{
    return getMap(fig)->begin();
  }
  const_iterator end(const char fig) const{
    return getMap(fig)->end();
  }
  iterator begin(const char fig){
    return getMap(fig)->begin();
  }
  iterator end(const char fig){
    return getMap(fig)->end();
  }
  GENERATE_HDF5_SERIALIZE_METHOD((C)(D)(R)(V)(Lt)(Nsample));
};

#ifdef HAVE_HDF5
template<typename C,typename E>
inline void write(HDF5writer &writer, const figureDataAllMomentaBase<C,E> &d, const std::string &tag){ d.write(writer,tag); }
template<typename C,typename E>
inline void read(HDF5reader &reader, figureDataAllMomentaBase<C,E> &d, const std::string &tag){ d.read(reader,tag); }
#endif


struct figureDataAllMomentaExtra{
  typedef figureDataAllMomentaBase<figureData, figureDataAllMomentaExtra> full_type;
  inline full_type & upcast(){ return *static_cast< full_type* >(this); }

  inline void bin(const int bin_size){
    full_type & me = upcast();
    static const char figs[4] = {'C','D','R','V'};
    for(int f=0;f<4;f++)
      for(full_type::iterator it = me.begin(figs[f]); it != me.end(figs[f]); it++)
	it->second.bin(bin_size);  
    me.Nsample /= bin_size;
  }
};

typedef figureDataAllMomentaBase<figureData, figureDataAllMomentaExtra> figureDataAllMomenta;
typedef figureDataAllMomentaBase<figureDataDoubleJack> figureDataDoubleJackAllMomenta;


template<typename _ContainerType,typename Extra=empty_t>
class bubbleDataAllMomentaBase: public Extra{
  friend Extra;
public:
  typedef _ContainerType ContainerType;
  typedef std::pair<SourceOrSink, threeMomentum> keyType;
  typedef typename std::map<keyType, ContainerType>::iterator iterator;  
  typedef typename std::map<keyType, ContainerType>::const_iterator const_iterator;  
private:
  std::map<keyType, ContainerType> B;
  int Lt;
  int Nsample;
  int tsep_pipi;
  
  ContainerType & get(const keyType &key, bool lock){
    auto it = B.find(key);
    if(it == B.end()){
      if(lock) error_exit(std::cout << "bubbleDataAllMomenta::get Could not find requested source/sink and momentum:" << (int)key.first << " " << key.second << std::endl );

      it = B.insert(std::make_pair(key, ContainerType())).first;
      it->second.setup(key.first,Lt,tsep_pipi,Nsample);
    }
    return it->second;
  }

public:
  bubbleDataAllMomentaBase(const int _Lt, const int tsep_pipi, const int _Nsample): Lt(_Lt), tsep_pipi(tsep_pipi), Nsample(_Nsample){}
  bubbleDataAllMomentaBase(){}

  void setup(const int _Lt, const int _tsep_pipi, const int _Nsample){
    Lt = _Lt; tsep_pipi = _tsep_pipi; Nsample = _Nsample;
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
  inline int getNsample() const{ return Nsample; }
  inline int getTsepPiPi() const{ return tsep_pipi; }

  inline const_iterator begin() const{ return B.begin(); }
  inline const_iterator end() const{ return B.end(); }

  inline iterator begin(){ return B.begin(); }
  inline iterator end(){ return B.end(); }

  GENERATE_HDF5_SERIALIZE_METHOD((B)(Lt)(Nsample)(tsep_pipi));
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
    me.Nsample /= bin_size;
  }
};


typedef bubbleDataAllMomentaBase<bubbleData, bubbleDataAllMomentaExtra<bubbleData> > bubbleDataAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataZ, bubbleDataAllMomentaExtra<bubbleDataZ> > bubbleDataAllMomentaZ;
typedef bubbleDataAllMomentaBase<bubbleDataDoubleJack> bubbleDataDoubleJackAllMomenta;

void saveHDF5checkpoint(const figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const std::string &file){
  (std::cout << "Saving HDF5 data checkpoint\n").flush(); boost::timer::auto_cpu_timer t("Report: Saved HDF5 data checkpoint in %w s\n");
  HDF5writer wr(file);
  write(wr,raw_data,"raw_data");
  write(wr,raw_bubble_data,"raw_bubble_data");
}
void loadHDF5checkpoint(figureDataAllMomenta &raw_data, bubbleDataAllMomenta &raw_bubble_data, const std::string &file){
  (std::cout << "Loading HDF5 data checkpoint\n").flush(); boost::timer::auto_cpu_timer t("Report: Loaded HDF5 data checkpoint in %w s\n");
  HDF5reader rd(file);
  read(rd,raw_data,"raw_data");
  read(rd,raw_bubble_data,"raw_bubble_data");
}

inline bubbleDataAllMomenta reIm(const bubbleDataAllMomentaZ &in, const int reim){
  bubbleDataAllMomenta out(in.getLt(), in.getTsepPiPi(), in.getNsample());
  for(auto it = in.begin(); it != in.end(); ++it)
    out(it->first) = reIm(it->second, reim);
  return out;
}

CPSFIT_END_NAMESPACE

#endif

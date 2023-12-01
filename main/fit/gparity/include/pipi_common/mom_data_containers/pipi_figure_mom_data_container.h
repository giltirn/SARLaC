#pragma once
#include "mom_data_container_common.h"
#include "../base_data_containers/figure_data_container.h"

SARLAC_START_NAMESPACE

template<typename _ContainerType, typename Extra = empty_t>
class figureDataAllMomentaBase: public Extra{
  friend Extra;  
public:
  typedef _ContainerType ContainerType;
  typedef typename std::map<sinkSourceMomenta, ContainerType>::const_iterator const_iterator;  
  typedef typename std::map<sinkSourceMomenta, ContainerType>::iterator iterator;
  typedef typename ContainerType::DistributionType DistributionType;
private:
  typedef std::map<sinkSourceMomenta, ContainerType> MapType;
  MapType C;
  MapType D;
  MapType R;
  MapType V;
  int Lt;

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
      it->second.setup(Lt);
    }
    return it->second;
  }

public:
  figureDataAllMomentaBase(const int _Lt): Lt(_Lt){}
  figureDataAllMomentaBase(){}

  void setup(const int _Lt){
    Lt = _Lt;
  }
  
  const ContainerType &operator()(const char fig, const sinkSourceMomenta &mom) const{
    figureDataAllMomentaBase<ContainerType,Extra> *t = const_cast<figureDataAllMomentaBase<ContainerType,Extra> *>(this);
    return const_cast<const ContainerType &>( t->get(fig,mom,true) );
  }
  ContainerType &operator()(const char fig, const sinkSourceMomenta &mom){
    return this->get(fig,mom,false);
  }
  inline int getLt() const{ return Lt; }

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
  GENERATE_HDF5_SERIALIZE_METHOD((C)(D)(R)(V)(Lt));
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
  }
};

typedef figureDataAllMomentaBase<figureData, figureDataAllMomentaExtra> figureDataAllMomenta;
typedef figureDataAllMomentaBase<figureDataJack> figureDataJackAllMomenta;
typedef figureDataAllMomentaBase<figureDataDoubleJack> figureDataDoubleJackAllMomenta;
typedef figureDataAllMomentaBase<figureDataBlockDoubleJack> figureDataBlockDoubleJackAllMomenta;

template<typename DistributionType>
struct _figureDataAllMomentaTypeSelector{ };
template<>
struct _figureDataAllMomentaTypeSelector<rawDataDistributionD>{ typedef figureDataAllMomenta type; };
template<>
struct _figureDataAllMomentaTypeSelector<jackknifeDistributionD>{ typedef figureDataJackAllMomenta type; };
template<>
struct _figureDataAllMomentaTypeSelector<doubleJackknifeDistributionD>{ typedef figureDataDoubleJackAllMomenta type; };
template<>
struct _figureDataAllMomentaTypeSelector<blockDoubleJackknifeDistributionD>{ typedef figureDataBlockDoubleJackAllMomenta type; };
template<>
struct _figureDataAllMomentaTypeSelector<bootstrapDistributionD>{ typedef figureDataAllMomentaBase<figureDataBoot>  type; };
template<>
struct _figureDataAllMomentaTypeSelector<bootJackknifeDistributionD>{ typedef figureDataAllMomentaBase<figureDataBootJack>  type; };



template<typename DistributionType>
using figureDataAllMomentaSelect = typename _figureDataAllMomentaTypeSelector<DistributionType>::type;


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

void zeroUnmeasuredSourceTimeslices(figureDataAllMomenta &data, const char fig, const int tstep_pipi){
  for(figureDataAllMomenta::iterator it = data.begin(fig); it != data.end(fig); it++){
    figureData &f = it->second;
    int Lt = f.getLt();
    for(int tsrc=0;tsrc<Lt;tsrc++)
      if(tsrc % tstep_pipi != 0)	  
	for(int tsep=0;tsep<Lt;tsep++) zeroit(f(tsrc,tsep));

  }    
}


SARLAC_END_NAMESPACE

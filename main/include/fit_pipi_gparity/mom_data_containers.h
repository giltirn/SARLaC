#ifndef PIPI_MOM_DATA_CONTAINERS_H
#define PIPI_MOM_DATA_CONTAINERS_H

#include <boost/serialization/map.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/binary_iarchive.hpp>


typedef std::array<int,3> threeMomentum;
typedef std::pair<threeMomentum, threeMomentum> sinkSourceMomenta;

std::ostream & operator<<(std::ostream &os, const threeMomentum &mom){
  os << "(" << mom[0] << ", " << mom[1] << ", " << mom[2] << ")";
  return os;
}
std::ostream & operator<<(std::ostream &os, const sinkSourceMomenta &mom){
  os << "Snk:" << mom.first << " Src:" << mom.second;
  return os;
}


inline threeMomentum operator-(const threeMomentum &p){
  return threeMomentum({-p[0],-p[1],-p[2]});
}
inline std::string momStr(const threeMomentum &p){
  std::ostringstream os;
  for(int i=0;i<3;i++)
    os << (p[i] < 0 ? "_" : "") << abs(p[i]);
  return os.str();
}
inline sinkSourceMomenta momComb(const int snkx, const int snky, const int snkz,
				 const int srcx, const int srcy, const int srcz){
  return sinkSourceMomenta( {snkx,snky,snkz}, {srcx,srcy,srcz} );
}
inline sinkSourceMomenta momComb(const threeMomentum &snk, const threeMomentum &src){
  return sinkSourceMomenta(snk,src);
}


template<typename _ContainerType>
class figureDataAllMomentaBase{
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
    return const_cast<MapType const*>( const_cast<figureDataAllMomentaBase<_ContainerType>* >(this)->getMap(fig) );
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

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & C & D & R & V & Lt & Nsample;
  }
public:
  figureDataAllMomentaBase(const int _Lt, const int _Nsample): Lt(_Lt), Nsample(_Nsample){}
  figureDataAllMomentaBase(){}

  void setup(const int _Lt, const int _Nsample){
    Lt = _Lt; Nsample = _Nsample;
  }
  
  const ContainerType &operator()(const char fig, const sinkSourceMomenta &mom) const{
    figureDataAllMomentaBase<ContainerType> *t = const_cast<figureDataAllMomentaBase<ContainerType> *>(this);
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
template<typename C>
inline void write(HDF5writer &writer, const figureDataAllMomentaBase<C> &d, const std::string &tag){ d.write(writer,tag); }
template<typename C>
inline void read(HDF5reader &reader, figureDataAllMomentaBase<C> &d, const std::string &tag){ d.read(reader,tag); }
#endif


typedef figureDataAllMomentaBase<figureData> figureDataAllMomenta;
typedef figureDataAllMomentaBase<figureDataDoubleJack> figureDataDoubleJackAllMomenta;




template<typename _ContainerType>
class bubbleDataAllMomentaBase{
public:
  typedef _ContainerType ContainerType;
  typedef typename std::map<threeMomentum, ContainerType>::iterator iterator;  
  typedef typename std::map<threeMomentum, ContainerType>::const_iterator const_iterator;  
private:
  std::map<threeMomentum, ContainerType> B;
  int Lt;
  int Nsample;
  
  ContainerType & get(const threeMomentum &mom, bool lock){
    typename std::map<threeMomentum, ContainerType>::iterator it = B.find(mom);
    if(it == B.end()){
      if(lock) error_exit(std::cout << "bubbleDataAllMomenta::get Could not find requested momentum\n");

      it = B.insert(std::make_pair(mom, ContainerType())).first;
      it->second.setup(Lt,Nsample);
    }
    return it->second;
  }

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & B & Lt & Nsample;
  }
public:
  bubbleDataAllMomentaBase(const int _Lt, const int _Nsample): Lt(_Lt), Nsample(_Nsample){}
  bubbleDataAllMomentaBase(){}

  void setup(const int _Lt, const int _Nsample){
    Lt = _Lt; Nsample = _Nsample;
  }
  
  const ContainerType &operator()(const threeMomentum &mom) const{
    bubbleDataAllMomentaBase<ContainerType> *t = const_cast<bubbleDataAllMomentaBase<ContainerType> *>(this);
    return const_cast<const ContainerType &>( t->get(mom,true) );
  }
  ContainerType &operator()(const threeMomentum &mom){
    return this->get(mom,false);
  }
  inline int getLt() const{ return Lt; }
  inline int getNsample() const{ return Nsample; }

  inline const_iterator begin() const{ return B.begin(); }
  inline const_iterator end() const{ return B.end(); }

  inline iterator begin(){ return B.begin(); }
  inline iterator end(){ return B.end(); }

  GENERATE_HDF5_SERIALIZE_METHOD((B)(Lt)(Nsample));
};
#ifdef HAVE_HDF5
template<typename C>
inline void write(HDF5writer &writer, const bubbleDataAllMomentaBase<C> &d, const std::string &tag){ d.write(writer,tag); }
template<typename C>
inline void read(HDF5reader &reader, bubbleDataAllMomentaBase<C> &d, const std::string &tag){ d.read(reader,tag); }
#endif

typedef bubbleDataAllMomentaBase<bubbleData> bubbleDataAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataDoubleJack> bubbleDataDoubleJackAllMomenta;

template<typename archiver>
struct archiveStream{};

template<>
struct archiveStream<boost::archive::binary_oarchive>{
  std::ofstream ofs;
  archiveStream(const std::string &file): ofs(file.c_str(),std::ios::binary){}
};
template<>
struct archiveStream<boost::archive::text_oarchive>{
  std::ofstream ofs;
  archiveStream(const std::string &file): ofs(file.c_str()){}
};
template<>
struct archiveStream<boost::archive::binary_iarchive>{
  std::ifstream ifs;
  archiveStream(const std::string &file): ifs(file.c_str(),std::ios::binary){}
};
template<>
struct archiveStream<boost::archive::text_iarchive>{
  std::ifstream ifs;
  archiveStream(const std::string &file): ifs(file.c_str()){}
};

template<typename archiver>
void saveCheckpoint(const figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const std::string &file){
  (std::cout << "Saving data checkpoint\n").flush(); boost::timer::auto_cpu_timer t("Report: Saved data checkpoint in %w s\n");
  archiveStream<archiver> st(file);
  archiver oa(st.ofs);    
  oa << raw_data;
  oa << raw_bubble_data;
}
template<typename archiver>
void loadCheckpoint(figureDataAllMomenta &raw_data, bubbleDataAllMomenta &raw_bubble_data, const std::string &file){
  (std::cout << "Loading data checkpoint\n").flush(); boost::timer::auto_cpu_timer t("Report: Loaded data checkpoint in %w s\n");
  archiveStream<archiver> st(file);
  archiver ia(st.ifs);
  ia >> raw_data;
  ia >> raw_bubble_data;
}



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

#endif

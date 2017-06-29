#ifndef PIPI_MOM_DATA_CONTAINERS_H
#define PIPI_MOM_DATA_CONTAINERS_H

#include <boost/serialization/map.hpp>
#include <boost/serialization/complex.hpp>
#include <boost/serialization/array.hpp>

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

    
};


typedef figureDataAllMomentaBase<figureData> figureDataAllMomenta;
typedef figureDataAllMomentaBase<figureDataDoubleJack> figureDataDoubleJackAllMomenta;




template<typename _ContainerType>
class bubbleDataAllMomentaBase{
public:
  typedef _ContainerType ContainerType;
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
};


typedef bubbleDataAllMomentaBase<bubbleData> bubbleDataAllMomenta;
typedef bubbleDataAllMomentaBase<bubbleDataDoubleJack> bubbleDataDoubleJackAllMomenta;



#endif

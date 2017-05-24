#ifndef _SUPERJACKKNIFE_H_
#define _SUPERJACKKNIFE_H_

#include <distribution.h>
#include <generic_ET.h>

class superJackknifeLayout{
  std::vector<std::string > ens_tags;
  std::vector<int> ens_sizes;
  std::vector< std::pair<int,int> > ens_sample_map;
  int total_size;
public:
  superJackknifeLayout(): total_size(0){}
  
  //returns -1 for non-existent
  inline int ensIdx(const std::string &e) const{
    for(int i=0;i<ens_tags.size();i++) if(e == ens_tags[i]) return i;
    return -1;
  }

  void addEnsemble(const std::string &tag, const int size){
    int ens_idx = ens_tags.size();
    ens_tags.push_back(tag);
    ens_sizes.push_back(size);
    for(int i=0;i<size;i++) ens_sample_map.push_back(std::pair<int,int>(ens_idx,i) );
    total_size += size;
  }

  inline int nSamplesTotal() const{ return total_size; }

  inline int nEnsembles() const{ return ens_tags.size(); }

  inline int nSamplesEns(const int ens_idx) const{ return ens_sizes[ens_idx]; }
  inline int nSamplesEns(const std::string &ens_tag) const{ return nSamplesEns(ensIdx(ens_tag)); }
  
  inline const std::pair<int,int> & sampleMap(int gsample) const{
    return ens_sample_map[gsample];
  }
};

template<typename _DataType>
class superJackknifeDistribution{
public:
  typedef _DataType DataType;
protected:
  std::vector< jackknifeDistribution<DataType> > ens_jacks;
  DataType cen;
  superJackknifeLayout const* layout;

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version){
    ar & ens_jacks;
    ar & cen;
  }

public:
  superJackknifeDistribution(const superJackknifeLayout &_layout, const DataType &_central): layout(&_layout), ens_jacks(_layout.nEnsembles()){
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i), _central ); 
    cen = _central;
  }

  superJackknifeDistribution(const superJackknifeLayout &_layout, const int first_ens_idx, const jackknifeDistribution<DataType> &first): superJackknifeDistribution(_layout, first.mean()) {
    assert(first.size() == layout->nSamplesEns(first_ens_idx));
    ens_jacks[first_ens_idx] = first;
  }
  
  superJackknifeDistribution(const superJackknifeDistribution &r): layout(r.layout), ens_jacks(r.ens_jacks), cen(r.cen){}
  superJackknifeDistribution(superJackknifeDistribution&& o) noexcept : layout(o.layout), ens_jacks(std::move(o.ens_jacks)), cen(o.cen){ }

  typedef superJackknifeDistribution<DataType> ET_tag;
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,superJackknifeDistribution<DataType> >::value, int>::type = 0>
  superJackknifeDistribution(U&& expr){
    layout = expr.common_properties();
    ens_jacks.resize(layout->nEnsembles());
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i) ); 
    
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
    cen = expr[-1];
  }
  
  superJackknifeDistribution & operator=(const superJackknifeDistribution &r){ layout = r.layout; ens_jacks = r.ens_jacks; cen = r.cen; return *this; }
  
  int size() const{ return layout->nSamplesTotal(); }
  
  const DataType & sample(const int idx) const{
    return ens_jacks[ layout->sampleMap(idx).first ].sample( layout->sampleMap(idx).second );
  }
  DataType & sample(const int idx){
    return ens_jacks[ layout->sampleMap(idx).first ].sample( layout->sampleMap(idx).second );
  }

  const DataType & best() const{ return cen; }
  DataType &best(){ return cen; }

  //This is the naive mean of the superjackknife samples. Use best for the propagated central value
  DataType mean(){
    DataType out(0.);
    int total_size = 0;    
    for(int i=0;i<ens_jacks.size();i++){
      int sz = ens_jacks[i].size();
      total_size += sz;
      out = out + ens_jacks[i]*double(sz);
    }
    out = out / double(total_size);
    return out;      
  }
  
  DataType standardError() const{
    DataType v = 0;
    for(int i=0;i<ens_jacks.size();i++){
      DataType vi = ens_jacks[i].standardError();
      v = v + vi*vi;
    }
    return sqrt(v);
  }

  void setEnsembleJackknife(const int ens_idx, const jackknifeDistribution<DataType> &j){
    assert(j.best() == cen);
    assert(j.size() == layout->nSamplesEns(ens_idx));
    ens_jacks[ens_idx] = j;    
  }
  inline void setEnsembleJackknife(const std::string & ens_tag, const jackknifeDistribution<DataType> &j){
    setEnsembleJackknife(layout->ensIdx(ens_tag),j);
  }

  inline const jackknifeDistribution<DataType> & getEnsembleJackknife(const int ens_idx) const{
    return ens_jacks[ens_idx];
  }
  inline const jackknifeDistribution<DataType> & getEnsembleJackknife(const std::string &ens_tag) const{
    return ens_jacks[layout->ensIdx(ens_tag)];
  }
  
  const superJackknifeLayout & getLayout() const{ return *layout; }
};

template<typename A>
struct getElem<superJackknifeDistribution<A> >{
  static inline auto elem(const superJackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return i == -1 ? v.best() : v.sample(i); }
  static inline auto elem(superJackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return i == -1 ? v.best() : v.sample(i); }
  static inline superJackknifeLayout const* common_properties(const superJackknifeDistribution<A> &v){ return &v.getLayout(); }
};


#endif

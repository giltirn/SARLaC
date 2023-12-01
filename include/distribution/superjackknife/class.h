#ifndef _SUPERJACKKNIFE_CLASS_H_
#define _SUPERJACKKNIFE_CLASS_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<distribution/jackknife.h>

#include<distribution/superjackknife/layout.h>

SARLAC_START_NAMESPACE


template<typename _DataType>
class superJackknifeDistribution{
public:
  typedef _DataType DataType;

  template<typename T>
  using rebase = superJackknifeDistribution<T>;
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
  inline static int checkEnsExistsPassthrough(const std::string &tag, const superJackknifeLayout &layout){
    int e = layout.ensIdx(tag);
    if(e==-1) error_exit(std::cout << "superJackknifeDistribution::checkEnsExistsPassthrough fail " << tag << std::endl);
    return e;
  }
  
public:
  superJackknifeDistribution(): layout(NULL){}
  
  superJackknifeDistribution(const superJackknifeLayout &_layout, const DataType &_central): layout(&_layout), ens_jacks(_layout.nEnsembles()){
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i), _central ); 
    cen = _central;
  }

  //Lambda-type initializer. Expect operator()(const int sample) that accepts sample=-1 for the central value
  template<typename Initializer>
  superJackknifeDistribution(const superJackknifeLayout &_layout, const Initializer &init): layout(&_layout), ens_jacks(_layout.nEnsembles()){
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i) );
    cen = init(-1);
    for(int i=0;i<layout->nSamplesTotal();i++) this->sample(i) = init(i);
  }			     

  superJackknifeDistribution(const superJackknifeLayout &_layout, const int first_ens_idx, const jackknifeDistribution<DataType> &first): superJackknifeDistribution(_layout, first.mean()) {
    assert(first.size() == layout->nSamplesEns(first_ens_idx));
    ens_jacks[first_ens_idx] = first;
  }
  superJackknifeDistribution(const superJackknifeLayout &_layout, const std::string &first_ens_tag, const jackknifeDistribution<DataType> &first): superJackknifeDistribution(_layout,checkEnsExistsPassthrough(first_ens_tag,_layout),first){}
  
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

  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,superJackknifeDistribution<DataType> >::value, int>::type = 0>
  superJackknifeDistribution<DataType> & operator=(U&& expr){
    layout = expr.common_properties();
    ens_jacks.resize(layout->nEnsembles());
    for(int i=0;i<ens_jacks.size();i++) ens_jacks[i].resize( layout->nSamplesEns(i) ); 
    
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
    cen = expr[-1];
    return *this;
  }
  
  inline superJackknifeDistribution & operator=(const superJackknifeDistribution &r){
    layout = r.layout; ens_jacks = r.ens_jacks; cen = r.cen; return *this;
  }
  inline superJackknifeDistribution & operator=(superJackknifeDistribution &&r){
    layout = r.layout; ens_jacks = std::move(r.ens_jacks); cen = r.cen; return *this;
  }
  
  inline int size() const{ return layout->nSamplesTotal(); }

  inline void resize(const int sz){ if(sz != size()) error_exit(std::cout << "Cannot resize a superJackknifeDistribution without changing the layout!\n"); }
  
  inline const DataType & sample(const int idx) const{
    return ens_jacks[ layout->sampleMap(idx).first ].sample( layout->sampleMap(idx).second );
  }
  inline DataType & sample(const int idx){
    return ens_jacks[ layout->sampleMap(idx).first ].sample( layout->sampleMap(idx).second );
  }

  inline const DataType & best() const{ return cen; }
  inline DataType &best(){ return cen; }

  //Versions of sample that map sample -1 to the central value, useful for looping
  inline const DataType & osample(const int idx) const{
    return idx == -1 ? this->best() : this->sample(idx);
  }
  inline DataType & osample(const int idx){
    return idx == -1 ? this->best() : this->sample(idx);
  }

  bool operator==(const superJackknifeDistribution &r) const{
    if(*layout != *r.layout) return false;
    if(cen != r.cen) return false;
    for(int i=0;i<ens_jacks.size();i++) if(ens_jacks[i] != r.ens_jacks[i]) return false;
    return true;
  }
  inline bool operator!=(const superJackknifeDistribution &r) const{ return !(*this == r); }
  
  //This is the naive mean of the superjackknife samples. Use best for the propagated central value
  DataType mean(){
    DataType out = best(); zeroit(out);
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
    DataType v = best(); zeroit(v);
    for(int i=0;i<ens_jacks.size();i++){
      DataType vi = ens_jacks[i].standardError();
      v = v + vi*vi;
    }
    return sqrt(v);
  }

  void setEnsembleJackknife(const int ens_idx, const jackknifeDistribution<DataType> &j){
    //assert(j.best() == cen); //this is actually not true for most jackknife
    if(j.size() != layout->nSamplesEns(ens_idx)) error_exit(std::cout << "superJackknifeDistribution::setEnsembleJackknife(int, jack) input jackknife has size " << j.size() << ", expected " << layout->nSamplesEns(ens_idx) <<std::endl);
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

  //Set the layout to that provided. New layout must contain all the original ensembles
  void setLayout(const superJackknifeLayout &to){
    if(layout == NULL){ //no layout set previously (i.e. default constructor)
      ens_jacks.resize(to.nEnsembles());
      for(int e=0;e<to.nEnsembles();e++) ens_jacks[e].resize(to.nSamplesEns(e));
      layout = &to;
      return;
    }
    if(&to == layout) return; //is identical
    if(to == *layout){ //no reorg needed
      layout = &to;
      return;
    }
    
    if(to.nEnsembles() < layout->nEnsembles()) error_exit(std::cout << "superJackknifeDistribution::setLayout new layout must contain an equal or larger number of ensembles\n");   

    std::vector<jackknifeDistribution<DataType> > result(to.nEnsembles());

    std::vector<bool> fill(to.nEnsembles(),true);
    
    for(int e=0;e<layout->nEnsembles();e++){
      int new_idx = to.ensIdx(layout->ensTag(e));
      if(new_idx == -1) error_exit(std::cout << "superJackknifeDistribution::setLayout ensemble '" << layout->ensTag(e) << "' does not exist in new layout\n");
      result[new_idx] = std::move(ens_jacks[e]);      
      fill[new_idx] = false;      
    }

    for(int e=0;e<to.nEnsembles();e++)
      if(fill[e])
	result[e] = jackknifeDistribution<DataType>(to.nSamplesEns(e),cen);

    ens_jacks = std::move(result);
    layout = &to;
  }

  //Allow generic initialization and resize
  typedef superJackknifeLayout initType;
  
  inline superJackknifeLayout getInitializer() const{ return this->getLayout(); }

  void resize(const superJackknifeLayout &layout){ this->setLayout(layout); }
  
  template<typename U=DataType, typename std::enable_if< is_std_complex<U>::value, int >::type = 0>
  superJackknifeDistribution<typename U::value_type> real() const{
    superJackknifeDistribution<typename U::value_type> out(*layout, this->best().real());
    for(int i=0;i<size();i++) out.sample(i) = this->sample(i).real();
    return out;
  }

  inline void zero(){
    zeroit(this->best());
    for(int i=0;i<this->size();i++) zeroit(this->sample(i));
  }

#ifdef HAVE_HDF5
  void write(HDF5writer &writer, const std::string &tag) const{
    writer.enter(tag);
    SARLaC::write(writer,ens_jacks,"ens_jacks",false); //don't flatten the internal vector
    SARLaC::write(writer,cen,"cen");
    SARLaC::write(writer,*layout,"layout");
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);
    SARLaC::read(reader,ens_jacks,"ens_jacks",false);
    SARLaC::read(reader,cen, "cen");

    static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));
    superJackknifeLayout* _layout = layouts.back().get();
    SARLaC::read(reader,*_layout,"layout");
    layout = _layout;
    reader.leave();
  }
#endif
  
};


template<typename A>
struct getElem<superJackknifeDistribution<A> >{
  static inline auto elem(const superJackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return i == -1 ? v.best() : v.sample(i); }
  static inline auto elem(superJackknifeDistribution<A> &v, const int i)->decltype(v.sample(i)){ return i == -1 ? v.best() : v.sample(i); }
  static inline superJackknifeLayout const* common_properties(const superJackknifeDistribution<A> &v){ return &v.getLayout(); }
};

template<typename T>
std::ostream & operator<<(std::ostream &os, const superJackknifeDistribution<T> &d){
  assert(distributionPrint<superJackknifeDistribution<T> >::printer() != NULL); distributionPrint<superJackknifeDistribution<T> >::printer()->print(os, d);
  return os;
}

SARLAC_END_NAMESPACE
#endif

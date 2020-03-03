#ifndef _SUPERMULTI_CLASS_H_
#define _SUPERMULTI_CLASS_H_

#include<config.h>
#include<ET/generic_ET.h>
#include<containers/general_container.h>
#include<distribution/jackknife.h>
#include<distribution/bootstrap.h>

#include<distribution/supermulti/layout.h>

CPSFIT_START_NAMESPACE


template<typename _DataType, template<typename> class _VectorType = basic_vector>
class superMultiDistribution: public distribution<_DataType, _VectorType>{
public:
  typedef _DataType DataType;
  template<typename T>
  using rebase = superMultiDistribution<T,_VectorType>;
protected:
  DataType cen;
  superMultiLayout const* layout;

  inline static int checkEnsExistsPassthrough(const std::string &tag, const superMultiLayout &layout){
    int e = layout.ensIdx(tag);
    if(e==-1) error_exit(std::cout << "superMultiDistribution::checkEnsExistsPassthrough fail " << tag << std::endl);
    return e;
  }
  
public:
  superMultiDistribution(): layout(NULL){}
  
  //Create with a given layout
  superMultiDistribution(superMultiLayout const *_layout): layout(_layout), distribution<_DataType, _VectorType>(_layout->nSamplesTotal()){  }
  superMultiDistribution(const superMultiLayout &_layout): superMultiDistribution(&_layout){  }

  //Create with a given layout and central value (copied to all samples so zero error)
  superMultiDistribution(superMultiLayout const *_layout, const DataType &_central): layout(_layout), cen(_central), 
										     distribution<_DataType, _VectorType>(_layout->nSamplesTotal(), _central){}
  superMultiDistribution(const superMultiLayout &_layout, const DataType &_central): superMultiDistribution(&_layout,_central){}


  //Lambda-type initializer. Expect operator()(const int sample) that accepts sample=-1 for the central value. Use osample accessor for superMultiDistribution samples internally to propagate this
  template<typename Initializer>
  superMultiDistribution(superMultiLayout const *_layout, const Initializer &init): layout(_layout), 
										    distribution<_DataType, _VectorType>(_layout->nSamplesTotal()){
    cen = init(-1);
    for(int i=0;i<this->size();i++) this->sample(i) = init(i);
  }			     
  template<typename Initializer>
  superMultiDistribution(const superMultiLayout &_layout, const Initializer &init): superMultiDistribution(&_layout, init){}


  //Create a multi-distribution from a jackknife placed on ensemble first_ens_idx
  superMultiDistribution(superMultiLayout const *_layout, const int first_ens_idx, const jackknifeDistribution<DataType> &first): superMultiDistribution(_layout, first.mean()) {
    setEnsembleDistribution(first_ens_idx, first);
  }
  superMultiDistribution(const superMultiLayout &_layout, const int first_ens_idx, const jackknifeDistribution<DataType> &first): superMultiDistribution(&_layout, first_ens_idx, first){}

  //Create from a jackknife placed on ensemble with given tag
  superMultiDistribution(superMultiLayout const *_layout, const std::string &first_ens_tag, const jackknifeDistribution<DataType> &first): superMultiDistribution(_layout,checkEnsExistsPassthrough(first_ens_tag,*_layout),first){}
  superMultiDistribution(const superMultiLayout &_layout, const std::string &first_ens_tag, const jackknifeDistribution<DataType> &first): superMultiDistribution(&_layout, first_ens_tag, first){}



  //Create a multi-distribution from a bootstrap placed on ensemble first_ens_idx
  superMultiDistribution(superMultiLayout const *_layout, const int first_ens_idx, const bootstrapDistribution<DataType> &first): superMultiDistribution(_layout, first.best()) {
    setEnsembleDistribution(first_ens_idx, first);
  }
  superMultiDistribution(const superMultiLayout &_layout, const int first_ens_idx, const bootstrapDistribution<DataType> &first): superMultiDistribution(&_layout, first_ens_idx, first){}


  //Create a multi-distribution from a bootstrap placed on ensemble with given tag
  superMultiDistribution(superMultiLayout const *_layout, const std::string &first_ens_tag, const bootstrapDistribution<DataType> &first): superMultiDistribution(_layout,checkEnsExistsPassthrough(first_ens_tag,*_layout),first){}
  superMultiDistribution(const superMultiLayout &_layout, const std::string &first_ens_tag, const bootstrapDistribution<DataType> &first): superMultiDistribution(&_layout, first_ens_tag, first){}



  //Create a multi-distribution from a wrapped jackknife/bootstrap placed on ensemble first_ens_idx
  superMultiDistribution(superMultiLayout const *_layout, const int first_ens_idx, const generalContainer &first): superMultiDistribution(_layout) {
    if(first.is<jackknifeDistribution<DataType> >()) cen = first.value<jackknifeDistribution<DataType> >().best();
    else if(first.is<bootstrapDistribution<DataType> >()) cen = first.value<bootstrapDistribution<DataType> >().best();
    else assert(0);
    setEnsembleDistribution(first_ens_idx, first);
  }
  superMultiDistribution(const superMultiLayout &_layout, const int first_ens_idx, const generalContainer &first): superMultiDistribution(&_layout, first_ens_idx, first){}

  //Create a multi-distribution from a wrapped jackknife/bootstrap placed on ensemble with given tag
  superMultiDistribution(superMultiLayout const *_layout, const std::string &first_ens_tag, const generalContainer &first): superMultiDistribution(_layout,checkEnsExistsPassthrough(first_ens_tag,*_layout),first){}
  superMultiDistribution(const superMultiLayout &_layout, const std::string &first_ens_tag, const generalContainer &first):  superMultiDistribution(&_layout, first_ens_tag, first){}


  superMultiDistribution(const superMultiDistribution &r): layout(r.layout), cen(r.cen), distribution<_DataType, _VectorType>(r){}
  superMultiDistribution(superMultiDistribution&& o) noexcept : layout(o.layout), cen(o.cen), distribution<_DataType, _VectorType>(std::move(o)){ }

  typedef superMultiDistribution<DataType> ET_tag;
  
  //Initialize distribution with a lambda of signature DataType [](const int s) . The lambda must return the central value for index s=-1
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,superMultiDistribution<DataType> >::value, int>::type = 0>
  superMultiDistribution(U&& expr): layout(NULL){
    this->setLayout(*expr.common_properties());
    cen = expr[-1];
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
  }

  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,superMultiDistribution<DataType> >::value, int>::type = 0>
  superMultiDistribution<DataType> & operator=(U&& expr){
    this->setLayout(*expr.common_properties());
    cen = expr[-1];
    for(int i=0;i<this->size();i++) this->sample(i) = expr[i];
    return *this;
  }

  inline superMultiDistribution & operator=(const superMultiDistribution &r){
    layout = r.layout; cen = r.cen; this->distribution<_DataType, _VectorType>::operator=(r);  return *this;
  }
  inline superMultiDistribution & operator=(superMultiDistribution &&r){
    layout = r.layout; cen = r.cen; this->distribution<_DataType, _VectorType>::operator=(std::move(r));  return *this;
  }
  
  inline int size() const{ return layout->nSamplesTotal(); }

  inline void resize(const int sz){ if(sz != size()) error_exit(std::cout << "Cannot resize a superMultiDistribution without changing the layout!\n"); }
  void resize(const int sz, const DataType &init_val){ if(sz != size()) error_exit(std::cout << "Cannot resize a superMultiDistribution without changing the layout!\n"); }

  inline const DataType & best() const{ return cen; }
  inline DataType & best(){ return cen; }
  
  inline DataType & osample(const size_t i){ return i==-1 ? this->best() : this->sample(i); }
  inline const DataType & osample(const size_t i) const{ return i==-1 ? this->best() : this->sample(i); }

  bool operator==(const superMultiDistribution &r) const{
    if(*layout != *r.layout) return false;
    if(cen != r.cen) return false;
    if(!(this->distribution<_DataType, _VectorType>::operator==(r))) return false;
    return true;

  }
  inline bool operator!=(const superMultiDistribution &r) const{ return !(*this == r); }
   
  DataType standardError() const{
    DataType v = best(); zeroit(v);
    size_t off = 0;

    for(int i=0;i<layout->nEnsembles();i++){
      DataType vi;
      size_t ens_sz = layout->nSamplesEns(i);
      if(layout->ensType(i) == MultiType::Jackknife){
	jackknifeDistribution<DataType> j(ens_sz);
	for(int s=0;s<ens_sz;s++) j.sample(s) = this->sample(s + off);
	vi = j.standardError();
      }else if(layout->ensType(i) == MultiType::Bootstrap){
	bootstrapInitType init(ens_sz);
	bootstrapDistribution<DataType> j(init);
	j.best() = cen;
	for(int s=0;s<ens_sz;s++) j.sample(s) = this->sample(s + off);
	vi = j.standardError();
      }else assert(0);
      
      v = v + vi*vi;
      off += ens_sz;
    }
    return sqrt(v);
  }

  static DataType covariance(const superMultiDistribution<_DataType,_VectorType> &A, const superMultiDistribution<_DataType,_VectorType> &B){
    assert(A.getLayout() == B.getLayout());
    
    DataType v = A.sample(0); zeroit(v);
    size_t off = 0;

    const superMultiLayout &layout = A.getLayout();    

    for(int i=0;i<layout.nEnsembles();i++){
      DataType vi;
      size_t ens_sz = layout.nSamplesEns(i);
      if(layout.ensType(i) == MultiType::Jackknife){
	jackknifeDistribution<_DataType> jA(ens_sz, [&](const int s){ return A.sample(s+off); });
	jackknifeDistribution<_DataType> jB(ens_sz, [&](const int s){ return B.sample(s+off); });
	vi = jackknifeDistribution<DataType>::covariance(jA, jB);
      }else if(layout.ensType(i) == MultiType::Bootstrap){
	bootstrapInitType init(ens_sz);
	bootstrapDistribution<_DataType> jA(init, [&](const int s){ return A.sample(s+off); });
	bootstrapDistribution<_DataType> jB(init, [&](const int s){ return B.sample(s+off); });
	jA.best() = A.best();
	jB.best() = B.best();
	vi = bootstrapDistribution<DataType>::covariance(jA, jB);
      }else assert(0);
      
      v = v + vi;
      off += ens_sz;
    }
    return v;
  }

  void setEnsembleDistribution(const int ens_idx, const jackknifeDistribution<DataType> &j){
    if(layout->ensType(ens_idx) != MultiType::Jackknife) error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, jack) wrong distribution type" <<std::endl);
    if(j.size() != layout->nSamplesEns(ens_idx)) error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, jack) input jackknife has size " << j.size() << ", expected " << layout->nSamplesEns(ens_idx) <<std::endl);
    size_t off = layout->offset(ens_idx);
    for(size_t i=0;i<j.size();i++) this->sample(i+off) = j.sample(i);
  }
  inline void setEnsembleDistribution(const std::string & ens_tag, const jackknifeDistribution<DataType> &j){
    setEnsembleDistribution(layout->ensIdx(ens_tag),j);
  }

  void setEnsembleDistribution(const int ens_idx, const bootstrapDistribution<DataType> &j){
    if(layout->ensType(ens_idx) != MultiType::Bootstrap) error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, boot) wrong distribution type" <<std::endl);
    if(j.size() != layout->nSamplesEns(ens_idx)) error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, boot) input bootstrap has size " << j.size() << ", expected " << layout->nSamplesEns(ens_idx) <<std::endl);
    size_t off = layout->offset(ens_idx);
    for(size_t i=0;i<j.size();i++) this->sample(i+off) = j.sample(i);
  }
  inline void setEnsembleDistribution(const std::string & ens_tag, const bootstrapDistribution<DataType> &j){
    setEnsembleDistribution(layout->ensIdx(ens_tag),j);
  }


  void setEnsembleDistribution(const int ens_idx, const generalContainer &j){
    std::vector<DataType> const* from;

    if(j.is<jackknifeDistribution<DataType> >()){
      const jackknifeDistribution<DataType> &jc = j.value<jackknifeDistribution<DataType> >();
      if(layout->ensType(ens_idx) != MultiType::Jackknife) 
	error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, boot) wrong distribution type" <<std::endl);
      if(jc.size() != layout->nSamplesEns(ens_idx)) 
	error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, boot) input bootstrap has size " << jc.size() << 
		   ", expected " << layout->nSamplesEns(ens_idx) <<std::endl);
      from = &jc.sampleVector();
    }else if(j.is<bootstrapDistribution<DataType> >()){
      const bootstrapDistribution<DataType> &jc = j.value<bootstrapDistribution<DataType> >();
      if(layout->ensType(ens_idx) != MultiType::Bootstrap) 
	error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, boot) wrong distribution type" <<std::endl);
      if(jc.size() != layout->nSamplesEns(ens_idx)) 
	error_exit(std::cout << "superMultiDistribution::setEnsembleDistribution(int, boot) input bootstrap has size " << jc.size() << 
		   ", expected " << layout->nSamplesEns(ens_idx) <<std::endl);
      from = &jc.sampleVector();
    }else assert(0);

    size_t off = layout->offset(ens_idx);
    for(size_t i=0;i<from->size();i++) this->sample(i+off) = from->operator[](i);
  }
  inline void setEnsembleDistribution(const std::string & ens_tag, const generalContainer &j){
    setEnsembleDistribution(layout->ensIdx(ens_tag),j);
  }


  inline generalContainer getEnsembleDistribution(const int ens_idx) const{
    size_t off = layout->offset(ens_idx);
    size_t ens_sz = layout->nSamplesEns(ens_idx);
    MultiType type = layout->ensType(ens_idx);
    if(type == MultiType::Jackknife){
      jackknifeDistribution<DataType> j(ens_sz);
      for(int s=0;s<ens_sz;s++) j.sample(s) = this->sample(s + off);
      return generalContainer(std::move(j));
    }else if(type == MultiType::Bootstrap){
      bootstrapInitType init(ens_sz);
      bootstrapDistribution<DataType> j(init);
      j.best() = cen;
      for(int s=0;s<ens_sz;s++) j.sample(s) = this->sample(s + off);
      return generalContainer(std::move(j));
    }else assert(0);
  }
  inline generalContainer getEnsembleDistribution(const std::string &ens_tag) const{
    return getEnsembleDistribution(layout->ensIdx(ens_tag));
  }

  const superMultiLayout & getLayout() const{ return *layout; }
  superMultiLayout const * getInitializer() const{ return layout; }
  
  //Set the layout to that provided. New layout must contain all the original ensembles
  void setLayout(const superMultiLayout &to){
    if(layout == NULL){ //no layout set previously (i.e. default constructor)
      this->distribution<_DataType,_VectorType>::resize(to.nSamplesTotal());
      layout = &to;
      return;
    }
    if(&to == layout) return; //is identical
    if(to == *layout){ //no reorg needed
      layout = &to;
      return;
    }
    
    if(to.nEnsembles() < layout->nEnsembles()) error_exit(std::cout << "superMultiDistribution::setLayout new layout must contain an equal or larger number of ensembles\n");   

    _VectorType<_DataType> result(to.nSamplesTotal(), cen);
    
    for(int e=0;e<layout->nEnsembles();e++){
      int new_idx = to.ensIdx(layout->ensTag(e));
      if(new_idx == -1) error_exit(std::cout << "superMultiDistribution::setLayout ensemble '" << layout->ensTag(e) << "' does not exist in new layout\n");
      
      size_t new_off = to.offset(new_idx);
      size_t cur_off = layout->offset(e);
      size_t ens_sz = layout->nSamplesEns(e);
      for(size_t s=0;s<ens_sz;s++) result[s + new_off] = this->sample(s + cur_off);
    } 
    this->sampleVector() = std::move(result);
    layout = &to;
  }
  
  inline void zero(){
    zeroit(cen);
    this->distribution<_DataType,_VectorType>::zero();
  }

#ifdef HAVE_HDF5
  void write(HDF5writer &writer, const std::string &tag) const{
    writer.enter(tag);
    CPSfit::write(writer,*layout,"layout");
    CPSfit::write(writer,cen,"cen");
    CPSfit::write(writer,this->_data,"data");
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);

    static std::vector<std::unique_ptr<superMultiLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superMultiLayout>(new superMultiLayout));
    superMultiLayout* _layout = layouts.back().get();
    CPSfit::read(reader,*_layout,"layout");
    layout = _layout;

    CPSfit::read(reader,cen, "cen");    
    CPSfit::read(reader,this->_data, "data");
    reader.leave();
  }
#endif
  
};


template<typename A, template<typename> class V>
struct getElem<superMultiDistribution<A,V> >{
  static inline auto elem(const superMultiDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return i == -1 ? v.best() : v.sample(i); }
  static inline auto elem(superMultiDistribution<A,V> &v, const int i)->decltype(v.sample(i)){ return i == -1 ? v.best() : v.sample(i); }
  static inline superMultiLayout const* common_properties(const superMultiDistribution<A,V> &v){ return &v.getLayout(); }
};

template<typename T, template<typename> class V>
std::ostream & operator<<(std::ostream &os, const superMultiDistribution<T,V> &d){
  auto const* printer = distributionPrint<superMultiDistribution<T,V> >::printer();
  assert(printer != NULL); printer->print(os, d);
  return os;
}

CPSFIT_END_NAMESPACE
#endif

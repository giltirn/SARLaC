#ifndef _SUPERJACKKNIFE_H_
#define _SUPERJACKKNIFE_H_

#include <distribution.h>
#include <generic_ET.h>

class superJackknifeLayout{
  std::vector<std::string > ens_tags;
  std::vector<int> ens_sizes;
  std::vector< std::pair<int,int> > ens_sample_map;
  int total_size;

  inline void clear(){
    ens_tags.resize(0);
    ens_sizes.resize(0);
    ens_sample_map.resize(0);
    total_size = 0;
  }
    
public:
  superJackknifeLayout(): total_size(0){}
  
  //returns -1 for non-existent
  inline int ensIdx(const std::string &e) const{
    for(int i=0;i<ens_tags.size();i++) if(e == ens_tags[i]) return i;
    return -1;
  }
  inline const std::string &ensTag(const int i) const{ return ens_tags[i]; }

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

  bool operator==(const superJackknifeLayout &r) const{
    if(&r == this) return true;
    if(total_size != r.total_size) return false;
    for(int e=0;e<ens_sizes.size();e++) if(ens_sizes[e] != r.ens_sizes[e]) return false;
    for(int e=0;e<ens_sample_map.size();e++) if(ens_sample_map[e] != r.ens_sample_map[e]) return false;
    for(int e=0;e<ens_tags.size();e++) if(ens_tags[e] != r.ens_tags[e]) return false;
    return true;
  }
  inline bool operator!=(const superJackknifeLayout &r) const{ return !(*this == r); }
  
#ifdef HAVE_HDF5
  void write(HDF5writer &writer, const std::string &tag) const{
    writer.enter(tag);
    ::write(writer,ens_tags,"ens_tags");
    ::write(writer,ens_sizes,"ens_sizes");
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);
    std::vector<std::string > tags;
    std::vector<int> sizes;
    
    ::read(reader,tags,"ens_tags");
    ::read(reader,sizes,"ens_sizes");

    clear();
    for(int i=0;i<tags.size();i++) addEnsemble(tags[i],sizes[i]);

    reader.leave();
  }
#endif
};
GENERATE_HDF5_SERIALIZE_FUNC(superJackknifeLayout);


superJackknifeLayout combine(const superJackknifeLayout &l, const superJackknifeLayout &r){
  superJackknifeLayout out(l);
  for(int i=0;i<r.nEnsembles();i++){
    const std::string &tag = r.ensTag(i);
    int e = out.ensIdx(tag);
    if(e == -1){ //doesn't currently exist
      out.addEnsemble(tag, r.nSamplesEns(i));
    }else if(out.nSamplesEns(e) != r.nSamplesEns(i)){
      error_exit(std::cout << "combine(const superJackknifeLayout &, const superJackknifeLayout &)  Ensemble " << tag << " exists in both but has different size: " << out.nSamplesEns(e) << " : " << r.nSamplesEns(i) << std::endl);
    }
  }
  return out;
}


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
    ::write(writer,ens_jacks,"ens_jacks",false); //don't flatten the internal vector
    ::write(writer,cen,"cen");
    ::write(writer,*layout,"layout");
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);
    ::read(reader,ens_jacks,"ens_jacks",false);
    ::read(reader,cen, "cen");

    static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));
    superJackknifeLayout* _layout = layouts.back().get();
    ::read(reader,*_layout,"layout");
    layout = _layout;
    reader.leave();
  }
#endif
  
};
#ifdef HAVE_HDF5
template<typename D>
inline void write(HDF5writer &writer, const superJackknifeDistribution<D> &sj, const std::string &tag){ sj.write(writer,tag); }
template<typename D>
inline void read(HDF5reader &reader, superJackknifeDistribution<D> &sj, const std::string &tag){ sj.read(reader,tag); }
#endif

template<typename T> //override for distributions with extra information
struct extraFlattenIO<superJackknifeDistribution<T> >{
  static inline void write(HDF5writer &writer, const std::vector<superJackknifeDistribution<T> > &value){
    assert(value.size() > 0);
    const superJackknifeLayout &layout = value[0].getLayout();
    if(value.size() > 1)
      for(int i=1;i<value.size();i++) assert(value[i].getLayout() == layout);
    ::write(writer,layout,"layout");
    
    std::vector<T> cen(value.size());
    for(int i=0;i<cen.size();i++) cen[i] = value[i].best();
    ::write(writer,cen,"central");
  }    
  static inline void read(HDF5reader &reader, std::vector<superJackknifeDistribution<T> > &value){
    static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));
    superJackknifeLayout* layout = layouts.back().get();
    ::read(reader,*layout,"layout");
    for(int i=0;i<value.size();i++) value[i].setLayout(*layout);
        
    std::vector<T> cen;
    ::read(reader,cen,"central");
    for(int i=0;i<cen.size();i++) value[i].best() = cen[i];
  }
  static inline void write(HDF5writer &writer, const std::vector<std::vector<superJackknifeDistribution<T> > > &value){
    assert(value.size() > 0 && value[0].size()>0);
    const superJackknifeLayout &layout = value[0][0].getLayout();
    for(int i=0;i<value.size();i++)
      for(int j= (i==0 ? 1:0); j<value[i].size(); j++)
	assert(value[i][j].getLayout() == layout);
    ::write(writer,layout,"layout");    
    
    int sz = 0; for(int i=0;i<value.size();i++) sz += value[i].size();
    
    std::vector<T> cen(sz);
    int off=0;
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++)
	cen[off++] = value[i][j].best();

    ::write(writer,cen,"central");
  }
  static inline void read(HDF5reader &reader, std::vector<std::vector<superJackknifeDistribution<T> > > &value){
    static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
    layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));
    superJackknifeLayout* layout = layouts.back().get();
    ::read(reader,*layout,"layout");
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++)
	value[i][j].setLayout(*layout);
    
    std::vector<T> cen;
    ::read(reader,cen,"central");

    int off=0;
    for(int i=0;i<value.size();i++) //has already been resized
      for(int j=0;j<value[i].size();j++)
	value[i][j].best() = cen[off++];
  } 
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


struct EnsembleData{
  std::string tag;
  std::string SampleType;
  int EnsembleSize;
  double avg;
  std::vector<double> values;
};
void read(XMLreader &reader, EnsembleData &v, const std::string &tag){
  //std::cout << "Reading EnsembleData with tag '" << tag << "'. Context contains:\n" << reader.printGroupEntries() << std::endl;
  reader.enter(tag);
  //std::cout << "Entered '" << tag << "'. Context now contains\n" << reader.printGroupEntries() << std::endl;
  read(reader,v.tag,"tag");
  read(reader,v.SampleType,"SampleType");
  read(reader,v.EnsembleSize,"EnsembleSize");
  read(reader,v.avg,"avg");
  read(reader,v.values,"values");
  reader.leave();
}

void read(XMLreader &reader, superJackknifeDistribution<double> &v, const std::string &tag){
  reader.enter(tag);
#define GETIT(type,name) type name; read(reader,name,#name)

  GETIT(std::string, SampleType);
  assert(SampleType == "SuperJackBoot");
  
  GETIT(int, Nmeas);
  GETIT(int, Nensembles);

  std::vector<EnsembleData> ens;
  read(reader,ens,"Ensembles");
  
  reader.leave();

  //Setup the layout and superjack
  assert(ens.size() == Nensembles);
  
  static std::vector<std::unique_ptr<superJackknifeLayout> > layouts; //all will be deleted at the end
  layouts.push_back(std::unique_ptr<superJackknifeLayout>(new superJackknifeLayout));  
  superJackknifeLayout* layout = layouts.back().get();

  int meas=0;
  double avg;
  for(int e=0;e<ens.size();e++){
    if(e==0) avg = ens[e].avg;
    else if(ens[e].avg != avg) error_exit(std::cout << "parseDoubleJackknifeXML expect average on ensemble " << ens[e].tag << ", " << ens[e].avg << " to be equal to average of other ensembles, " << avg << std::endl);
    
    meas += ens[e].EnsembleSize;
    if(ens[e].SampleType != "Jackknife") error_exit(std::cout << "read(XMLreader &reader, superJackknifeDistribution<double> &v, const std::string &tag) does not presently support sub-distributions that aren't jackknife\n");

    layout->addEnsemble(ens[e].tag,ens[e].EnsembleSize);    
  }
  assert(meas == Nmeas);

  v = superJackknifeDistribution<double>(*layout, avg);
  for(int e=0;e<ens.size();e++){
    jackknifeDistribution<double> tmp(ens[e].EnsembleSize);
    for(int s=0;s<ens[e].EnsembleSize;s++) tmp.sample(s) = ens[e].values[s];
    v.setEnsembleJackknife(e,tmp);
  }
#undef GETIT
}


#endif

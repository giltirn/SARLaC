#ifndef _SUPERMULTI_LAYOUT_H_
#define _SUPERMULTI_LAYOUT_H_

//To keep track of what sub-distribution belongs to what source, we define layout objects shared by instances of the supermulti class

#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>
#include<parser/parser.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(MultiType, (Jackknife)(Bootstrap));
GENERATE_HDF5_ENUM_SERIALIZE(MultiType);

class superMultiLayout{
  std::vector<std::string > ens_tags;
  std::vector<MultiType> ens_types;
  std::vector<int> ens_sizes;
  std::vector< std::pair<int,int> > ens_sample_map; //map a global sample index to a pair { ens_idx, sub_idx }   where sub_idx is the index within the subensemble
  int total_size;

  inline void clear(){
    ens_tags.resize(0);
    ens_types.resize(0);
    ens_sizes.resize(0);
    ens_sample_map.resize(0);
    total_size = 0;
  }
    
public:
  superMultiLayout(): total_size(0){}
  
  //returns -1 for non-existent
  inline int ensIdx(const std::string &e) const{
    for(int i=0;i<ens_tags.size();i++) if(e == ens_tags[i]) return i;
    return -1;
  }
  inline const std::string &ensTag(const int i) const{ return ens_tags[i]; }

  void setEnsTag(const int i, const std::string &to){ ens_tags[i] = to; } //use only if you know what you are doing!

  void setEnsSize(const int i, const int new_sz){  //use only if you know what you are doing!
    ens_sizes[i] = new_sz; 

    //Regenerate the total size and ensemble map
    total_size = 0;
    ens_sample_map.clear();
    for(int e=0;e<ens_sizes.size();e++){
      total_size += ens_sizes[e];
      for(int i=0;i<ens_sizes[e];i++){
	ens_sample_map.push_back(std::pair<int,int>(e,i) );
      }
    }
  }

  void addEnsemble(const std::string &tag, const MultiType type, const int size){
    int ens_idx = ens_tags.size();
    ens_tags.push_back(tag);
    ens_types.push_back(type);
    ens_sizes.push_back(size);
    for(int i=0;i<size;i++) ens_sample_map.push_back(std::pair<int,int>(ens_idx,i) );
    total_size += size;
  }

  inline size_t offset(const int ens_idx) const{ size_t r = 0; for(int i=0;i<ens_idx;i++) r += ens_sizes[i]; return r; }

  inline int nSamplesTotal() const{ return total_size; }

  inline int nEnsembles() const{ return ens_tags.size(); }

  inline int nSamplesEns(const int ens_idx) const{ return ens_sizes[ens_idx]; }
  inline int nSamplesEns(const std::string &ens_tag) const{ return nSamplesEns(ensIdx(ens_tag)); }
  
  inline MultiType ensType(const int ens_idx) const{ return ens_types[ens_idx]; }
  inline MultiType ensType(const std::string &ens_tag) const{ return ens_types[ensIdx(ens_tag)]; }

  inline const std::pair<int,int> & sampleMap(int gsample) const{
    return ens_sample_map[gsample];
  }

  bool operator==(const superMultiLayout &r) const{
    if(&r == this) return true;
    if(total_size != r.total_size) return false;
    for(int e=0;e<ens_sizes.size();e++) if(ens_sizes[e] != r.ens_sizes[e]) return false;
    for(int e=0;e<ens_sample_map.size();e++) if(ens_sample_map[e] != r.ens_sample_map[e]) return false;
    for(int e=0;e<ens_tags.size();e++) if(ens_tags[e] != r.ens_tags[e]) return false;
    for(int e=0;e<ens_types.size();e++) if(ens_types[e] != r.ens_types[e]) return false;
    return true;
  }
  inline bool operator!=(const superMultiLayout &r) const{ return !(*this == r); }
  
  
  //r is considered 'equivalent' to *this if all the subensembles inside r are contained inside *this
  //Note this allows *this to contain other ensembles
  bool equiv(const superMultiLayout &r) const{
    if(&r == this) return true;

    std::map<std::string, std::pair<MultiType, int> > this_subens;
    for(int e=0;e< nEnsembles() ;e++){
      this_subens[ ens_tags[e] ] = { ens_types[e], ens_sizes[e] };
    }

    for(int e=0; e< r.nEnsembles(); e++){
      auto it = this_subens.find( r.ens_tags[e] );
      if(it == this_subens.end()) return false;
      if(it->second.first != r.ens_types[e]) return false;
      if(it->second.second != r.ens_sizes[e]) return false;
    }
    return true;
  }


#ifdef HAVE_HDF5
  void write(HDF5writer &writer, const std::string &tag) const{
    writer.enter(tag);
    CPSfit::write(writer,ens_tags,"ens_tags");
    CPSfit::write(writer,ens_types,"ens_types");
    CPSfit::write(writer,ens_sizes,"ens_sizes");
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);
    std::vector<std::string > tags;
    std::vector<int> sizes;
    std::vector<MultiType> types;
    
    CPSfit::read(reader,tags,"ens_tags");
    CPSfit::read(reader,types,"ens_types");
    CPSfit::read(reader,sizes,"ens_sizes");

    clear();
    for(int i=0;i<tags.size();i++) addEnsemble(tags[i],types[i],sizes[i]);

    reader.leave();
  }
#endif
};
GENERATE_HDF5_SERIALIZE_FUNC(superMultiLayout);

inline std::ostream & operator<<(std::ostream &os, const superMultiLayout &l){
  os << "{ samples:"<<l.nSamplesTotal() << ", #ens:" << l.nEnsembles() << ", ensembles:";
  for(int i=0;i<l.nEnsembles();i++) os << "(" << l.ensTag(i) << ", type: " << toString(l.ensType(i)) << ", size:" << l.nSamplesEns(i) << ")";
  os << "}";
  return os;
}

superMultiLayout combine(const superMultiLayout &l, const superMultiLayout &r){
  superMultiLayout out(l);
  for(int i=0;i<r.nEnsembles();i++){
    const std::string &tag = r.ensTag(i);
    int e = out.ensIdx(tag);
    if(e == -1){ //doesn't currently exist
      out.addEnsemble(tag, r.ensType(i), r.nSamplesEns(i));
    }else if(out.nSamplesEns(e) != r.nSamplesEns(i)){
      error_exit(std::cout << "combine(const superMultiLayout &, const superMultiLayout &)  Ensemble " << tag << " exists in both but has different size: " << out.nSamplesEns(e) << " : " << r.nSamplesEns(i) << std::endl);
    }else if(out.ensType(e) != r.ensType(i)){
      error_exit(std::cout << "combine(const superMultiLayout &, const superMultiLayout &)  Ensemble " << tag << " exists in both but has different types: " << toString(out.ensType(e)) << " : " << toString(r.ensType(i)) << std::endl);
    }
  }
  return out;
}

//A class that acts as a global container for layouts with ownership semantics
class superMultiLayoutManagerDef{
  std::set<superMultiLayout*> layouts;
  std::map<std::string, std::set<superMultiLayout*>::iterator > tag_map;
public:

  std::set<superMultiLayout*>::const_iterator begin() const{ return layouts.begin(); }
  std::set<superMultiLayout*>::const_iterator end() const{ return layouts.end(); }
  
  void addLayout(superMultiLayout *l, const std::string &tag){
    auto r = layouts.insert(l);
    if(!r.second){ //already exists; check tags / pointers match in tag_map
      auto it = tag_map.find(tag);
      if(it == tag_map.end() || *(it->second) != l) error_exit(std::cout << "Attempt to insert a layout twice with different tags!" << std::endl); 
    }else{
      auto it = tag_map.find(tag);
      if(it != tag_map.end()) error_exit(std::cout << "Attempt to insert two different layouts with the same tag" << std::endl);
      tag_map[tag] = r.first;
    }
  }

  //If an existing layout exists which is 'equivalent' to l(see above), return the pointer to it
  //Otherwise return nullptr
  superMultiLayout* findEquivalent(superMultiLayout const* l) const{
    for(auto it = layouts.begin(); it != layouts.end(); it++){
      if( (*it)->equiv(*l) ){
	return *it;
      }
    }
    return nullptr;
  }   

  superMultiLayout* getLayout(const std::string &tag) const{
    auto it = tag_map.find(tag);
    if(it == tag_map.end())  error_exit(std::cout << "getLayout cannot find tag: " << tag << std::endl);
    return *(it->second);
  }
  const std::string & getTag(superMultiLayout const* ptr){
    for(auto e = tag_map.begin(); e != tag_map.end(); e++){
      if(*(e->second) == ptr){
	return e->first;
      }
    }
    error_exit(std::cout << "getTag cannot find ptr: " << ptr << std::endl);
  }
  
  //Release from ownership
  superMultiLayout* releaseLayout(const std::string &tag){
    auto it = tag_map.find(tag);
    if(it == tag_map.end())  error_exit(std::cout << "releaseLayout cannot find tag: " << tag << std::endl);
    superMultiLayout* ptr = *(it->second);
    layouts.erase(it->second);
    tag_map.erase(it);
    return ptr;
  }
  void releaseLayout(superMultiLayout const* ptr){
    for(auto e = tag_map.begin(); e != tag_map.end(); e++){
      if(*(e->second) == ptr){
	auto it = layouts.find(const_cast<superMultiLayout*>(ptr));
	if(it == layouts.end()) error_exit(std::cout << "releaseLayout corrupted state" << std::endl);
	tag_map.erase(e);
	layouts.erase(it);
	return;
      }
    }
    error_exit(std::cout << "releaseLayout cannot find ptr: " << ptr << std::endl);
  }

  ~superMultiLayoutManagerDef(){
    for(auto e = layouts.begin(); e!= layouts.end(); e++){
      delete *e;
    }
  }

};

superMultiLayoutManagerDef &  superMultiLayoutManager(){ static superMultiLayoutManagerDef d; return d; }

CPSFIT_END_NAMESPACE
#endif

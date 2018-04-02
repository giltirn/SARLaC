#ifndef _SUPERJACKKNIFE_LAYOUT_H_
#define _SUPERJACKKNIFE_LAYOUT_H_

//To keep track of what sub-distribution belongs to what source, we define layout objects shared by instances of the superjackknife class

#include<vector>

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>

CPSFIT_START_NAMESPACE

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
    CPSfit::write(writer,ens_tags,"ens_tags");
    CPSfit::write(writer,ens_sizes,"ens_sizes");
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);
    std::vector<std::string > tags;
    std::vector<int> sizes;
    
    CPSfit::read(reader,tags,"ens_tags");
    CPSfit::read(reader,sizes,"ens_sizes");

    clear();
    for(int i=0;i<tags.size();i++) addEnsemble(tags[i],sizes[i]);

    reader.leave();
  }
#endif
};
GENERATE_HDF5_SERIALIZE_FUNC(superJackknifeLayout);

inline std::ostream & operator<<(std::ostream &os, const superJackknifeLayout &l){
  os << "{ samples:"<<l.nSamplesTotal() << ", #ens:" << l.nEnsembles() << ", ensembles:";
  for(int i=0;i<l.nEnsembles();i++) os << "(" << l.ensTag(i) << ", size:" << l.nSamplesEns(i) << ")";
  os << "}";
  return os;
}

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

CPSFIT_END_NAMESPACE
#endif

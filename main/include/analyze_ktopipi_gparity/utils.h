#ifndef _ANALYZE_KTOPIPI_UTILS_H_
#define _ANALYZE_KTOPIPI_UTILS_H_

template<typename DistributionType>
void readFromXML(DistributionType &into, const int param, const std::string &file, const std::string &descr){
  UKvalenceDistributionContainer<DistributionType> tmp;
  XMLreader rd(file);
  read(rd, tmp, "data_in_file");
  into = std::move(tmp.list[param]);
  std::cout << "Read " << descr << " = " << into << std::endl;
}
template<typename DistributionType>
inline void readFromXML(DistributionType &into, const FileIdxPair &p, const std::string &descr){ readFromXML(into, p.idx, p.file, descr); }

template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type = 0>
void readFromHDF5(DistributionType &into, const int param, const std::string &file, const std::string &descr){  
  std::vector<DistributionType> tmp;
  readParamsStandard(tmp,file);
  into = std::move(tmp[param]);
  std::cout << "Read " << descr << " = " << into << std::endl;
}
template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type = 0>
inline void readFromHDF5(DistributionType &into, const FileIdxPair &p, const std::string &descr){ readFromHDF5(into, p.idx, p.file, descr); }

inline void writeParamsStandard(const NumericTensor<superMultiDistribution<double>,1> &params, const std::string &filename){
  writeParamsStandard(params.internalVector(),filename);
}

static std::vector<std::unique_ptr<superMultiLayout> > layout_store;

superMultiLayout* getLayout(bool &is_new, const std::string &ens_tag){
  superMultiLayout* layout = NULL;
  for(int i=0;i<layout_store.size();i++){
    if(layout_store[i]->ensIdx(ens_tag) != -1){
      layout = layout_store[i].get();
      is_new = false;
      break;
    }
  }
  if(layout == NULL){ 
    layout_store.push_back( std::unique_ptr<superMultiLayout>(new superMultiLayout()) ); 
    layout = layout_store.back().get(); 
    is_new = true;
  }
  return layout;
}


void readFromHDF5(superMultiDistribution<double> &into, const int param, const std::string &file, const std::string &descr, const std::string &ens_tag){  
  DistributionTypeEnum type;
  int vector_depth;
  getTypeInfo(type,vector_depth,file);
  assert(vector_depth == 1);

  bool is_new; 
  superMultiLayout* layout = getLayout(is_new, ens_tag);

#define CD(TYPE, TYPENM) \
  TYPE<double> tmp; \
  readFromHDF5(tmp, param, file, descr); \
  if(is_new) layout->addEnsemble(ens_tag, MultiType::TYPENM, tmp.size()); \
  into = superMultiDistribution<double>(*layout, ens_tag, tmp)

  if(type == DistributionTypeEnum::Jackknife){
    CD(jackknifeDistribution, Jackknife);
  }else if(type == DistributionTypeEnum::Bootstrap){
    CD(bootstrapDistribution, Bootstrap);
  }else assert(0);
  
#undef CD
}
inline void readFromHDF5(superMultiDistribution<double> &into, const FileIdxPair &p, const std::string &descr, const std::string &ens_tag){ readFromHDF5(into, p.idx, p.file, descr, ens_tag); }


void readFromHDF5(std::vector<superMultiDistribution<double> > &into, const std::string &file, const std::string &ens_tag){  
  DistributionTypeEnum type;
  int vector_depth;
  getTypeInfo(type,vector_depth,file);
  assert(vector_depth == 1);

  bool is_new; 
  superMultiLayout* layout = getLayout(is_new, ens_tag);

#define CD(TYPE, TYPENM) \
  std::vector<TYPE<double> > tmp;	   \
  readParamsStandard(tmp, file);					\
  if(is_new) layout->addEnsemble(ens_tag, MultiType::TYPENM, tmp[0].size()); \
  into.resize(tmp.size());						\
  for(int i=0;i<tmp.size();i++){					\
    into[i] = superMultiDistribution<double>(*layout, ens_tag, tmp[i]); \
  }

  if(type == DistributionTypeEnum::Jackknife){
    CD(jackknifeDistribution, Jackknife);
  }else if(type == DistributionTypeEnum::Bootstrap){
    CD(bootstrapDistribution, Bootstrap);
  }else assert(0);
  
#undef CD;
}

  
void readFromHDF5(std::vector<std::vector<superMultiDistribution<double> > > &into, const std::string &file, const std::string &ens_tag){  
  DistributionTypeEnum type;
  int vector_depth;
  getTypeInfo(type,vector_depth,file);
  assert(vector_depth == 2);

  bool is_new; 
  superMultiLayout* layout = getLayout(is_new, ens_tag);

#define CD(TYPE, TYPENM) \
  std::vector<std::vector<TYPE<double> > > tmp;				\
  readParamsStandard(tmp, file);					\
  if(is_new) layout->addEnsemble(ens_tag, MultiType::TYPENM, tmp[0][0].size()); \
  into.resize(tmp.size());						\
  for(int i=0;i<tmp.size();i++){					\
    into[i].resize(tmp[i].size());					\
    for(int j=0;j<tmp[i].size();j++){					\      
      into[i][j] = superMultiDistribution<double>(*layout, ens_tag, tmp[i][j]); \
    }									\
  }

  if(type == DistributionTypeEnum::Jackknife){
    CD(jackknifeDistribution, Jackknife);
  }else if(type == DistributionTypeEnum::Bootstrap){
    CD(bootstrapDistribution, Bootstrap);
  }else assert(0);
  
#undef CD;
}

  
#endif

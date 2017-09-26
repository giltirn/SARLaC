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

template<typename DistributionType>
void readFromHDF5(DistributionType &into, const int param, const std::string &file, const std::string &descr){  
  std::vector<DistributionType> tmp;
  readParamsStandard(tmp,file);
  into = std::move(tmp[param]);
  std::cout << "Read " << descr << " = " << into << std::endl;
}
template<typename DistributionType>
inline void readFromHDF5(DistributionType &into, const FileIdxPair &p, const std::string &descr){ readFromHDF5(into, p.idx, p.file, descr); }

inline void writeParamsStandard(const NumericTensor<superJackknifeDistribution<double>,1> &params, const std::string &filename){
  writeParamsStandard(params.internalVector(),filename);
}
  
#endif

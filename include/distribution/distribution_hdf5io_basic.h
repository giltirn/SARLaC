#ifndef _DISTRIBUTION_HDF5IO_BASIC_H
#define _DISTRIBUTION_HDF5IO_BASIC_H

//Serialize general distribution types to HDF5 format

#include<config.h>
#include<serialize/hdf5_serialize.h>

#ifdef HAVE_HDF5

CPSFIT_START_NAMESPACE

//Generic HDF5 IO for distributions. Overwrite if you want to save extra data
template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type = 0>
void write(HDF5writer &writer, const DistributionType &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.sampleVector(),"data");
  writer.leave();
}
template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type = 0>
void read(HDF5reader &reader, DistributionType &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.sampleVector(),"data");  
  reader.leave();
}


//For vector<distribution> and vector<vector<distribution> > of POD types optionally avoid the large metadata overheads by flattening
//This requires sufficient metadata to be written prior to the data block to reconstruct the input

//Specialize this class to add IO of appropriate metadata for flattened data structures
template<typename T>
struct extraFlattenIO{
  static inline void write(HDF5writer &writer, const std::vector<T> &value){}
  static inline void read(HDF5reader &reader, std::vector<T> &value){}
  static inline void write(HDF5writer &writer, const std::vector<std::vector<T> > &value){}
  static inline void read(HDF5reader &reader, std::vector<std::vector<T> > &value){} 
};

//Vectors of distributions
template<typename DistributionType, IF_DISTRIBUTION_NATIVE(DistributionType)>
void write(HDF5writer &writer, const std::vector<DistributionType> &value, const std::string &tag, bool flatten = true){
  writer.enter(tag);
  unsigned long size = value.size();
  writer.write(size,"size");  
  
  if(flatten){  
    int nsample = value.size() > 0 ? value[0].size() : 0;
    write(writer, nsample, "nsample");

    extraFlattenIO<DistributionType>::write(writer, value); //extra info
    
    std::vector<typename DistributionType::DataType> tmp(value.size() * nsample);
    for(int i=0;i<value.size();i++){
      assert(value[i].size() == nsample);
      for(int s=0;s<nsample;s++)
	tmp[s + nsample*i] = value[i].sample(s);
    }
    write(writer,tmp,"unrolled_data");
  }else{
    for(int i=0;i<value.size();i++){
      std::ostringstream os;
      os << "elem_" << i;
      write(writer,value[i],os.str());
    }
  }    
  writer.leave();
}
template<typename DistributionType, IF_DISTRIBUTION_NATIVE(DistributionType)>
void read(HDF5reader &reader, std::vector<DistributionType> &value, const std::string &tag, bool flattened = true){
  reader.enter(tag);
  unsigned long size;
  reader.read(size,"size");
  value.resize(size);
    
  if(flattened){
    int nsample;
    read(reader, nsample, "nsample");

    extraFlattenIO<DistributionType>::read(reader, value); 
    
    std::vector<typename DistributionType::DataType> tmp;
    read(reader,tmp,"unrolled_data");
    
    for(int i=0;i<size;i++){
      value[i].resize(nsample);      
      for(int s=0;s<nsample;s++)
	value[i].sample(s) = tmp[s + nsample*i];
    }
  }else{
    for(int i=0;i<value.size();i++){
      std::ostringstream os;
      os << "elem_" << i;
      read(reader,value[i],os.str());
    }    
  }
  reader.leave();
}

//Vector-vector of distributions
template<typename T, IF_DISTRIBUTION_NATIVE(T)>
void write(HDF5writer &writer, const std::vector<std::vector<T> > &value, const std::string &tag, bool flatten = true){
  writer.enter(tag); //enter a group

  unsigned long size1 = value.size();
  write(writer, size1,"size1");
  
  std::vector<unsigned long> size2(value.size());
  for(int i=0;i<value.size();i++) size2[i] = value[i].size();  
  writeCompact(writer, size2,"size2");
  
  if(flatten){
    int nsample = value.size() != 0 && value[0].size() != 0 ? value[0][0].size() : 0;    

    write(writer, nsample, "nsample");

    extraFlattenIO<T>::write(writer, value);
    
    unsigned long total_size = 0;
    for(int i=0;i<value.size();i++)
      total_size += value[i].size() * nsample;

    std::vector<typename T::DataType> tmp(total_size);
    
    int off=0;
    for(int i=0;i<value.size();i++){
      for(int j=0;j<value[i].size();j++){
	assert(value[i][j].size() == nsample );
	for(int s=0;s<nsample;s++)
	  tmp[off++] = value[i][j].sample(s);
      }
    }
    write(writer,tmp,"unrolled_data");
  }else{
    for(int i=0;i<value.size();i++){
      for(int j=0;j<value[i].size();j++){    
	std::ostringstream os;
	os << "elem_" << i << "_" << j;
	write(writer,value[i][j],os.str());
      }
    }
  }
  writer.leave();
}

template<typename T, IF_DISTRIBUTION_NATIVE(T)>
void read(HDF5reader &reader, std::vector<std::vector<T> > &value, const std::string &tag, bool flattened = true){
  reader.enter(tag); //enter a group
  unsigned long size1;
  read(reader, size1,"size1");
  
  std::vector<unsigned long> size2;
  readCompact(reader, size2,"size2");

  value.resize(size1);
  for(int i=0;i<value.size();i++)
    value[i].resize(size2[i]);
  
  if(flattened){
    int nsample;
    read(reader, nsample, "nsample");
    
    std::vector<typename T::DataType> tmp;
    read(reader,tmp,"unrolled_data");

    extraFlattenIO<T>::read(reader, value);
    
    int off = 0;
    for(int i=0;i<value.size();i++){
      for(int j=0;j<value[i].size();j++){
	value[i][j].resize(nsample);
	for(int s=0;s<nsample;s++)
	  value[i][j].sample(s) = tmp[off++];
      }
    }
  }else{
    for(int i=0;i<value.size();i++){
      for(int j=0;j<value[i].size();j++){    
	std::ostringstream os;
	os << "elem_" << i << "_" << j;
	read(reader,value[i][j],os.str());
      }
    }
  }
  reader.leave();
}

CPSFIT_END_NAMESPACE

#endif
#endif

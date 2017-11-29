#ifndef _DISTRIBUTION_IO_H
#define _DISTRIBUTION_IO_H

#include<distribution.h>
#include<parser.h>
#include<hdf5_serialize.h>

#ifdef HAVE_HDF5

//Generic HDF5 IO for distributions. Overwrite if you want to save extra data
template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type = 0>
void write(HDF5writer &writer, const DistributionType &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value._data,"data");
  writer.leave();
}
template<typename DistributionType, typename std::enable_if<hasSampleMethod<DistributionType>::value, int>::type = 0>
void read(HDF5reader &reader, DistributionType &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value._data,"data");  
  reader.leave();
}
//Specialization for jackknifeCdistribution
template<typename T>
void write(HDF5writer &writer, const jackknifeCdistribution<T> &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.cen,"cen");
  write(writer,value._data,"data");
  writer.leave();
}
template<typename T>
void read(HDF5reader &reader, jackknifeCdistribution<T> &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.cen,"cen");
  read(reader,value._data,"data");  
  reader.leave();
}

//For vector<distribution> and vector<vector<distribution> > of POD types avoid the large metadata overheads by optionally flattening
template<typename T>
struct extraFlattenIO{
  static inline void write(HDF5writer &writer, const std::vector<T> &value){}
  static inline void read(HDF5reader &reader, std::vector<T> &value){}
  static inline void write(HDF5writer &writer, const std::vector<std::vector<T> > &value){}
  static inline void read(HDF5reader &reader, std::vector<std::vector<T> > &value){} 
};
template<typename T> //override for distributions with extra information
struct extraFlattenIO<jackknifeCdistribution<T> >{
  static inline void write(HDF5writer &writer, const std::vector<jackknifeCdistribution<T> > &value){
    std::vector<T> cen(value.size());
    for(int i=0;i<cen.size();i++) cen[i] = value[i].best();
    ::write(writer,cen,"central");
  }    
  static inline void read(HDF5reader &reader, std::vector<jackknifeCdistribution<T> > &value){
    std::vector<T> cen;
    ::read(reader,cen,"central");
    for(int i=0;i<cen.size();i++) value[i].best() = cen[i];
  }
  static inline void write(HDF5writer &writer, const std::vector<std::vector<jackknifeCdistribution<T> > > &value){
    int sz = 0; for(int i=0;i<value.size();i++) sz += value[i].size();
    
    std::vector<T> cen(sz);
    int off=0;
    for(int i=0;i<value.size();i++)
      for(int j=0;j<value[i].size();j++)
	cen[off++] = value[i][j].best();

    ::write(writer,cen,"central");
  }
  static inline void read(HDF5reader &reader, std::vector<std::vector<jackknifeCdistribution<T> > > &value){
    std::vector<T> cen;
    ::read(reader,cen,"central");

    int off=0;
    for(int i=0;i<value.size();i++) //has already been resized
      for(int j=0;j<value[i].size();j++)
	value[i][j].best() = cen[off++];
  } 
};


template<typename DistributionType, typename std::enable_if<isDistributionOfHDF5nativetype<DistributionType>::value, int>::type = 0>
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
template<typename DistributionType, typename std::enable_if<isDistributionOfHDF5nativetype<DistributionType>::value, int>::type = 0>
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

template<typename T, typename std::enable_if<isDistributionOfHDF5nativetype<T>::value, int>::type = 0>
void write(HDF5writer &writer, const std::vector<std::vector<T> > &value, const std::string &tag, bool flatten = true){
  writer.enter(tag); //enter a group

  unsigned long size1 = value.size();
  writer.write(size1,"size1");
  
  std::vector<unsigned long> size2(value.size());
  for(int i=0;i<value.size();i++) size2[i] = value[i].size();  
  writer.write(size2,"size2");
  
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
template<typename T, typename std::enable_if<isDistributionOfHDF5nativetype<T>::value, int>::type = 0>
void read(HDF5reader &reader, std::vector<std::vector<T> > &value, const std::string &tag, bool flattened = true){
  reader.enter(tag); //enter a group
  unsigned long size1;
  reader.read(size1,"size1");
  
  std::vector<unsigned long> size2;
  reader.read(size2,"size2");

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








//Convenience functions for HDF5 write of fit parameters into standard std::vector<DistributionType<double> > format
template<typename T, int hasResize>
struct _actionResize{
  static inline void doit(T &v, const int n){}
};

template<typename T>
struct _actionResize<T,1>{
  static inline void doit(T &v, const int n){
    v.resize(n);
  }
};

template<typename DistributionOfPODtype, typename DistributionOfStructType>
struct standardIOhelper{
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    const int nsample = from.size();
    to.resize(from.size());
    for(int s=0;s<nsample;s++) to.sample(s) = from.sample(s)(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample){
    to.resize(nsample);
    typedef typename DistributionOfStructType::DataType StructType;
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<nsample;s++) to.sample(s) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    const int nsample = from.size();
    for(int s=0;s<nsample;s++){
      to.sample(s)(param_idx) = from.sample(s);
    }
  }
};

template<typename PODtype, typename StructType>
struct standardIOhelper<jackknifeCdistribution<PODtype>, jackknifeCdistribution<StructType> >{
  typedef jackknifeCdistribution<PODtype> DistributionOfPODtype;
  typedef jackknifeCdistribution<StructType> DistributionOfStructType;

  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    const int nsample = from.size();
    to.resize(from.size());
    for(int s=0;s<nsample;s++) to.sample(s) = from.sample(s)(param_idx);
    to.best() = from.best()(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample){
    to.resize(nsample);
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<nsample;s++) to.sample(s) = base;
    to.best() = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    const int nsample = from.size();
    for(int s=0;s<nsample;s++){
      to.sample(s)(param_idx) = from.sample(s);
    }
    to.best()(param_idx) = from.best();
  }
};


template<typename PODtype, typename StructType>
struct standardIOhelper<doubleJackknifeDistribution<PODtype>, doubleJackknifeDistribution<StructType> >{
  typedef doubleJackknifeDistribution<PODtype> DistributionOfPODtype;
  typedef doubleJackknifeDistribution<StructType> DistributionOfStructType;
  
  static inline void extractStructEntry(DistributionOfPODtype &to, const DistributionOfStructType &from, const int param_idx){
    const int nsample = from.size();
    to.resize(from.size());
    for(int s=0;s<nsample;s++)
      for(int t=0;t<nsample-1;t++)      
	to.sample(s).sample(t) = from.sample(s).sample(t)(param_idx);
  }
  static inline void setupDistributionOfStruct(DistributionOfStructType &to, const int nparams, const int nsample){
    to.resize(nsample);
    StructType base;
    _actionResize<StructType, hasResizeMethod<StructType>::value>::doit(base, nparams);
    for(int s=0;s<nsample;s++)
      for(int t=0;t<nsample-1;t++)
	to.sample(s).sample(t) = base;
  }
      
  static inline void insertStructEntry(DistributionOfStructType &to, const DistributionOfPODtype &from, const int param_idx){
    const int nsample = from.size();
    for(int s=0;s<nsample;s++){
      for(int t=0;t<nsample-1;t++){
	to.sample(s).sample(t)(param_idx) = from.sample(s).sample(t);
      }
    }
  }
};

template<typename T, int>
struct getDataType{};

template<typename T>
struct getDataType<T,0>{ typedef void type; };

template<typename T>
struct getDataType<T,1>{ typedef typename T::DataType type; };

template<typename T>
struct getDataType<doubleJackknifeDistribution<T>,1>{ typedef T type; };

#define IO_ENABLE_IF_POD(D) typename std::enable_if< hasSampleMethod<D>::value && std::is_arithmetic<typename getDataType<D,hasDataType<D>::value>::type>::value, int>::type = 0
#define IO_ENABLE_IF_NOT_POD(D) typename std::enable_if<hasSampleMethod<D>::value && !std::is_arithmetic<typename getDataType<D,hasDataType<D>::value>::type>::value, int>::type = 0

template<typename T>
struct standardIOformat{
  inline static void write(HDF5writer &writer, const std::vector<T> &value, const std::string &tag){ ::write(writer,value,tag,false); } //no flattening
  inline static void read(HDF5reader &reader, std::vector<T> &value, const std::string &tag){ ::read(reader,value,tag,false); }
  inline static void write(HDF5writer &writer, const std::vector<std::vector<T> > &value, const std::string &tag){ ::write(writer,value,tag,false); } //no flattening
  inline static void read(HDF5reader &reader, std::vector<std::vector<T> > &value, const std::string &tag){ ::read(reader,value,tag,false); }  
};
template<typename T>
struct standardIOformat<doubleJackknifeDistribution<T> >{
  inline static void write(HDF5writer &writer, const std::vector<doubleJackknifeDistribution<T> > &value, const std::string &tag){ ::write(writer,value,tag); } //no option to flatten for double jack
  inline static void read(HDF5reader &reader, std::vector<doubleJackknifeDistribution<T> > &value, const std::string &tag){ ::read(reader,value,tag); }
  inline static void write(HDF5writer &writer, const std::vector<std::vector<doubleJackknifeDistribution<T> > > &value, const std::string &tag){ ::write(writer,value,tag); }
  inline static void read(HDF5reader &reader, std::vector<std::vector<doubleJackknifeDistribution<T> > > &value, const std::string &tag){ ::read(reader,value,tag); }  
};


//Single distributions get written to single-element vectors
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const DistributionOfPODtype &params, const std::string &filename){
  //std::cout << "writeParamsStandard: POD - Writing " << printType<DistributionOfPODtype>() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");
  std::vector<DistributionOfPODtype> out(1,params);
  standardIOformat<DistributionOfPODtype>::write(writer, out , "value");
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(DistributionOfPODtype &params, const std::string &filename){  
  std::vector<DistributionOfPODtype> in;
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  standardIOformat<DistributionOfPODtype>::read(reader, in, "value");
  assert(in.size() == 1);
  params = in[0];
}
//vectors of distributions just written directly to disk
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const std::vector<DistributionOfPODtype> &params, const std::string &filename){
  //std::cout << "writeParamsStandard: POD - Writing " << printType<std::vector<DistributionOfPODtype> >() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  standardIOformat<DistributionOfPODtype>::write(writer, params , "value");
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(std::vector<DistributionOfPODtype> &params, const std::string &filename){  
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  standardIOformat<DistributionOfPODtype>::read(reader, params , "value");
}
//As are vectors of vectors
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const std::vector<std::vector<DistributionOfPODtype> > &params, const std::string &filename){
  //std::cout << "writeParamsStandard: POD - Writing " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << " as " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << std::endl;
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");
  standardIOformat<DistributionOfPODtype>::write(writer, params , "value");   //don't flatten
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(std::vector<std::vector<DistributionOfPODtype> > &params, const std::string &filename){  
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  standardIOformat<DistributionOfPODtype>::read(reader, params , "value");
}


//Distributions of structs get written to vectors of length equal to the number of struct elements. Assumes data type of struct elements can be converted to double
template<typename DistributionOfStructType, IO_ENABLE_IF_NOT_POD(DistributionOfStructType)>   
void writeParamsStandard(const DistributionOfStructType &params, const std::string &filename){
  typedef typename DistributionOfStructType::DataType ParamsType; //assumes ParamsType has methods size() and operator(const int) which is standard for fit parameters
  typedef typename DistributionOfStructType::template rebase<double> DistributionOfPODtype;

  //std::cout << "writeParamsStandard: Non-POD - Writing " << printType<DistributionOfStructType>() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  
  const int nsample = params.size();
  assert(nsample > 0);
  const int np = params.sample(0).size();

  std::vector<DistributionOfPODtype> out(np, DistributionOfPODtype(nsample));
  for(int p=0;p<np;p++)
    standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::extractStructEntry(out[p], params, p);
  
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  standardIOformat<DistributionOfPODtype>::write(writer, out , "value");
}


template<typename DistributionOfStructType, IO_ENABLE_IF_NOT_POD(DistributionOfStructType)>
void readParamsStandard(DistributionOfStructType &params, const std::string &filename){
  typedef typename DistributionOfStructType::DataType ParamsType;  //assumes ParamsType has methods size() and operator(const int) which is standard for fit parameters.
  //If the ParamsType is variable-length then it must have a resize(const int) method.
  typedef typename DistributionOfStructType::template rebase<double> DistributionOfPODtype;
  
  std::vector<DistributionOfPODtype> in;
  HDF5reader reader(filename);

  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());
  
  standardIOformat<DistributionOfPODtype>::read(reader, in , "value");

  assert(in.size() > 0);
  const int np = in.size();
  const int nsample = in[0].size();
  standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::setupDistributionOfStruct(params, np, nsample);  
  
  for(int p=0;p<np;p++)
    standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::insertStructEntry(params, in[p], p);
}




//Vectors of distributions of structs get written to vector<vector<Distribution> >
template<typename DistributionOfStructType, IO_ENABLE_IF_NOT_POD(DistributionOfStructType)>
void writeParamsStandard(const std::vector<DistributionOfStructType> &params, const std::string &filename){
  typedef typename DistributionOfStructType::DataType ParamsType; //assumes ParamsType has methods size() and operator(const int) which is standard for fit parameters
  typedef typename DistributionOfStructType::template rebase<double> DistributionOfPODtype;

  //std::cout << "writeParamsStandard: Non-POD - Writing " << printType<std::vector<DistributionOfStructType> >() << " as " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << std::endl;
  
  const int nouter = params.size();
  assert(nouter > 0);
  const int nsample = params[0].size();
  assert(nsample > 0);
  const int np = params[0].sample(0).size();

  std::vector<std::vector<DistributionOfPODtype> > out(nouter, std::vector<DistributionOfPODtype>(np, DistributionOfPODtype(nsample)));
  for(int o=0;o<nouter;o++)
    for(int p=0;p<np;p++)
      standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::extractStructEntry(out[o][p], params[o], p);
        
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  standardIOformat<DistributionOfPODtype>::write(writer, out , "value");
}
template<typename DistributionOfStructType, IO_ENABLE_IF_NOT_POD(DistributionOfStructType)>
void readParamsStandard(std::vector<DistributionOfStructType> &params, const std::string &filename){
  typedef typename DistributionOfStructType::DataType ParamsType;  //assumes ParamsType has methods size() and operator(const int) which is standard for fit parameters.
  //If the ParamsType is variable-length then it must have a resize(const int) method.
  typedef typename DistributionOfStructType::template rebase<double> DistributionOfPODtype;
  
  std::vector<std::vector<DistributionOfPODtype> > in;
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  standardIOformat<DistributionOfPODtype>::read(reader, in , "value");

  const int nouter = in.size();
  assert(nouter > 0);
  const int np = in[0].size();
  assert(np > 0);
  const int nsample = in[0][0].size();
  
  params.resize(nouter);      
  for(int o=0;o<nouter;o++){
    standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::setupDistributionOfStruct(params[o], np, nsample);    
    for(int p=0;p<np;p++)
       standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::insertStructEntry(params[o], in[o][p], p);
  }
}

#undef IO_ENABLE_IF_POD
#undef IO_ENABLE_IF_NOT_POD

//Open a standard format hdf5 file and query the distribution type and vector depth (i.e. vector<dist>=1, vector<vector<dist> >=2)
GENERATE_ENUM_AND_PARSER( DistributionTypeEnum, (Jackknife)(JackknifeC)(Raw)(DoubleJackknife)(SuperJackknife)  );

void getTypeInfo(DistributionTypeEnum &type, int & vector_depth, const std::string &filename){
  HDF5reader rd(filename);
  std::string typestr;
  read(rd, typestr, "distributionType");
  rd.enter("value");
  if(rd.contains("size2")) vector_depth = 2;
  else vector_depth =1;

  if(typestr == "rawDataDistribution<double>"){
    type = Raw;
  }else if(typestr == "jackknifeDistribution<double>"){
    type = Jackknife;
  }else if(typestr == "jackknifeCdistribution<double>"){
    type = JackknifeC;
  }else if(typestr == "doubleJackknifeDistribution<double>"){
    type = DoubleJackknife;
  }else if(typestr == "superJackknifeDistribution<double>"){
    type = SuperJackknife;    
  }else error_exit(std::cout << "getTypeInfo type " << typestr << " unimplemented\n");
}


#endif  //HAVE_HDF5



//For XML IO compatible with UKfit "boot_print" results
template<typename DistributionType>
struct UKvalenceDistributionContainer{
  int Nentries;
  std::vector<DistributionType> list;
};
template<typename DistributionType>
void read(XMLreader &reader, UKvalenceDistributionContainer<DistributionType> &v, const std::string &tag){
  reader.enter(tag);
  read(reader,v.Nentries,"Nentries");
  read(reader,v.list,"list");
  reader.leave();
}

template<typename T>
void read(XMLreader &reader, jackknifeDistribution<T> &v, const std::string &tag){
  reader.enter(tag);
#define GETIT(type,name) type name; read(reader,name,#name)

  GETIT(std::string, SampleType);
  assert(SampleType == "Jackknife");
  
  GETIT(int, Nmeas);
  GETIT(double, avg);
  GETIT(std::vector<T>, values);
  
  if(values.size() != Nmeas) error_exit(std::cout << "read(XMLreader &, jackknifeDistribution<T> &, const std::string &) file states Nmeas=" << Nmeas << " but read " << values.size() << " samples!\n");
  v.resize(Nmeas);
  for(int i=0;i<Nmeas;i++) v.sample(i) = values[i];
  
#undef GETIT
  reader.leave();
}












#endif

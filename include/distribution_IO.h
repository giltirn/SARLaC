#ifndef _DISTRIBUTION_IO_H
#define _DISTRIBUTION_IO_H

#include<distribution.h>

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

//Single distributions get written to single-element vectors
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const DistributionOfPODtype &params, const std::string &filename){
  std::cout << "writeParamsStandard: POD - Writing " << printType<DistributionOfPODtype>() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");
  std::vector<DistributionOfPODtype> out(1,params);
  write(writer, out , "value");
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(DistributionOfPODtype &params, const std::string &filename){  
  std::vector<DistributionOfPODtype> in;
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  read(reader, in, "value");
  assert(in.size() == 1);
  params = in[0];
}
//vectors of distributions just written directly to disk
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const std::vector<DistributionOfPODtype> &params, const std::string &filename){
  std::cout << "writeParamsStandard: POD - Writing " << printType<std::vector<DistributionOfPODtype> >() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  write(writer, params , "value");
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(std::vector<DistributionOfPODtype> &params, const std::string &filename){  
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  read(reader, params , "value");
}
//As are vectors of vectors
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const std::vector<std::vector<DistributionOfPODtype> > &params, const std::string &filename){
  std::cout << "writeParamsStandard: POD - Writing " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << " as " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << std::endl;
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  write(writer, params , "value");
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(std::vector<std::vector<DistributionOfPODtype> > &params, const std::string &filename){  
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == printType<DistributionOfPODtype>());  
  
  read(reader, params , "value");
}


//Distributions of structs get written to vectors of length equal to the number of struct elements. Assumes data type of struct elements can be converted to double
template<typename DistributionOfStructType, IO_ENABLE_IF_NOT_POD(DistributionOfStructType)>   
void writeParamsStandard(const DistributionOfStructType &params, const std::string &filename){
  typedef typename DistributionOfStructType::DataType ParamsType; //assumes ParamsType has methods size() and operator(const int) which is standard for fit parameters
  typedef typename DistributionOfStructType::template rebase<double> DistributionOfPODtype;

  std::cout << "writeParamsStandard: Non-POD - Writing " << printType<DistributionOfStructType>() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  
  const int nsample = params.size();
  assert(nsample > 0);
  const int np = params.sample(0).size();

  std::vector<DistributionOfPODtype> out(np, DistributionOfPODtype(nsample));
  for(int p=0;p<np;p++)
    standardIOhelper<DistributionOfPODtype, DistributionOfStructType>::extractStructEntry(out[p], params, p);
  
  HDF5writer writer(filename);
  std::string type = printType<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  write(writer, out , "value");
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
  
  read(reader, in , "value");

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

  std::cout << "writeParamsStandard: Non-POD - Writing " << printType<std::vector<DistributionOfStructType> >() << " as " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << std::endl;
  
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
  write(writer, out , "value");
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
  
  read(reader, in , "value");

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

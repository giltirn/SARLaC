#ifndef _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_COMPOUND_IO_H_
#define _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_COMPOUND_IO_H_

//IO routines for distributions of compound (struct) types

#include<config.h>
#include<utils/macros.h>
#include<distribution/distribution_hdf5io_conventional/helper.h>
#include<distribution/distribution_hdf5io_conventional/io_format.h>

CPSFIT_START_NAMESPACE

#ifdef HAVE_HDF5

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
  std::string type = getDistributionTypeString<DistributionOfPODtype>();
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
  assert(type == getDistributionTypeString<DistributionOfPODtype>());
  
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
  std::string type = getDistributionTypeString<DistributionOfPODtype>();
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
  assert(type == getDistributionTypeString<DistributionOfPODtype>());  
  
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

#endif

CPSFIT_END_NAMESPACE
#endif

#ifndef _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_POD_IO_H_
#define _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_POD_IO_H_

//IO routines for distributions of POD types

#include<config.h>
#include<utils/macros.h>
#include<distribution/distribution_hdf5io_conventional/io_format.h>

SARLAC_START_NAMESPACE

#ifdef HAVE_HDF5

//Single distributions get written to single-element vectors
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const DistributionOfPODtype &params, const std::string &filename){
  //std::cout << "writeParamsStandard: POD - Writing " << printType<DistributionOfPODtype>() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  
  HDF5writer writer(filename);
  std::string type = getDistributionTypeString<DistributionOfPODtype>();
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
  assert(type == getDistributionTypeString<DistributionOfPODtype>());  
  
  standardIOformat<DistributionOfPODtype>::read(reader, in, "value");
  assert(in.size() == 1);
  params = in[0];
}


//vectors of distributions just written directly to disk
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const std::vector<DistributionOfPODtype> &params, const std::string &filename){
  //std::cout << "writeParamsStandard: POD - Writing " << printType<std::vector<DistributionOfPODtype> >() << " as " << printType<std::vector<DistributionOfPODtype> >() << std::endl;
  HDF5writer writer(filename);
  std::string type = getDistributionTypeString<DistributionOfPODtype>();
  write(writer, type, "distributionType");  
  standardIOformat<DistributionOfPODtype>::write(writer, params , "value");
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(std::vector<DistributionOfPODtype> &params, const std::string &filename){  
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  if(type != getDistributionTypeString<DistributionOfPODtype>()) error_exit(std::cout << "readParamsStandard(std::vector<DistributionOfPODtype> &params, const std::string &filename) distribution type mismatch. Got "
									    << type << ", my type " << getDistributionTypeString<DistributionOfPODtype>() << std::endl);
  
  standardIOformat<DistributionOfPODtype>::read(reader, params , "value");
}


//As are vectors of vectors
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void writeParamsStandard(const std::vector<std::vector<DistributionOfPODtype> > &params, const std::string &filename){
  //std::cout << "writeParamsStandard: POD - Writing " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << " as " << printType<std::vector<std::vector<DistributionOfPODtype> > >() << std::endl;
  HDF5writer writer(filename);
  std::string type = getDistributionTypeString<DistributionOfPODtype>();
  write(writer, type, "distributionType");
  standardIOformat<DistributionOfPODtype>::write(writer, params , "value");   //don't flatten
}
template<typename DistributionOfPODtype, IO_ENABLE_IF_POD(DistributionOfPODtype)>
void readParamsStandard(std::vector<std::vector<DistributionOfPODtype> > &params, const std::string &filename){  
  HDF5reader reader(filename);
  
  std::string type;
  read(reader, type, "distributionType");
  assert(type == getDistributionTypeString<DistributionOfPODtype>());  
  
  standardIOformat<DistributionOfPODtype>::read(reader, params , "value");
}

#endif

SARLAC_END_NAMESPACE
#endif

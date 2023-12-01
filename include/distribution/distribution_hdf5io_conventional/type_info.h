#ifndef _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_TYPE_INFO_H_
#define _SARLAC_DISTRIBUTION_HDF5IO_CONVENTIONAL_TYPE_INFO_H_

//Parse the type info metadata in advance to allow setup of appropriate containers

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>
#include<parser/parser.h>

SARLAC_START_NAMESPACE

#ifdef HAVE_HDF5

//Open a standard format hdf5 file and query the distribution type and vector depth (i.e. vector<dist>=1, vector<vector<dist> >=2)
GENERATE_ENUM_AND_PARSER( DistributionTypeEnum, (Jackknife)(JackknifeC)(Raw)(DoubleJackknife)(SuperJackknife)(Bootstrap)(SuperMulti)  );

DistributionTypeEnum getTypeEnum(const std::string &typestr){
  if(typestr == "rawDataDistribution<double>"){
    return DistributionTypeEnum::Raw;
  }else if(typestr == "jackknifeDistribution<double>"){
    return DistributionTypeEnum::Jackknife;
  }else if(typestr == "jackknifeCdistribution<double>"){
    return DistributionTypeEnum::JackknifeC;
  }else if(typestr == "doubleJackknifeDistribution<double>"){
    return DistributionTypeEnum::DoubleJackknife;
  }else if(typestr == "superJackknifeDistribution<double>"){
    return DistributionTypeEnum::SuperJackknife;    
  }else if(typestr == "superMultiDistribution<double>"){
    return DistributionTypeEnum::SuperMulti;    
  }else if(typestr == "bootstrapDistribution<double>"){
    return DistributionTypeEnum::Bootstrap;
  }else error_exit(std::cout << "getTypeInfo type " << typestr << " unimplemented\n");
}

void getTypeInfo(DistributionTypeEnum &type, int & vector_depth, const std::string &filename){
  HDF5reader rd(filename);
  std::string typestr;
  read(rd, typestr, "distributionType");
  rd.enter("value");
  if(rd.contains("size2")) vector_depth = 2;
  else vector_depth =1;
  type = getTypeEnum(typestr);
}

#endif

SARLAC_END_NAMESPACE

#endif

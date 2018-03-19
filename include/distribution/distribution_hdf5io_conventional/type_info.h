#ifndef _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_TYPE_INFO_H_
#define _CPSFIT_DISTRIBUTION_HDF5IO_CONVENTIONAL_TYPE_INFO_H_

//Parse the type info metadata in advance to allow setup of appropriate containers

#include<config.h>
#include<utils/macros.h>
#include<serialize/hdf5_serialize.h>
#include<parser/parser.h>

CPSFIT_START_NAMESPACE

#ifdef HAVE_HDF5

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

#endif

CPSFIT_END_NAMESPACE
#endif

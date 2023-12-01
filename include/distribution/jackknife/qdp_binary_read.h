#ifndef _SARLAC_JACKKNIFE_QDP_BINARY_READ_H_
#define _SARLAC_JACKKNIFE_QDP_BINARY_READ_H_

#include<config.h>
#include<utils/macros.h>
#include<serialize/qdp_binary_read.h>
#include<distribution/jackknife/class.h>

SARLAC_START_NAMESPACE

template<typename T>
void read(QDPbinaryReader &rd, jackknifeDistribution<T> &into){
  int sample_type;
  read(rd, sample_type);
  if(sample_type != 2) error_exit(std::cout << "read(QDPbinaryReader &, jackknifeDistribution<T> &) failed: Expected sample_type=2, got " << sample_type << std::endl);
  int nmeas;
  read(rd, nmeas);
  into.resize(nmeas);
  for(int i=0;i<nmeas;i++) read(rd,into.sample(i));
  T avg;
  read(rd,avg);
}

SARLAC_END_NAMESPACE

#endif

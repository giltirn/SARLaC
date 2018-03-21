#ifndef _CPSFIT_JACKKNIFE_QDP_BINARY_READ_H_
#define _CPSFIT_JACKKNIFE_QDP_BINARY_READ_H_

#include<config.h>
#include<utils/macros.h>
#include<serialize/qdp_binary_read.h>
#include<distribution/jackknife/class.h>

CPSFIT_START_NAMESPACE

template<typename T>
void read(QDPbinaryReader &rd, jackknifeDistribution<T> &into){
  int sample_type;
  read(rd, sample_type);
  assert(sample_type == 2);
  int nmeas;
  read(rd, nmeas);
  into.resize(nmeas);
  for(int i=0;i<nmeas;i++) read(rd,into.sample(i));
  T avg;
  read(rd,avg);
}

CPSFIT_END_NAMESPACE

#endif

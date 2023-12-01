#ifndef _SARLAC_BOOTSTRAP_QDP_BINARY_READ_H_
#define _SARLAC_BOOTSTRAP_QDP_BINARY_READ_H_

#include<config.h>
#include<utils/macros.h>
#include<serialize/qdp_binary_read.h>
#include<distribution/bootstrap/class.h>

SARLAC_START_NAMESPACE

template<typename T>
void read(QDPbinaryReader &rd, bootstrapDistribution<T> &into){
  int sample_type;
  read(rd, sample_type);
  assert(sample_type == 3);
  int boots;
  read(rd, boots);
  into.resize(boots);
  for(int i=0;i<boots;i++) read(rd,into.sample(i));
  read(rd,into.best());
}

SARLAC_END_NAMESPACE

#endif

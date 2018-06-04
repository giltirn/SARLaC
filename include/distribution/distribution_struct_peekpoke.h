#ifndef _DISTRIBUTION_STRUCT_PEEKPOKE
#define _DISTRIBUTION_STRUCT_PEEKPOKE

#include<utility>
#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

CPSFIT_START_NAMESPACE

template<typename DistributionOfStruct>
struct _distribution_struct_peekpoke_helper{
  typedef typename std::decay<decltype(std::declval<DistributionOfStruct>().sample(std::declval<int>()))>::type StructType;
  typedef typename std::decay<decltype(std::declval<StructType>()(std::declval<int>()))>::type BaseType;
  typedef BaseType StructType::* StructMemberPointer;
  typedef typename DistributionOfStruct::template rebase<BaseType> DistributionOfBaseType;
};

#define STRUCT_TYPE typename _distribution_struct_peekpoke_helper<DistributionOfStruct>::StructType
#define BASE_DATA_TYPE typename _distribution_struct_peekpoke_helper<DistributionOfStruct>::BaseType
#define STRUCT_MEMBER_POINTER typename _distribution_struct_peekpoke_helper<DistributionOfStruct>::StructMemberPointer
#define DISTRIBUTION_OF_BASETYPE typename _distribution_struct_peekpoke_helper<DistributionOfStruct>::DistributionOfBaseType


template<typename DistributionOfStruct>
DISTRIBUTION_OF_BASETYPE distributionStructPeek(const DistributionOfStruct &in, const int idx){
  DISTRIBUTION_OF_BASETYPE out(in.size());
  for(int s=0;s<in.size();s++) out.sample(s) = in.sample(s)(idx);
  return out;
}

template<typename DistributionOfStruct>
void distributionStructPoke(DistributionOfStruct &out, const DISTRIBUTION_OF_BASETYPE &in, const int idx){
  assert(out.size() == in.size());
  for(int s=0;s<out.size();s++) out.sample(s)(idx) = in.sample(s);
}


template<typename DistributionOfStruct>
DISTRIBUTION_OF_BASETYPE distributionStructPeek(const DistributionOfStruct &in, STRUCT_MEMBER_POINTER elem){
  DISTRIBUTION_OF_BASETYPE out(in.size());
  for(int s=0;s<in.size();s++) out.sample(s) = in.sample(s).*elem;
  return out;
}

template<typename DistributionOfStruct>
void distributionStructPoke(DistributionOfStruct &out, const DISTRIBUTION_OF_BASETYPE &in, STRUCT_MEMBER_POINTER elem){
  assert(out.size() == in.size());
  for(int s=0;s<out.size();s++) out.sample(s).*elem = in.sample(s);
}
  

#undef STRUCT_TYPE
#undef BASE_DATA_TYPE
#undef STRUCT_MEMBER_POINTER
#undef DISTRIBUTION_OF_BASETYPE


CPSFIT_END_NAMESPACE

#endif

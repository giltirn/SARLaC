#ifndef _DISTRIBUTION_STRUCT_PEEKPOKE
#define _DISTRIBUTION_STRUCT_PEEKPOKE

#include<utility>
#include<config.h>
#include<utils/macros.h>
#include<utils/template_wizardry.h>

SARLAC_START_NAMESPACE

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

#define DSTRUCT_IT iterate<DistributionOfStruct>
#define DBASE_IT iterate<DISTRIBUTION_OF_BASETYPE>

template<typename DistributionOfStruct>
DISTRIBUTION_OF_BASETYPE distributionStructPeek(const DistributionOfStruct &in, const int idx){
  DISTRIBUTION_OF_BASETYPE out(getElem<DistributionOfStruct>::common_properties(in));
  for(int s=0;s<DSTRUCT_IT::size(in);s++) DBASE_IT::at(s,out) = DSTRUCT_IT::at(s, in)(idx);
  return out;
}

template<typename DistributionOfStruct>
void distributionStructPoke(DistributionOfStruct &out, const DISTRIBUTION_OF_BASETYPE &in, const int idx){
  assert(out.size() == in.size());
  for(int s=0;s<DSTRUCT_IT::size(out);s++) DSTRUCT_IT::at(s, out)(idx) = DBASE_IT::at(s, in);
}


template<typename DistributionOfStruct>
DISTRIBUTION_OF_BASETYPE distributionStructPeek(const DistributionOfStruct &in, STRUCT_MEMBER_POINTER elem){
  DISTRIBUTION_OF_BASETYPE out(getElem<DistributionOfStruct>::common_properties(in));
  for(int s=0;s<DSTRUCT_IT::size(in);s++) DBASE_IT::at(s,out) = DSTRUCT_IT::at(s, in).*elem;
  return out;
}

template<typename DistributionOfStruct>
void distributionStructPoke(DistributionOfStruct &out, const DISTRIBUTION_OF_BASETYPE &in, STRUCT_MEMBER_POINTER elem){
  assert(out.size() == in.size());
  for(int s=0;s<DSTRUCT_IT::size(out);s++) DSTRUCT_IT::at(s, out).*elem = DBASE_IT::at(s, in);
}
  

#undef STRUCT_TYPE
#undef BASE_DATA_TYPE
#undef STRUCT_MEMBER_POINTER
#undef DISTRIBUTION_OF_BASETYPE
#undef DSTRUCT_IT
#undef DBASE_IT

SARLAC_END_NAMESPACE

#endif

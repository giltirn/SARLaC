#ifndef _GLOBAL_DATA_MAP_ENUMS_H
#define _GLOBAL_DATA_MAP_ENUMS_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

//An arbitrary data type can be treated by dividing the configurations up into the set of data with only symmetric quark momenta, the set with only asymmetric quark momenta,
//and the set containing both upon which we compute the sampleAMA correction.
//The user can divide up the available data how they wish, but the sets must not overlap

enum class SubensTag { AsymmAll, SymmAll, Correction, SymmOnly, AsymmOnly }; 

std::ostream & operator<<(std::ostream &os, SubensTag tag){
  switch(tag){
  case SubensTag::AsymmAll:
    os << "asymm all"; break;
  case SubensTag::SymmAll:
    os << "symm all"; break;
  case SubensTag::Correction:
    os << "correction"; break;
  case SubensTag::SymmOnly:
    os << "symm only"; break;
  case SubensTag::AsymmOnly:
    os << "asymm only"; break;
  }  
  return os;
}

enum class DataTag { AsymmOnly, SymmOnly, AsymmCorr, SymmCorr }; 

CPSFIT_END_NAMESPACE

#endif

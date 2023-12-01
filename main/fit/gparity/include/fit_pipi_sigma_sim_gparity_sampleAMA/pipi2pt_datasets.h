#ifndef PIPI2PT_DATASETS_H_
#define PIPI2PT_DATASETS_H_

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

//Pipi operators with pion in momentum set {+-1, +-1, +-1}
std::map<SubensTag, std::set<int> > getGroundPiPiSubsets(const DataMap &dmap){
  std::map<SubensTag, std::set<int> > subens;

  subens[SubensTag::AsymmAll] = setUnion( dmap.getOuterConfigsWithTag("orig"), dmap.getOuterConfigsWithTag("correction") );
  subens[SubensTag::SymmAll] = setUnion( dmap.getOuterConfigsWithTag("extended"), dmap.getOuterConfigsWithTag("correction") );

  //subens[SubensTag::Correction] = setIntersection(subens[SubensTag::AsymmAll], subens[SubensTag::SymmAll]);
  //Fix the configurations used for the correction to just those in the "correction" set
  subens[SubensTag::Correction] = dmap.getOuterConfigsWithTag("correction");

  subens[SubensTag::SymmOnly] = setComplement(subens[SubensTag::SymmAll], subens[SubensTag::Correction]);
  subens[SubensTag::AsymmOnly] = setComplement(subens[SubensTag::AsymmAll], setUnion(subens[SubensTag::Correction], subens[SubensTag::SymmOnly]) );

  assert( setIntersection(subens[SubensTag::Correction], subens[SubensTag::SymmOnly]).size() == 0);
  assert( setIntersection(subens[SubensTag::Correction], subens[SubensTag::AsymmOnly]).size() == 0);
  assert( setIntersection(subens[SubensTag::AsymmOnly], subens[SubensTag::SymmOnly]).size() == 0);

  for(auto it = subens.begin(); it != subens.end(); ++it)
    std::cout << "pipi 2pt " << it->first << " " << it->second.size() << std::endl;
  
  return subens;
}

std::map<DataTag, std::map<int, DataLocationInfo const*> > getGroundPiPiDataSubsets(const std::map<SubensTag, std::set<int> > &subens, const DataMap &dmap){
  std::map<DataTag, std::map<int, DataLocationInfo const*> > out;
  
  out[DataTag::AsymmOnly] = dmap.getInfoPtrMapping(subens.find(SubensTag::AsymmOnly)->second,
						   { "orig", "correction" } );
  out[DataTag::SymmOnly] = dmap.getInfoPtrMapping(subens.find(SubensTag::SymmOnly)->second,
						   { "extended", "correction" } );

  out[DataTag::AsymmCorr] = dmap.getInfoPtrMapping(subens.find(SubensTag::Correction)->second,
						   { "orig", "correction" } );

  out[DataTag::SymmCorr] = dmap.getInfoPtrMapping(subens.find(SubensTag::Correction)->second,
						  { "extended", "correction" } );
  return out;
}

SARLAC_END_NAMESPACE

#endif

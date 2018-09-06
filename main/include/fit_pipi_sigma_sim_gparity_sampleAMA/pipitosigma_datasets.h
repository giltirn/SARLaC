#ifndef PIPITOSIGMA_DATASETS_H_
#define PIPITOSIGMA_DATASETS_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

std::map<SubensTag, std::set<int> > getPiPiToSigmaSubsets(const DataMap &dmap){
  std::map<SubensTag, std::set<int> > subens;
  
  subens[SubensTag::AsymmAll] = dmap.getOuterConfigsWithTag("sigma");
  subens[SubensTag::SymmAll] = dmap.getOuterConfigsWithTag("extended");
  subens[SubensTag::Correction] = setIntersection(subens[SubensTag::AsymmAll], subens[SubensTag::SymmAll]);
  subens[SubensTag::SymmOnly] = setComplement(subens[SubensTag::SymmAll], subens[SubensTag::Correction]);
  subens[SubensTag::AsymmOnly] = setComplement(subens[SubensTag::AsymmAll], setUnion(subens[SubensTag::Correction], subens[SubensTag::SymmOnly]) );

  assert( setIntersection(subens[SubensTag::Correction], subens[SubensTag::SymmOnly]).size() == 0);
  assert( setIntersection(subens[SubensTag::Correction], subens[SubensTag::AsymmOnly]).size() == 0);
  assert( setIntersection(subens[SubensTag::AsymmOnly], subens[SubensTag::SymmOnly]).size() == 0);

  for(auto it = subens.begin(); it != subens.end(); ++it)
    std::cout << "pipi->sigma " << it->first << " " << it->second.size() << std::endl;

  return subens;
}

std::map<DataTag, std::map<int, DataLocationInfo const*> > getPiPiToSigmaDataSubsets(const std::map<SubensTag, std::set<int> > &subens, const DataMap &dmap){
  std::map<DataTag, std::map<int, DataLocationInfo const*> > out;
  
  out[DataTag::AsymmOnly] = dmap.getInfoPtrMapping(subens.find(SubensTag::AsymmOnly)->second,
						   { "sigma" } );
  out[DataTag::SymmOnly] = dmap.getInfoPtrMapping(subens.find(SubensTag::SymmOnly)->second,
						   { "extended" } );

  out[DataTag::AsymmCorr] = dmap.getInfoPtrMapping(subens.find(SubensTag::Correction)->second,
						   { "sigma" } );

  out[DataTag::SymmCorr] = dmap.getInfoPtrMapping(subens.find(SubensTag::Correction)->second,
						  { "extended" } );
  return out;
}

CPSFIT_END_NAMESPACE

#endif

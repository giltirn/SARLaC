#ifndef SIGMA2PT_DATASETS_H_
#define SIGMA2PT_DATASETS_H_

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

std::map<SubensTag, std::set<int> > getSigma2ptSubsets(const DataMap &dmap){
  std::map<SubensTag, std::set<int> > subens;

  //No AMA correction needed
  subens[SubensTag::AsymmAll] = subens[SubensTag::Correction] = subens[SubensTag::AsymmOnly] = std::set<int>();
  subens[SubensTag::SymmAll] = subens[SubensTag::SymmOnly] = setUnion( dmap.getOuterConfigsWithTag("sigma"), dmap.getOuterConfigsWithTag("extended") );
  
  for(auto it = subens.begin(); it != subens.end(); ++it)
    std::cout << "sigma 2pt " << it->first << " " << it->second.size() << std::endl;
  
  return subens;
}

std::map<DataTag, std::map<int, DataLocationInfo const*> > getSigma2ptDataSubsets(const std::map<SubensTag, std::set<int> > &subens, const DataMap &dmap){
  std::map<DataTag, std::map<int, DataLocationInfo const*> > out;
  
  out[DataTag::AsymmOnly] = std::map<int, DataLocationInfo const*>();
			
  out[DataTag::SymmOnly] = dmap.getInfoPtrMapping(subens.find(SubensTag::SymmOnly)->second,
						  { "sigma", "extended" } );

  out[DataTag::AsymmCorr] = std::map<int, DataLocationInfo const*>();
			

  out[DataTag::SymmCorr] = std::map<int, DataLocationInfo const*>();
			
  return out;
}

CPSFIT_END_NAMESPACE

#endif

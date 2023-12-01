#ifndef _FIT_KTOSIGMA_KTOPIPI_GPARITY_TIANLE_COMPARE_H_
#define _FIT_KTOSIGMA_KTOPIPI_GPARITY_TIANLE_COMPARE_H_

#include<config.h>
#include<utils/macros.h>

//Compare K->sigma data with Tianle

SARLAC_START_NAMESPACE

//ktosigma_pss_vs_*jk indexing is [qidx][tsep_idx]   tsep_idx=  10:0,12:1,14:2,16:3,18:4

void readTianle(std::array<std::vector<jackknifeCorrelationFunctionD>,10> &ktosigma_pss_vs_jk,
		std::array<std::vector<doubleJackknifeCorrelationFunctionD>,10>  &ktosigma_pss_vs_2jk,
		const std::string &file){
  HDF5reader rd(file);
  for(int i=0; i<10; i++){
    std::ostringstream os[2];
    os[0] << "ktosigma_pss_vs_jk" << i;
    read(rd, ktosigma_pss_vs_jk[i], os[0].str().c_str());
    os[1] << "ktosigma_pss_vs_2jk" << i;
    read(rd, ktosigma_pss_vs_2jk[i], os[1].str().c_str());
  }
}


void compareTianle(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &me_j,
		   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &me_dj,
		   const std::string &tianle_file){
  
  std::array<std::vector<jackknifeCorrelationFunctionD>,10> tianle_j;
  std::array<std::vector<doubleJackknifeCorrelationFunctionD>,10>  tianle_dj;
  readTianle(tianle_j, tianle_dj, tianle_file);
  
  static const std::map<int,int> tsep_map = {  {10,0}, {12,1}, {14,2}, {16,3}, {18,4} };
  
  for(int q=0;q<10;q++){
    for(int i=0;i<me_j[q].size();i++){
      int t = int(me_j[q].coord(i).t);
      int tsep_k_sigma = int(me_j[q].coord(i).tsep_k_pi);

      const jackknifeDistributionD &me = me_j[q].value(i);

      auto tsepit = tsep_map.find(tsep_k_sigma);
      if(tsepit == tsep_map.end()) error_exit(std::cout << "Tianle check could not find map entry for tsep " << tsep_k_sigma << std::endl);
      
      const int tianle_t = tianle_j[q][tsepit->second].coord(t);
      if(tianle_t != t) error_exit(std::cout << "Tianle check t disparity " << t << " " << tianle_t << std::endl);

      const jackknifeDistributionD &tianle = tianle_j[q][tsepit->second].value(t);
      
      jackknifeDistributionD diff = tianle - me;
      
      std::cout << q << " " << t << " " << tsep_k_sigma << " " << me << " " << tianle << " " << diff << std::endl;
    }
  } 
}


SARLAC_END_NAMESPACE

#endif

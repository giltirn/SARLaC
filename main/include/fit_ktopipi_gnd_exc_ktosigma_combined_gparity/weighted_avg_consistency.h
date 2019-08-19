#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_WEIGHTED_AVG_CONSISTENCY_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_WEIGHTED_AVG_CONSISTENCY_H

#include<config.h>
#include<utils/macros.h>

#include "resampled_data.h"

CPSFIT_START_NAMESPACE

template<typename DistributionType>
void checkWeightedAvgConsistency(const ResampledData<DistributionType> &data,
				 const InputParamArgs &args, const std::vector<PiPiOperator> &operators,
				 const int tmin_k_op){
  DistributionType mK, cK;
  {				 
    std::vector<DistributionType> p;
    readParamsStandard(p,  args.kaon2pt_fit_result);
    mK = p[args.idx_mK];
    cK = sqrt( p[args.idx_cK] );
  }
  
  int nq = data(operators[0]).size();

  typedef correlationFunction<amplitudeDataCoord, DistributionType> CorrFuncType;

  std::cout << "Examining consistency of data with weighted avg after dividing out kaon time dependence" << std::endl;

  for(int q=0;q<nq;q++){
    for(int o=0;o<operators.size();o++){
      const CorrFuncType &corr = data(operators[o])[q];
   
      std::map<int,std::vector<int> > equal_top_snk;
      for(int i=0;i<corr.size();i++){
	int tk_op = (int)corr.coord(i).t;
	int top_snk = corr.coord(i).tsep_k_pi - tk_op;

	if(tk_op >= tmin_k_op) equal_top_snk[top_snk].push_back(i);
      }
      
      for(auto it=equal_top_snk.begin(); it != equal_top_snk.end(); ++it){
	int top_snk = it->first;
	const std::vector<int> &didx = it->second;
	
	std::vector<DistributionType> dnrm(didx.size());
	std::vector<DistributionType const*> tonrm(didx.size());

	//Divide out cK * e^{-mK * t} so only a function of top_snk
	for(int i=0;i<didx.size();i++){
	  dnrm[i] = corr.value(didx[i])/cK/exp(-mK* corr.coord(didx[i]).t );
	  tonrm[i] = &dnrm[i];
	}
	DistributionType wavg = weightedAvg(tonrm);
	
	for(int i=0;i<didx.size();i++){
	  DistributionType diff = dnrm[i] - wavg;
	  double significance = diff.best()/diff.standardError();

	  std::cout << "Q" << q+1 << " " << operators[o] << " top_snk " << top_snk << " " << corr.coord(didx[i]) 
		    << " value " << dnrm[i] << " wavg ( " << didx.size() << " data )" << wavg << " diff " << diff << " = " << significance << " sigma" << std::endl;
	}
      }
    }
  }
}
      

CPSFIT_END_NAMESPACE

#endif

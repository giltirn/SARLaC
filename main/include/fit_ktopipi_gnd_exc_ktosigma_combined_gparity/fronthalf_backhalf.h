#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_FRONTHALF_BACKHALF_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_FRONTHALF_BACKHALF_H

#include<config.h>
#include<utils/macros.h>

#include "raw_data.h"

CPSFIT_START_NAMESPACE

//Do a front-half / back-half analysis of each data point

struct GetSubset{
  int start;
  int lessthan;

  GetSubset(const int start, const int lessthan): start(start), lessthan(lessthan){}
  
  void operator()(rawDataDistributionD &d) const{
    rawDataDistributionD out(lessthan - start);
    int j=0;
    for(int i=start;i<lessthan;i++) out.sample(j++) = d.sample(i);
    d = out;
  }
};

template<typename DistributionType, typename Resampler, typename ArgsType, typename CMDlineType>
void frontHalfBackHalfAnalysis(const RawData &raw,
			       const Resampler &resampler_fh, const Resampler &resampler_bh, //resampler may care about number of configurations etc (eg bootstrap)
			       const ArgsType &args, const CMDlineType &cmdline){
  std::cout << "Starting front-half/back-half analysis" << std::endl;

  int nsample = raw.nsample();
  int bh_start= nsample/2;
  RawData fh_raw(raw), bh_raw(raw);
  fh_raw.applyFunction(GetSubset(0, bh_start));
  bh_raw.applyFunction(GetSubset(bh_start, nsample));

  ResampledData<DistributionType> fh_r, bh_r;
  fh_r.resample(fh_raw, args, cmdline, "front half", resampler_fh);
  bh_r.resample(bh_raw, args, cmdline, "back half", resampler_bh);

  int nq = fh_r(args.operators[0]).size();
  
  std::cout << "Front-half/back-half analysis" << std::endl;

  for(int q=0;q<nq;q++){
    for(int o=0;o<args.operators.size();o++){
      typedef correlationFunction<amplitudeDataCoord, DistributionType> CorrFuncType;

      const CorrFuncType &fh_corr = fh_r(args.operators[o])[q];
      const CorrFuncType &bh_corr = bh_r(args.operators[o])[q];
      assert(fh_corr.size() == bh_corr.size());
      
      for(int i=0;i<fh_corr.size();i++){
	double diff = bh_corr.value(i).best() -  fh_corr.value(i).best();
	double derr = sqrt( pow(bh_corr.value(i).standardError(),2) + pow(fh_corr.value(i).standardError(),2) ); //formally uncorrelated

	double err_ratio = bh_corr.value(i).standardError() / fh_corr.value(i).standardError();

	std::cout << "Q" << q+1 << " " << args.operators[o] << " " << bh_corr.coord(i) 
		  << " FH " << fh_corr.value(i) << " BH " <<  bh_corr.value(i) 
		  << " error ratio BH/FH " << err_ratio 
		  << " diff BH-FH " << diff << " +- " << derr << " = " << diff/derr << " sigma" << std::endl; 
      }
    }
  }

  std::cout << "Finished front-half/back-half analysis" << std::endl;
}

			       



CPSFIT_END_NAMESPACE

#endif

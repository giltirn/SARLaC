#ifndef _KTOPIPI_SAMPLE_AMA_RESAMPLE_CORRECT_H_
#define _KTOPIPI_SAMPLE_AMA_RESAMPLE_CORRECT_H_

class sampleAMA_resampler{
  char ens;
  int nS;
  int nC;
public:
  sampleAMA_resampler(){}
  sampleAMA_resampler(char _ens, int _nS, int _nC): ens(_ens), nS(_nS), nC(_nC){}

  template<typename DistributionType>
  inline void resample(DistributionType &out, const rawDataDistributionD &in) const{ 
    out = sampleAMAresample<DistributionType>::resample(in,ens,nS,nC);
  }
};


template<typename resampledDistributionType>
resampledDistributionType resampleCorrect(const rawDataDistributionD &sloppy_S, const rawDataDistributionD &sloppy_C, const rawDataDistributionD &exact_C, 
					  const sampleAMA_resampler &resampler_S, const sampleAMA_resampler &resampler_C, const std::string &descr = ""){
  resampledDistributionType out_r, sloppy_S_r, sloppy_C_r, exact_C_r;
  resampler_S.resample(sloppy_S_r, sloppy_S);
  resampler_C.resample(sloppy_C_r, sloppy_C);
  resampler_C.resample(exact_C_r, exact_C);
  
  out_r = sloppy_S_r + exact_C_r - sloppy_C_r;

#ifdef PRINT_CORRECTION
  if(descr != ""){
    resampledDistributionType diff = out_r - sloppy_S_r;
    resampledDistributionType reldiff = (out_r - sloppy_S_r)/sloppy_S_r;
    std::cout << descr << " corrected:" << out_r << " sloppy:" << sloppy_S_r << " diff:" << diff << " reldiff:" << reldiff << std::endl;
  }
#endif
  return out_r;
}

//Utility for adding descriptions/printing for only a particular distribution type
template<typename T, typename U>
struct printOnlyIfType{
  inline static std::string str(const std::string &str){ return ""; }
};
template<typename T>
struct printOnlyIfType<T,T>{
  inline static std::string str(const std::string &str){ return str; }
};


#endif

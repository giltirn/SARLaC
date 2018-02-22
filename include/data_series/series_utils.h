#ifndef _SERIES_UTILS_H_
#define _SERIES_UTILS_H_

#include<utils/macros.h>

CPSFIT_START_NAMESPACE

template<typename ResampledDataSeriesType, typename RawDataSeriesType>
inline void resample(ResampledDataSeriesType &out, const RawDataSeriesType &in){
  out.resize(in.size());
  for(int i=0;i<in.size();i++){
    out.coord(i) = in.coord(i);
    out.value(i).resample(in.value(i));    
  }  
}

CPSFIT_END_NAMESPACE
#endif

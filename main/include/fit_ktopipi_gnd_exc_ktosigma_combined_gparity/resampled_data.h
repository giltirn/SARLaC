#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RESAMPLED_DATA_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RESAMPLED_DATA_H

#include<ktopipi_common/basis_convert.h>

CPSFIT_START_NAMESPACE

template<typename DistributionType>
class ResampledData{
  typedef std::vector<correlationFunction<amplitudeDataCoord, DistributionType> >  CorrFuncAllQ;

  std::map<PiPiOperator, CorrFuncAllQ> data;
public:
  bool contains(const PiPiOperator op) const{ return data.find(op) != data.end(); }
  
  CorrFuncAllQ & operator()(const PiPiOperator op){ 
    auto it = data.find(op); 
    if(it == data.end()){
      auto &v = data[op];
      v.resize(10);
      return v;
    }else{
      return it->second;
    }
  }
  CorrFuncAllQ const & operator()(const PiPiOperator op) const{ auto it = data.find(op); assert(it != data.end()); return it->second;}

  void convertBasis10to7(){
    for(auto it = data.begin(); it != data.end(); it++){
      CorrFuncAllQ conv(7);
      convert10to7(conv, it->second);
      it->second = std::move(conv);
    }
  }

};

bool doOp(const PiPiOperator op, const std::vector<PiPiOperator> &ops){
  return std::find(ops.begin(),ops.end(),op) != ops.end();
}
  
CPSFIT_END_NAMESPACE

#endif

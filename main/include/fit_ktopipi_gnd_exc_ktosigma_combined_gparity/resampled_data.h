#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RESAMPLED_DATA_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RESAMPLED_DATA_H

#include<ktopipi_common/basis_convert.h>

#include "raw_data.h"

CPSFIT_START_NAMESPACE

template<typename DistributionType>
class ResampledData{
  typedef std::vector<correlationFunction<amplitudeDataCoord, DistributionType> >  CorrFuncAllQ;

  std::map<PiPiOperator, CorrFuncAllQ> data;
public:
  GENERATE_HDF5_SERIALIZE_METHOD( (data) );

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
  template<typename ArgsType, typename CMDlineType>
  inline void resample(const RawData &raw, const PiPiOperator op, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr){
    raw.resample( (*this)(op), op, args, cmdline, descr);
  }
  template<typename ArgsType, typename CMDlineType>
  inline void resample(const RawData &raw, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr){
    for(int i=0;i<args.operators.size();i++)
      raw.resample( (*this)(args.operators[i]), args.operators[i], args, cmdline, descr);
  }

  template<typename ArgsType, typename CMDlineType, typename BinResampler>
  inline void resample(const RawData &raw, const PiPiOperator op, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const BinResampler &bin_resampler){
    raw.resample( (*this)(op), op, args, cmdline, descr, bin_resampler);
  }
  template<typename ArgsType, typename CMDlineType, typename BinResampler>
  inline void resample(const RawData &raw, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const BinResampler &bin_resampler){
    for(int i=0;i<args.operators.size();i++)
      raw.resample( (*this)(args.operators[i]), args.operators[i], args, cmdline, descr, bin_resampler);
  }


};
template<typename DistributionType>
inline void write(CPSfit::HDF5writer &writer, const ResampledData<DistributionType> &d, const std::string &tag){ d.write(writer,tag); }
template<typename DistributionType>
inline void read(CPSfit::HDF5reader &reader, ResampledData<DistributionType> &d, const std::string &tag){ d.read(reader,tag); }


  
CPSFIT_END_NAMESPACE

#endif

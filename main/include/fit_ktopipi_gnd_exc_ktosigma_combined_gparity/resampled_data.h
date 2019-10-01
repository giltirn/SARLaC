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

  void constrainSourceSinkSep(const PiPiOperator op, const std::vector<int> tsep_k_snk_keep){
    auto it = data.find(op);
    if(it != data.end()){
      CorrFuncAllQ &corr = it->second;
      int nq = corr.size();
      for(int q=0;q<nq;q++){
	auto cit = corr[q].begin();
	while(cit != corr[q].end()){
	  int tsep_k_snk = cit->first.tsep_k_pi;
	  bool not_in_list = std::find(tsep_k_snk_keep.begin(), tsep_k_snk_keep.end(), tsep_k_snk) == tsep_k_snk_keep.end();
	  if(not_in_list){
	    std::cout << op << " removing " << cit->first << " as tsep_k_pi not in input list" << std::endl;
	    cit = corr[q].remove(cit);
	  }else ++cit;
	}
      }
    }
  }
  
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
  inline void resample(const RawData &raw, const PiPiOperator op, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const double alpha_coeff=1.){
    raw.resample( (*this)(op), op, args, cmdline, descr, alpha_coeff);
  }
  template<typename ArgsType, typename CMDlineType>
  inline void resample(const RawData &raw, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const double alpha_coeff=1.){
    for(int i=0;i<args.operators.size();i++)
      raw.resample( (*this)(args.operators[i]), args.operators[i], args, cmdline, descr, alpha_coeff);
  }

  template<typename ArgsType, typename CMDlineType, typename BinResampler>
  inline void resample(const RawData &raw, const PiPiOperator op, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const BinResampler &bin_resampler, const double alpha_coeff=1.){
    raw.resample( (*this)(op), op, args, cmdline, descr, bin_resampler, alpha_coeff);
  }
  template<typename ArgsType, typename CMDlineType, typename BinResampler>
  inline void resample(const RawData &raw, const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const BinResampler &bin_resampler, const double alpha_coeff=1.){
    for(int i=0;i<args.operators.size();i++)
      raw.resample( (*this)(args.operators[i]), args.operators[i], args, cmdline, descr, bin_resampler, alpha_coeff);
  }


};
template<typename DistributionType>
inline void write(CPSfit::HDF5writer &writer, const ResampledData<DistributionType> &d, const std::string &tag){ d.write(writer,tag); }
template<typename DistributionType>
inline void read(CPSfit::HDF5reader &reader, ResampledData<DistributionType> &d, const std::string &tag){ d.read(reader,tag); }


  
CPSFIT_END_NAMESPACE

#endif

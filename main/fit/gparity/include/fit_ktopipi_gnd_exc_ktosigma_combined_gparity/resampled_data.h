#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RESAMPLED_DATA_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RESAMPLED_DATA_H

#include<ktopipi_common/basis_convert.h>

#include "raw_data.h"

CPSFIT_START_NAMESPACE

#define PRUNE_ELEM_MEM (PiPiOperator, op)(int, q)(amplitudeDataCoord, c)
struct pruneElem{
  _GENERATE_MEMBERS(PRUNE_ELEM_MEM);
  pruneElem(): op(PiPiOperator::PiPiGnd), q(-1), c(0,10){} //note  q=-1  applies to all q
};
_GENERATE_PARSER(pruneElem, PRUNE_ELEM_MEM);

#define PRUNE_ARGS_MEM (std::vector<pruneElem>, prune)
struct pruneArgs{
  _GENERATE_MEMBERS(PRUNE_ARGS_MEM);
  pruneArgs(): prune(1){}
};
_GENERATE_PARSER(pruneArgs, PRUNE_ARGS_MEM);


template<typename DistributionType>
class ResampledData{
public:
  typedef std::vector<correlationFunction<amplitudeDataCoord, DistributionType> >  CorrFuncAllQ;
private:
  std::map<PiPiOperator, CorrFuncAllQ> data;
public:
  GENERATE_HDF5_SERIALIZE_METHOD( (data) );

  void prune(const pruneArgs &pargs){  //remove particular data points from data set
    for(int p=0;p<pargs.prune.size();p++){
      auto op_it = data.find(pargs.prune[p].op);
      if(op_it != data.end()){
	int qstart=0, qlessthan=op_it->second.size();
	if(pargs.prune[p].q != -1){
	  qstart = pargs.prune[p].q;
	  qlessthan = qstart+1;
	}
	for(int q=qstart; q<qlessthan; q++){
	  auto data_it = op_it->second[q].begin();
	  while(data_it->first != pargs.prune[p].c && data_it != op_it->second[q].end()) ++data_it;
	  if(data_it != op_it->second[q].end()){
	    std::cout << "Pruning data " << pargs.prune[p] << " with q=" << q << std::endl;
	    op_it->second[q].remove(data_it);
	  }
	}
      }
    }
  }

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

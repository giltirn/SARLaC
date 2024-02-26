#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RAW_DATA_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RAW_DATA_H

#include<ktopipi_common/basis_convert.h>
#include "enums.h"

SARLAC_START_NAMESPACE

bool doOp(const PiPiOperator op, const std::vector<PiPiOperator> &ops){
  return std::find(ops.begin(),ops.end(),op) != ops.end();
}

#define COPYOPTS(INTO, NM) \
  INTO.load_amplitude_data = cmdline.load_##NM##_amplitude_data; \
  INTO.load_amplitude_data_file = cmdline.load_##NM##_amplitude_data_file; \
  INTO.save_amplitude_data = cmdline.save_##NM##_amplitude_data; \
  INTO.save_amplitude_data_file = cmdline.save_##NM##_amplitude_data_file; \
  INTO.read_opts.load_data_checkpoint = cmdline.load_##NM##_data_checkpoint; \
  INTO.read_opts.load_data_checkpoint_stub = cmdline.load_##NM##_data_checkpoint_stub; \
  INTO.read_opts.save_data_checkpoint = cmdline.save_##NM##_data_checkpoint; \
  INTO.read_opts.save_data_checkpoint_stub = cmdline.save_##NM##_data_checkpoint_stub

struct RawData{
  ProjectedBubbleData *bubble_data_gnd;
  ProjectedBubbleData *bubble_data_exc;
  ProjectedSigmaBubbleData *bubble_data_sigma;

  std::vector<RawKtoPiPiData *> raw_ktopipi_gnd; //[tsep_k_pi]
  std::vector<RawKtoPiPiData *> raw_ktopipi_exc;
  std::vector<RawKtoSigmaData *> raw_ktosigma;
      
  int nsample() const{
    if(bubble_data_gnd) return bubble_data_gnd->bubble({0}).size();
    if(bubble_data_exc) return bubble_data_exc->bubble({0}).size();
    if(bubble_data_sigma) return bubble_data_sigma->bubble({0}).size();
    assert(0); return -1;
  }

  template<typename Functor> //Functor should act on rawDataDistributionD
  void applyFunction(const Functor &func){
    if(bubble_data_gnd) bubble_data_gnd->applyFunction(func);
    if(bubble_data_exc) bubble_data_exc->applyFunction(func);
    if(bubble_data_sigma) bubble_data_sigma->applyFunction(func);

    for(int i=0;i<raw_ktopipi_gnd.size();i++) if(raw_ktopipi_gnd[i]) raw_ktopipi_gnd[i]->applyFunction(func);
    for(int i=0;i<raw_ktopipi_exc.size();i++) if(raw_ktopipi_exc[i]) raw_ktopipi_exc[i]->applyFunction(func);
    for(int i=0;i<raw_ktosigma.size();i++) if(raw_ktosigma[i]) raw_ktosigma[i]->applyFunction(func);
  }

  template<typename DistributionType>
  inline static DistributionType removeSamplesInRange(const DistributionType &raw, const int start, const int lessthan){
    int remove_num_samples = lessthan - start;
    assert(remove_num_samples < raw.size());
    DistributionType out(raw.size() - remove_num_samples);
    int s=0;
    for(int p=0;p<start;p++)
      out.sample(s++) = raw.sample(p);
    for(int p=lessthan;p<raw.size();p++)
      out.sample(s++) = raw.sample(p);
    return out;
  }

  void removeSamplesInRange(const int start, const int lessthan){
    std::cout << "RawData removing samples in range [" << start << ", " << lessthan << ")" << std::endl;
    applyFunction(
		  [&](auto &v){ 
		    v = removeSamplesInRange(v,start,lessthan);				  
		  }
		  );
    std::cout << "RawData samples removed" << std::endl;
  }


  template<typename ArgsType, typename CMDlineType>
  void read(const ArgsType &args, const CMDlineType &cmdline){
    readKtoPiPiAllDataOptions read_opts;
    read_opts.read_opts.include_V_diagram = !cmdline.exclude_V_diagram;

    if(doOp(PiPiOperator::PiPiGnd, args.operators)){
      COPYOPTS(read_opts, ktopipi);
      
      std::cout << "Reading K->pipi(111) data" << std::endl;

      //Read the bubble data
      bubble_data_gnd = new ProjectedBubbleData(args.data_dir, 
						args.pipi_bubble_file_fmt,
						args.traj_start, args.traj_inc, args.traj_lessthan, 
						args.Lt, args.tsep_pipi, args.pipi_bubble_pimom_proj, read_opts.read_opts);
      //Read the K->pipi data
      raw_ktopipi_gnd.resize(args.tsep_k_pi.size(),NULL);
      for(int i=0;i<args.tsep_k_pi.size();i++)
	raw_ktopipi_gnd[i] = new RawKtoPiPiData(args.tsep_k_pi[i], *bubble_data_gnd, args.data_dir, 
						args.ktopipi_type_file_fmt, args.ktopipi_type1_pimom_proj, 
						args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, args.tsep_pipi, read_opts.read_opts);
    }
    if(doOp(PiPiOperator::PiPiExc, args.operators)){
      COPYOPTS(read_opts, ktopipi_exc);
      
      std::cout << "Reading K->pipi(311) data" << std::endl;

      //Read the bubble data
      bubble_data_exc = new ProjectedBubbleData(args.data_dir, 
						args.pipi_bubble_file_fmt,
						args.traj_start, args.traj_inc, args.traj_lessthan, 
						args.Lt, args.tsep_pipi, args.pipi_exc_bubble_pimom_proj, read_opts.read_opts);
      //Read the K->pipi data
      raw_ktopipi_exc.resize(args.tsep_k_pi.size(),NULL);
      for(int i=0;i<args.tsep_k_pi.size();i++)
	raw_ktopipi_exc[i] = new RawKtoPiPiData(args.tsep_k_pi[i], *bubble_data_exc, args.data_dir, 
						args.ktopipi_exc_type_file_fmt, args.ktopipi_exc_type1_pimom_proj, 
						args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, args.tsep_pipi, read_opts.read_opts);
    }
    if(doOp(PiPiOperator::Sigma, args.operators)){
      COPYOPTS(read_opts, ktosigma);
      
      std::cout << "Reading K->sigma data" << std::endl;

      //Read the bubble data
      static const std::vector<std::pair<threeMomentum, double> > sigma_bub_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
											      { {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
											      { {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
											      { {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       };
      
      bubble_data_sigma = new ProjectedSigmaBubbleData(args.data_dir, 
						       args.sigma_bubble_file_fmt,
						       args.traj_start, args.traj_inc, args.traj_lessthan, 
						       args.Lt, sigma_bub_quarkmom_proj, read_opts.read_opts);
      //Read the K->pipi data
      raw_ktosigma.resize(args.tsep_k_pi.size(),NULL);
      for(int i=0;i<args.tsep_k_pi.size();i++)
	raw_ktosigma[i] = new RawKtoSigmaData(args.tsep_k_pi[i], *bubble_data_sigma, args.data_dir, 
					      args.ktosigma_type_file_fmt, 
					      args.traj_start, args.traj_inc, args.traj_lessthan, args.Lt, read_opts.read_opts);
    }
  }

  //Generate resampled data for each four-quark operator
  template<typename DistributionType, typename ArgsType, typename CMDlineType, typename BinResampler>
  void resample(std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &A0_all, const PiPiOperator op, 
		const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const BinResampler &bin_resampler, const double alpha_coeff = 1.) const{

    assert(doOp(op, args.operators));
    NumericTensor<DistributionType,1> A0_full_srcavg;
    computeQamplitudeOpts opt;
    opt.alpha_scale = alpha_coeff;
    
    if(cmdline.disable_vacuum_subtraction)
      opt.do_vacuum_subtraction = false;

    for(int x=0;x<args.tsep_k_pi.size();x++){
      int tsep_k_pi = args.tsep_k_pi[x];
      
      for(int q=0;q<10;q++){
	switch(op){
	case PiPiOperator::PiPiGnd:
	  A0_full_srcavg = computeQamplitude<DistributionType>(q, tsep_k_pi, *raw_ktopipi_gnd[x], *bubble_data_gnd, args.Lt, descr, bin_resampler, opt);
	  break;
	case PiPiOperator::PiPiExc:
	  A0_full_srcavg = computeQamplitude<DistributionType>(q, tsep_k_pi, *raw_ktopipi_exc[x], *bubble_data_exc, args.Lt, descr, bin_resampler, opt);
	  break;
	case PiPiOperator::Sigma:
	  A0_full_srcavg = computeQamplitude<DistributionType>(q, tsep_k_pi, *raw_ktosigma[x], *bubble_data_sigma, args.Lt, descr, bin_resampler, opt);
	  break;
	}
	for(int t=0;t<=tsep_k_pi;t++)
	  A0_all[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg(&t));
      }
    }
  }

  template<typename DistributionType, typename ArgsType, typename CMDlineType>
  void resample(std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &A0_all, const PiPiOperator op, 
		const ArgsType &args, const CMDlineType &cmdline, const std::string &descr, const double alpha_coeff = 1.) const{
    basicBinResampler bin_resampler(args.bin_size);  
    resample(A0_all, op, args, cmdline, descr, bin_resampler,alpha_coeff);
  }

  //Extract <0|P|K> averaged over all kaon source timeslices
  template<typename DistributionType, typename ArgsType, typename CMDlineType, typename BinResampler>
  void computeKtoVacuumMatrixElem(correlationFunction<double, DistributionType> &vac_P_K,  //[t]
				  const ArgsType &args, const CMDlineType &cmdline, 
				  const std::string &descr, const BinResampler &bin_resampler) const{
    const int Lt = args.Lt;

    const RawKtoPiPiData &raw = *raw_ktopipi_gnd[0]; //tsep_K_pi irrelevant, no pions!

    //source avg
    int nt = raw.mix4_alltK_nobub.size(1);
    int ntK=raw.mix4_alltK_nobub.size(0);
    std::vector<DistributionType> src_avg(nt);
    DistributionType tmp;

    for(int tKidx=0;tKidx<ntK;tKidx++){
      for(int t=0;t<nt;t++){
	bin_resampler.binResample(tmp, raw.mix4_alltK_nobub({tKidx,t}));
	double nrm = 1./ntK;
	if(tKidx==0) src_avg[t] = tmp * nrm;
	else src_avg[t] = src_avg[t] + tmp * nrm;
      }
    }

    vac_P_K.resize(nt);
    for(int t=0;t<nt;t++){
      vac_P_K.coord(t) = t;
      vac_P_K.value(t) = -3./sqrt(6.)*src_avg[t]; //CHECKME: is the normalization right?
    }
  }


  //Extract the <op|P|K> matrix element and the coefficients alpha (10-basis) for external analysis
  template<typename DistributionType, typename ArgsType, typename CMDlineType, typename BinResampler>
  void computeAlphaAndPseudoscalarMatrixElem(std::vector<std::vector<std::vector<DistributionType> > > &alpha_sep_q_t, //[tsep_k_snk_idx][q][t]
					     std::vector<correlationFunction<double, DistributionType> > &op_P_K,  //[tsep_k_snk_idx][t]
					     const PiPiOperator op, 
					     const ArgsType &args, const CMDlineType &cmdline, 
					     const std::string &descr, const BinResampler &bin_resampler) const{

    const int Lt = args.Lt;

    if(op == PiPiOperator::PiPiGnd || op == PiPiOperator::PiPiExc){
      alpha_sep_q_t.resize(args.tsep_k_pi.size());
      op_P_K.resize(args.tsep_k_pi.size());

      NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
      ProjectedBubbleData const* bub = op == PiPiOperator::PiPiGnd ? bubble_data_gnd : bubble_data_exc;
      NumericTensor<DistributionType,1> resampled_bubble = bub->binResample<DistributionType>(bin_resampler);

      for(int x=0;x<args.tsep_k_pi.size();x++){
	int tsep_k_pi = args.tsep_k_pi[x];
	const RawKtoPiPiData &raw = op == PiPiOperator::PiPiGnd ? *raw_ktopipi_gnd[x] : *raw_ktopipi_exc[x];
	
	alpha_sep_q_t[x].resize(10);

	for(int q=0;q<10;q++){
	  computeAlphaAndVacuumSubtractions(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
					    raw.A0_type4_alltK_nobub, raw.mix4_alltK_nobub, resampled_bubble ,q, raw.nonzerotK(4),tsep_k_pi,Lt,bin_resampler);      

	  alpha_sep_q_t[x][q].resize(tsep_k_pi+1); //0<=t<=tsep_k_pi
	  for(int t=0;t<=tsep_k_pi;t++)
	    alpha_sep_q_t[x][q][t] = -sqrt(6.)/3 * alpha_r(&t); //this is the *true* normalization for alpha (internally we use a more convenient normalization to avoid a similar inverse normalization of the matrix elements it multiplies
	}
      
	IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
	for(int i=3;i<=4;i++) mix_srcavg_r(i) = binResampleAverageMixDiagram<DistributionType>(raw.mix_alltK(i), bin_resampler);
  
	op_P_K[x].resize(tsep_k_pi+1);
	for(int t=0;t<=tsep_k_pi;t++){
	  op_P_K[x].coord(t) = t;
	  op_P_K[x].value(t) = -3./sqrt(6.)*(mix_srcavg_r(3)(&t) + mix_srcavg_r(4)(&t) - mix4_srcavg_vacsub_r(&t)); //note the inverse normalization
	}
      }
    }else{
      assert(op == PiPiOperator::Sigma);

      alpha_sep_q_t.resize(args.tsep_k_sigma.size());
      op_P_K.resize(args.tsep_k_sigma.size());

      NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
      NumericTensor<DistributionType,1> resampled_bubble = bubble_data_sigma->binResample<DistributionType>(bin_resampler);

      std::vector<int> nonzero_tK(Lt); for(int t=0;t<Lt;t++) nonzero_tK[t] = t;

      for(int x=0;x<args.tsep_k_sigma.size();x++){
	int tsep_k_sigma = args.tsep_k_sigma[x];
	const RawKtoSigmaData & raw = *raw_ktosigma[x];
	
	alpha_sep_q_t[x].resize(10);

	for(int q=0;q<10;q++){
	  computeAlphaAndVacuumSubtractions(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
					    raw.A0_type4_alltK_nobub, raw.mix4_alltK_nobub, resampled_bubble ,q, nonzero_tK,tsep_k_sigma,Lt,bin_resampler);      

	  alpha_sep_q_t[x][q].resize(tsep_k_sigma+1);
	  for(int t=0;t<=tsep_k_sigma;t++)
	    alpha_sep_q_t[x][q][t] = 2. * alpha_r(&t); //this is the *true* normalization for alpha (internally we use a more convenient normalization to avoid a similar inverse normalization of the matrix elements it multiplies
	}
      
	IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
	for(int i=3;i<=4;i++) mix_srcavg_r(i) = binResampleAverageMixDiagram<DistributionType>(raw.mix_alltK(i), bin_resampler);
  
	op_P_K[x].resize(tsep_k_sigma+1);
	for(int t=0;t<=tsep_k_sigma;t++){
	  op_P_K[x].coord(t) = t;
	  op_P_K[x].value(t) = 1./2*(-mix_srcavg_r(3)(&t) + mix_srcavg_r(4)(&t) - mix4_srcavg_vacsub_r(&t)); //note - sign on the mix3 diag is intentional
	}
      }
    }
  }      

  template<typename DistributionType, typename ArgsType, typename CMDlineType>
  void computeAlphaAndPseudoscalarMatrixElem(std::vector<std::vector<std::vector<DistributionType> > > &alpha_sep_q_t, //[tsep_k_snk_idx][q][t]
					     std::vector<correlationFunction<double, DistributionType> > &op_P_K,  //[tsep_k_snk_idx][t]
					     const PiPiOperator op, 
					     const ArgsType &args, const CMDlineType &cmdline, 
					     const std::string &descr) const{
    basicBinResampler bin_resampler(args.bin_size);  
    computeAlphaAndPseudoscalarMatrixElem(alpha_sep_q_t, op_P_K, op, args, cmdline, descr, bin_resampler);
  }

  void write(HDF5writer &writer, const std::string &tag) const{
    writer.enter(tag);
#define WRITE(T) writePointer(writer, T, #T)
    WRITE(bubble_data_gnd);
    WRITE(bubble_data_exc);
    WRITE(bubble_data_sigma);
    WRITE(raw_ktopipi_gnd);
    WRITE(raw_ktopipi_exc);
    WRITE(raw_ktosigma);
#undef WRITE
    writer.leave();
  }
  void read(HDF5reader &reader, const std::string &tag){
    reader.enter(tag);
#define READ(T) readPointer(reader, T, #T)
    READ(bubble_data_gnd);
    READ(bubble_data_exc);
    READ(bubble_data_sigma);
    READ(raw_ktopipi_gnd);
    READ(raw_ktopipi_exc);
    READ(raw_ktosigma);
#undef READ
    reader.leave();
  }


  RawData(): bubble_data_gnd(NULL), bubble_data_exc(NULL), bubble_data_sigma(NULL){}
  
  template<typename ArgsType, typename CMDlineType>
  RawData(const ArgsType &args, const CMDlineType &cmdline): RawData(){
    read(args,cmdline);
  }

  RawData(const RawData &r): RawData(){
    if(r.bubble_data_gnd) bubble_data_gnd = new ProjectedBubbleData(*r.bubble_data_gnd);
    if(r.bubble_data_exc) bubble_data_exc = new ProjectedBubbleData(*r.bubble_data_exc);
    if(r.bubble_data_sigma) bubble_data_sigma = new ProjectedSigmaBubbleData(*r.bubble_data_sigma);

    raw_ktopipi_gnd.resize(r.raw_ktopipi_gnd.size(), NULL);
    for(int i=0;i<r.raw_ktopipi_gnd.size();i++) 
      if(r.raw_ktopipi_gnd[i]) raw_ktopipi_gnd[i] = new RawKtoPiPiData(*r.raw_ktopipi_gnd[i]);

    raw_ktopipi_exc.resize(r.raw_ktopipi_exc.size(), NULL);
    for(int i=0;i<r.raw_ktopipi_exc.size();i++) 
      if(r.raw_ktopipi_exc[i]) raw_ktopipi_exc[i] = new RawKtoPiPiData(*r.raw_ktopipi_exc[i]);

    raw_ktosigma.resize(r.raw_ktosigma.size(), NULL);
    for(int i=0;i<r.raw_ktosigma.size();i++) 
      if(r.raw_ktosigma[i]) raw_ktosigma[i] = new RawKtoSigmaData(*r.raw_ktosigma[i]);
  }

  ~RawData(){
#define DEL(T) if(T) delete T
#define FORDEL(T) for(int i=0;i<T.size();i++) if(T[i]) delete T[i]    

    DEL(bubble_data_gnd); DEL(bubble_data_exc); DEL(bubble_data_sigma);
    FORDEL(raw_ktopipi_gnd); FORDEL(raw_ktopipi_exc); FORDEL(raw_ktosigma);

#undef DEL
#undef FORDEL
  }

  RawData & operator=(const RawData &r) = delete;
};

inline void write(SARLaC::HDF5writer &writer, const RawData &d, const std::string &tag){ d.write(writer,tag); }
inline void read(SARLaC::HDF5reader &reader, RawData &d, const std::string &tag){ d.read(reader,tag); }






  
SARLAC_END_NAMESPACE

#endif

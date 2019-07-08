#ifndef _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RAW_DATA_H
#define _FIT_KTOPIPI_GND_EXC_KTOSIGMA_GPARITY_RAW_DATA_H

#include<ktopipi_common/basis_convert.h>

CPSFIT_START_NAMESPACE

bool doOp(const PiPiOperator op, const std::vector<PiPiOperator> &ops){
  return std::find(ops.begin(),ops.end(),op) != ops.end();
}

struct RawData{
  ProjectedBubbleData *bubble_data_gnd;
  ProjectedBubbleData *bubble_data_exc;
  ProjectedSigmaBubbleData *bubble_data_sigma;

  std::vector<RawKtoPiPiData *> raw_ktopipi_gnd; //[tsep_k_pi]
  std::vector<RawKtoPiPiData *> raw_ktopipi_exc;
  std::vector<RawKtoSigmaData *> raw_ktosigma;
  
  void read(const Args &args, const CMDline &cmdline){
    readKtoPiPiAllDataOptions read_opts;

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
						args.ktopipi_type_file_fmt, args.ktopipi_exc_type1_pimom_proj, 
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

  template<typename DistributionType>
  void resample(std::vector<correlationFunction<amplitudeDataCoord, DistributionType> > &A0_all, const PiPiOperator op, 
		const Args &args, const CMDline &cmdline, const std::string &descr) const{
    assert(doOp(op, args.operators));

    basic_resampler resampler;    
    NumericTensor<DistributionType,1> A0_full_srcavg;

    for(int x=0;x<args.tsep_k_pi.size();x++){
      int tsep_k_pi = args.tsep_k_pi[x];
      
      for(int q=0;q<10;q++){
	switch(op){
	case PiPiOperator::PiPiGnd:
	  A0_full_srcavg = computeQamplitude<DistributionType>(q, tsep_k_pi, *raw_ktopipi_gnd[x], *bubble_data_gnd, args.Lt, descr, args.bin_size, resampler);
	  break;
	case PiPiOperator::PiPiExc:
	  A0_full_srcavg = computeQamplitude<DistributionType>(q, tsep_k_pi, *raw_ktopipi_exc[x], *bubble_data_exc, args.Lt, descr, args.bin_size, resampler);
	  break;
	case PiPiOperator::Sigma:
	  A0_full_srcavg = computeQamplitude<DistributionType>(q, tsep_k_pi, *raw_ktosigma[x], *bubble_data_sigma, args.Lt, descr, args.bin_size, resampler);
	  break;
	}
	for(int t=0;t<=tsep_k_pi;t++)
	  A0_all[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg(&t));
      }
    }
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
  
  RawData(const Args &args, const CMDline &cmdline): RawData(){
    read(args,cmdline);
  }

  ~RawData(){
#define DEL(T) if(T) delete T
#define FORDEL(T) for(int i=0;i<T.size();i++) if(T[i]) delete T[i]    

    DEL(bubble_data_gnd); DEL(bubble_data_exc); DEL(bubble_data_sigma);
    FORDEL(raw_ktopipi_gnd); FORDEL(raw_ktopipi_exc); FORDEL(raw_ktosigma);

#undef DEL
#undef FORDEL
  }
};

inline void write(CPSfit::HDF5writer &writer, const RawData &d, const std::string &tag){ d.write(writer,tag); }
inline void read(CPSfit::HDF5reader &reader, RawData &d, const std::string &tag){ d.read(reader,tag); }






  
CPSFIT_END_NAMESPACE

#endif

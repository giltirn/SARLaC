#include<ktopipi_common/ktopipi_common.h>

#include<fit_ktosigma_gparity/args.h>
#include<fit_ktosigma_gparity/cmdline.h>

CPSFIT_START_NAMESPACE

//read_data
//---------------------------------------------------------------------------


//type can be 2,3,4
type1234Data readKtoSigmaType(const int type, const int traj_start, const int traj_inc, const int traj_lessthan,
		      const int tsep_k_sigma, const int Lt, 
		      const std::string &data_dir, const std::string &file_fmt){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  if(type == 2 || type == 3){
    //File format expected <TRAJ> <TSEP_K_SIGMA>
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>"), subStringSpecify("<TSEP_K_SIGMA>") };
    subStringReplace repl(file_fmt, keys);

    type1234Data typedata(type,Lt,nsample,Sigma);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;

      std::ostringstream fname;
      fname << data_dir << '/';
      repl.replace(fname, { anyToStr(traj), anyToStr(tsep_k_sigma) });
      typedata.parse(fname.str(),i);
    }
    return typedata;
  }else if(type == 4){
    //File format expects <TRAJ>
    static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>") };
    subStringReplace repl(file_fmt, keys);

    type1234Data typedata(type,Lt,nsample,Sigma);
#pragma omp parallel for
    for(int i=0;i<nsample;i++){
      const int traj = traj_start + i*traj_inc;

      std::ostringstream fname;
      fname << data_dir << '/';
      repl.replace(fname, { anyToStr(traj) });
      typedata.parse(fname.str(),i);
    }
    return typedata;
  }else{
    error_exit(std::cout << "readKtoSigmaType invalid type " << type << std::endl);
  }
}

//For original and extended data, use bubble_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
//                                               	       { {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
//                      				       { {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
//                      				       { {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       }
NumericTensor<rawDataDistributionD,1> readProjectedSigmaBubble(const std::string &data_dir, const std::string &file_fmt,
							       const int traj_start, const int traj_inc, const int traj_lessthan, 
							       const int Lt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj){
  int nsample = (traj_lessthan - traj_start)/traj_inc;
  NumericTensor<rawDataDistributionD,1> out({Lt}, rawDataDistributionD(nsample,0.));
  
  //File format expects <TRAJ>  <MOM>
  static std::vector<subStringSpecify> keys { subStringSpecify("<TRAJ>"),  subStringSpecify("<MOM>") };
  subStringReplace repl(file_fmt, keys);
  
  for(int p=0;p<bubble_quarkmom_proj.size();p++){
    sigmaSelfContraction self(Lt,nsample);

#pragma omp parallel for
    for(int sample=0; sample < nsample; sample++){
      int traj = traj_start + sample * traj_inc;
      std::ostringstream filename; filename << data_dir << '/';
      repl.replace(filename, { anyToStr(traj), momStr(bubble_quarkmom_proj[p].first) });

      std::cout << "Parsing " << filename.str() << std::endl;
      self.parse(filename.str(), sample);
    }

    for(int t=0;t<Lt;t++)
      out({t}) = out({t}) + self(t) * bubble_quarkmom_proj[p].second; //project
  }

  return out;
}


//-------------------------------------------------------
//amplitude_data
//-----------------------------------------------------------

class ProjectedSigmaBubbleData{
  static NumericTensor<rawDataDistributionD,1> getProjectedSigmaBubble(const std::string &data_dir, const std::string &file_fmt,
								       const int traj_start, const int traj_inc, const int traj_lessthan, 
								       const int Lt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj, 
								       const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    NumericTensor<rawDataDistributionD,1> bubble;
    if(opt.load_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.load_data_checkpoint_stub << "_sigmabubble.hdf5";
      std::cout << "Loading checkpoint data for bubble from " << file.str() << std::endl;
      HDF5reader rd(file.str());
      readProjectedBubble(rd,bubble,"bubble");
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }else{
      bubble = readProjectedSigmaBubble(data_dir, file_fmt, traj_start,traj_inc,traj_lessthan, Lt, bubble_quarkmom_proj);
    }
    if(opt.save_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.save_data_checkpoint_stub << "_sigmabubble.hdf5";
      std::cout << "Saving checkpoint data for bubble to " << file.str() << std::endl;
      HDF5writer wr(file.str());
      writeProjectedBubble(wr,bubble,"bubble");
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }
    return bubble;
  }
public:
  NumericTensor<rawDataDistributionD,1> bubble;
  NumericTensor<rawDataDistributionD,1> bubble_binned;
  NumericTensor<jackknifeDistributionD,1> bubble_j;
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj;

  ProjectedSigmaBubbleData(){}

  template<typename Resampler>
  ProjectedSigmaBubbleData(const std::string &data_dir, const std::string &file_fmt, \
		  const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size,
		  const int Lt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj,
		  const Resampler &resampler, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    bubble = getProjectedSigmaBubble(data_dir,file_fmt,
				     traj_start,traj_inc,traj_lessthan,
				     Lt, bubble_quarkmom_proj, opt);
    bubble_binned = bin(bubble,bin_size);
    bubble_j = bubble_binned.transform(resampleFunctorGeneral<jackknifeDistributionD,rawDataDistributionD,Resampler>(resampler));
    bubble_dj = bubble_binned.transform(resampleFunctorGeneral<doubleJackknifeDistributionD,rawDataDistributionD,Resampler>(resampler));
  }
};
template<typename DistributionType> struct getResampledSigmaBubble{};
template<> struct getResampledSigmaBubble<jackknifeDistributionD>{ static inline const NumericTensor<jackknifeDistributionD,1> &get(const ProjectedSigmaBubbleData &bubble_data){ return bubble_data.bubble_j; }  };
template<> struct getResampledSigmaBubble<doubleJackknifeDistributionD>{ static inline const NumericTensor<doubleJackknifeDistributionD,1> &get(const ProjectedSigmaBubbleData &bubble_data){ return bubble_data.bubble_dj; }  };


//Compute and store the raw amplitude data and data necessary for mix and vacuum subtractions
//Note: type index 1 is unused for K->sigma
class RawKtoSigmaData{
public:
  //Type data
  IndexedContainer<NumericTensor<rawDataDistributionD,3>, 3, 2> A0_alltK; //(type idx){q,tK,t}

  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_nobub; //{q,tK,t}
  
  //Mix data
  IndexedContainer<NumericTensor<rawDataDistributionD,2>, 2, 3> mix_alltK; //(mix idx){tK,t}

  NumericTensor<rawDataDistributionD,2> mix4_alltK_nobub;  

private:
  //Data dir requires *3* formats: type1/2,  type3, type4
  void getTypeData(IndexedContainer<type1234Data, 3, 2> &type_data, const int tsep_k_sigma, 
		   const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
		   const int traj_start, const int traj_inc, const int traj_lessthan,
		   const int Lt, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    //Read the data
    if(opt.load_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.load_data_checkpoint_stub << "_ktosigma_tsepksigma" << tsep_k_sigma << ".hdf5";
      std::cout << "Loading K->sigma checkpoint data for tsep_k_sigma = " << tsep_k_sigma << " from " << file.str() << std::endl;
      HDF5reader rd(file.str());
      for(int i=2;i<=4;i++) read(rd,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }else{
      for(int i=2;i<=4;i++) type_data(i) = readKtoSigmaType(i, traj_start, traj_inc, traj_lessthan, 
							    tsep_k_sigma, Lt, data_dir, data_file_fmt[i-2]);
    }

    //Write checkpoint if necessary
    if(opt.save_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.save_data_checkpoint_stub << "_ktosigma_tsepksigma" << tsep_k_sigma << ".hdf5";
      std::cout << "Saving checkpoint data for tsep_k_sigma = " << tsep_k_sigma << " to " << file.str() << std::endl;
      HDF5writer wr(file.str());
      for(int i=1;i<=4;i++) write(wr,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }
  }
 
public:
  RawKtoSigmaData(){}

  RawKtoSigmaData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		  const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
		  const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		  const int Lt, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    std::cout << "Reading K->sigma data with tsep_k_sigma=" << tsep_k_sigma << std::endl;

    IndexedContainer<type1234Data, 3, 2> type_data;
    getTypeData(type_data, tsep_k_sigma, data_dir, data_file_fmt, traj_start, traj_inc, traj_lessthan, Lt, opt);
    
    //Get the type1-4 and mix3/mix4 components of the data.
    //Note that we have not yet multiplied the type4/mix4 data by the pipi bubble, hence these  are just the K->vac component which are needed for computing alpha
    std::cout << "Computing raw data diagrams prior to multiplying by bubble\n";
    for(int i=2;i<=3;i++) A0_alltK(i) = computeKtoSigmaAmplitudeType<computeAmplitudeAlltKtensorControls>(i,type_data(i)); //[Qidx][tK][t]
    A0_type4_alltK_nobub = computeKtoSigmaAmplitudeType<computeAmplitudeAlltKtensorControls>(4,type_data(4)); //[Qidx][tK][t]
    
    mix_alltK(3) = NumericTensor<rawDataDistributionD,2>({Lt,Lt}); //[tK][t]
    mix4_alltK_nobub = NumericTensor<rawDataDistributionD,2>({Lt,Lt}); //[tK][t]
    for(int tK=0;tK<Lt;tK++)
      for(int t=0;t<Lt;t++){
	mix_alltK(3)({tK,t}) = type_data(3)(tK,t).mix();
	mix4_alltK_nobub({tK,t}) = type_data(4)(tK,t).mix();
      }
    
    printMem("Pre typedata free");
    
    //Data is no longer needed, so free it
    for(int i=2;i<=4;i++) type_data(i).freeData();
    
    printMem("Post typedata free");
    
    //Compute the type4/mix4 data with the bubble included
    std::cout << "Computing raw type4/mix4 data with bubble included" << std::endl;
    A0_alltK(4) = NumericTensor<rawDataDistributionD,3>({10,Lt,Lt});
    mix_alltK(4) = NumericTensor<rawDataDistributionD,2>({Lt,Lt});
    for(int tK=0;tK<Lt;tK++) for(int t=0;t<Lt;t++){
	int tB = (tK + tsep_k_sigma) % Lt;
	mix_alltK(4)({tK,t}) = mix4_alltK_nobub({tK,t})*bubble_data.bubble(&tB);
	for(int q=0;q<10;q++)
	  A0_alltK(4)({q,tK,t}) = A0_type4_alltK_nobub({q,tK,t})*bubble_data.bubble(&tB);
      }
    
    //Bin everything we are going to use henceforth
    std::cout << "Binning data\n";
    for(int i=2;i<=4;i++) A0_alltK(i) = bin(A0_alltK(i), bin_size);
    for(int i=3;i<=4;i++) mix_alltK(i) = bin(mix_alltK(i), bin_size);
    A0_type4_alltK_nobub = bin(A0_type4_alltK_nobub, bin_size);
    mix4_alltK_nobub = bin(mix4_alltK_nobub, bin_size);  
    
    std::cout << "Finished reading raw K->sigma data with tsep_k_sigma=" << tsep_k_sigma << std::endl;
  }
};


//------------------------------------------------------
//get_data
//-----------------------------------------------------

template<typename DistributionType, typename Resampler>
NumericTensor<DistributionType,1> computeQamplitude(const int q, const int tsep_k_sigma, const RawKtoSigmaData &raw, const ProjectedSigmaBubbleData &bubble_data, const int Lt, const std::string &descr, const Resampler &resampler){
  std::vector<int> nonzero_tK(Lt); for(int t=0;t<Lt;t++) nonzero_tK[t] = t;

  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
  computeAlphaAndVacuumSubtractions(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
				    raw.A0_type4_alltK_nobub, raw.mix4_alltK_nobub, getResampledSigmaBubble<DistributionType>::get(bubble_data),q, nonzero_tK,tsep_k_sigma,Lt,resampler);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 3, 2> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=2;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeData<DistributionType>(raw.A0_alltK(i), q, nonzero_tK, Lt, resampler); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = resampleAverageMixDiagram<DistributionType>(raw.mix_alltK(i), nonzero_tK, Lt, resampler);

  //Subtract the pseudoscalar operators and mix4 vacuum term
  std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under " << descr << "\n";
  A0_srcavg_r(3) = A0_srcavg_r(3).transform([&](int const* t, const DistributionType &from){ return DistributionType(from - alpha_r(t)*mix_srcavg_r(3)(t)); }); 
  A0_srcavg_r(4) = A0_srcavg_r(4).transform([&](int const* t, const DistributionType &from){
      return DistributionType(from - alpha_r(t)*( mix_srcavg_r(4)(t) - mix4_srcavg_vacsub_r(t) ) );
    }); 

  //Perform the type 4 vacuum subtraction
  std::cout << "Performing type-4 vacuum subtraction\n";
  A0_srcavg_r(4) = A0_srcavg_r(4) - A0_type4_srcavg_vacsub_r;

  //Get the full double-jackknife amplitude
  std::cout << "Computing full amplitudes\n";
  NumericTensor<DistributionType,1> A0_full_srcavg_r = A0_srcavg_r(2) + A0_srcavg_r(3) + A0_srcavg_r(4);

  return A0_full_srcavg_r;
}

//Read and prepare the data for a particular tsep_k_pi_idx
template<typename Resampler>
void getKtoSigmaData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		     const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		     const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
		     const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		     const int Lt, const Resampler &resampler, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

  std::cout << "Getting K->sigma data for tsep_k_sigma = " <<  tsep_k_sigma << std::endl;
  printMem("getData called");
  
  RawKtoSigmaData raw(tsep_k_sigma, bubble_data, data_dir, data_file_fmt, traj_start, traj_inc, traj_lessthan, bin_size, Lt, opt);

  for(int q=0;q<10;q++){
    std::cout << "Starting K->sigma Q" << q+1 << std::endl;

    printMem("Starting new K->sigma Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitude<doubleJackknifeDistributionD>(q, tsep_k_sigma, raw, bubble_data, Lt, "double jackknife", resampler);
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitude<jackknifeDistributionD>(q, tsep_k_sigma, raw, bubble_data, Lt, "single jackknife", resampler);

    //Insert data into output containers    
    for(int t=0;t<Lt;t++){
      A0_all_j[q].push_back(amplitudeDataCoord(t,tsep_k_sigma), A0_full_srcavg_j(&t));
      A0_all_dj[q].push_back(amplitudeDataCoord(t,tsep_k_sigma), A0_full_srcavg_dj(&t));
    }
  }
}


template<typename Resampler>
void getKtoSigmaData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		     const std::vector<int> &tsep_k_sigma,
		     const std::string &data_dir,  
		     const std::vector<std::string> &data_file_fmt, 
		     const std::string &bubble_file_fmt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj,
		     const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		     const int Lt, const Resampler &resampler, const readKtoPiPiAllDataOptions &opt = readKtoPiPiAllDataOptions()){

  if(opt.load_amplitude_data){
#ifdef HAVE_HDF5
    HDF5reader reader(opt.load_amplitude_data_file);
    read(reader, A0_all_j, "A0_all_j");
    read(reader, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    //Read the bubble data
    ProjectedSigmaBubbleData bubble_data(data_dir, bubble_file_fmt, traj_start, traj_inc, traj_lessthan, bin_size, Lt, bubble_quarkmom_proj, resampler, opt.read_opts);
  
    //Read and prepare the amplitude data for fitting
    scratch scratch_store(opt.use_scratch, opt.use_scratch_stub, opt.use_existing_scratch_files, tsep_k_sigma); //Setup scratch space if in use

    for(int tsep_k_sigma_idx=0;tsep_k_sigma_idx<tsep_k_sigma.size();tsep_k_sigma_idx++){
      if(scratch_store.doSkipLoad(tsep_k_sigma_idx)) continue;
      
      getKtoSigmaData(A0_all_j,A0_all_dj,tsep_k_sigma[tsep_k_sigma_idx],bubble_data,
		      data_dir, data_file_fmt, traj_start, traj_inc, traj_lessthan, bin_size,
		      Lt, resampler, opt.read_opts);

      scratch_store.writeScratch(A0_all_j, A0_all_dj, tsep_k_sigma_idx);
    }	

    scratch_store.reloadScratch(A0_all_j, A0_all_dj);      
  }

  if(opt.save_amplitude_data){
#ifdef HAVE_HDF5
    HDF5writer writer(opt.save_amplitude_data_file);
    write(writer, A0_all_j, "A0_all_j");
    write(writer, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }

}


void getKtoSigmaData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		     const std::vector<int> &tsep_k_sigma,
		     const std::string &data_dir,  
		     const std::vector<std::string> &data_file_fmt, 
		     const std::string &bubble_file_fmt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj,
		     const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		     const int Lt, const readKtoPiPiAllDataOptions &opt = readKtoPiPiAllDataOptions()){
  basic_resampler resampler;
  getKtoSigmaData(A0_all_j, A0_all_dj, tsep_k_sigma,
		  data_dir, 
		  data_file_fmt,
		  bubble_file_fmt, bubble_quarkmom_proj, 
		  traj_start, traj_inc, traj_lessthan, bin_size, 
		  Lt, resampler, opt);
}




CPSFIT_END_NAMESPACE





using namespace CPSfit;

int main(const int argc, const char* argv[]){
  printMem("Beginning of execution");
  
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    of << args;
    std::cout << "Wrote template argument file to template.args\n";
    return 0;
  }
  const std::string arg_file = argv[1];
  parse(args, arg_file);

  CMDline cmdline(argc,argv,2);

  const std::vector<std::pair<threeMomentum, double> > bubble_quarkmom_proj = { { {1,1,1}, 1./8 }, { {-1,-1,-1}, 1./8 },
										{ {-3,1,1}, 1./8 }, { {3,-1,-1}, 1./8 },
										{ {1,-3,1}, 1./8 }, { {-1,3,-1}, 1./8 },
										{ {1,1,-3}, 1./8 }, { {-1,-1,3}, 1./8 }       };
  readKtoPiPiAllDataOptions read_opt;
  read_opt.import(cmdline);

  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);
  getKtoSigmaData(A0_all_j, A0_all_dj, args.tsep_k_sigma, args.data_dir, args.type_file_fmt,  args.bubble_file_fmt, 
		  bubble_quarkmom_proj, args.traj_start, args.traj_inc, args.traj_lessthan, args.bin_size, args.Lt, read_opt);
  
  printMem("Prior to fitting");
  fitAndPlot(A0_all_j,A0_all_dj, args.Lt, args.tmin_k_op, args.tmin_op_sigma, args.fitfunc, args.correlated, cmdline.load_freeze_data, cmdline.freeze_data);
  
  std::cout << "Done" << std::endl;
  
  return 0;
}

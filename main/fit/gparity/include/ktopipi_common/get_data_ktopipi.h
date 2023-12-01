#ifndef _FIT_KTOPIPI_GETDATA_H
#define _FIT_KTOPIPI_GETDATA_H

#include<config.h>
#include<utils/macros.h>

#include "scratch.h"
#include "amplitude_data_ktopipi.h"
#include "computeQ_amplitude_opts.h"

CPSFIT_START_NAMESPACE

template<typename resampledDistributionType, typename BinResampler>
void computeAlphaAndVacuumSubtractions(NumericTensor<resampledDistributionType,1> &alpha,
				       NumericTensor<resampledDistributionType,1> &A0_type4_srcavg_vacsub,
				       NumericTensor<resampledDistributionType,1> &mix4_srcavg_vacsub,
				       const NumericTensor<rawDataDistributionD,3> &A0_type4_nobub_alltK,
				       const NumericTensor<rawDataDistributionD,2> &mix4_nobub_alltK,
				       const NumericTensor<resampledDistributionType,1> &bubble_rs,
				       const int q,
				       const std::vector<int> &type4_nonzerotK,
				       const int tsep_k_pi,
				       const int Lt, const BinResampler &bin_resampler){
  //Compute mix4 double-jackknife and tK average
  resampledDistributionType zro = bubble_rs({0}); zeroit(zro);
  
#pragma omp parallel for
  for(int t=0;t<=tsep_k_pi;t++){
    mix4_srcavg_vacsub(&t) = zro;
    A0_type4_srcavg_vacsub(&t) = zro;

    resampledDistributionType mix4_nobub_srcavg = zro;
    resampledDistributionType A0_type4_nobub_srcavg = zro;
    resampledDistributionType mix4_nobub_rs, mix4_vacsub_rs, A0_type4_nobub_rs, A0_type4_vacsub_rs;

    for(int tK_idx=0;tK_idx<type4_nonzerotK.size();tK_idx++){
      const int tK = type4_nonzerotK[tK_idx];
      const int tB = (tK + tsep_k_pi) % Lt;      

      bin_resampler.binResample(mix4_nobub_rs, mix4_nobub_alltK({tK_idx,t}) );
      mix4_vacsub_rs = mix4_nobub_rs * bubble_rs(&tB);

      bin_resampler.binResample(A0_type4_nobub_rs, A0_type4_nobub_alltK({q,tK_idx,t}) );
      A0_type4_vacsub_rs = A0_type4_nobub_rs * bubble_rs(&tB);

      mix4_srcavg_vacsub(&t) = mix4_srcavg_vacsub(&t) + mix4_vacsub_rs;
      mix4_nobub_srcavg = mix4_nobub_srcavg + mix4_nobub_rs;

      A0_type4_srcavg_vacsub(&t) = A0_type4_srcavg_vacsub(&t) + A0_type4_vacsub_rs;
      A0_type4_nobub_srcavg = A0_type4_nobub_srcavg + A0_type4_nobub_rs;
    }
    double n(type4_nonzerotK.size());
    
    mix4_srcavg_vacsub(&t) = mix4_srcavg_vacsub(&t)/n;
    A0_type4_srcavg_vacsub(&t) = A0_type4_srcavg_vacsub(&t)/n;
    mix4_nobub_srcavg = mix4_nobub_srcavg/n;
    A0_type4_nobub_srcavg = A0_type4_nobub_srcavg/n;

    alpha(&t) = A0_type4_nobub_srcavg/mix4_nobub_srcavg;
  }
}


template<typename DistributionType, typename BinResampler>
NumericTensor<DistributionType,1> binResampleAverageTypeData(const NumericTensor<rawDataDistributionD,3> &typedata_alltK,
							     const int q,
							     const BinResampler &bin_resampler){
  const int ntK = typedata_alltK.size(1);
  const int nt = typedata_alltK.size(2);

  NumericTensor<DistributionType,1> out({nt}); //[t]
  for(int t=0;t<nt;t++)
    binResampleAverage(out(&t), bin_resampler, [&](const int tK_idx)->rawDataDistributionD const& { return typedata_alltK({q,tK_idx,t}); }, ntK);
  return out;
}

template<typename DistributionType, typename BinResampler>
NumericTensor<DistributionType,1> binResampleAverageMixDiagram(const NumericTensor<rawDataDistributionD,2> &mixdata_alltK,
							       const BinResampler &bin_resampler){
  const int ntK = mixdata_alltK.size(0);
  const int nt = mixdata_alltK.size(1);

  NumericTensor<DistributionType,1> out({nt}); //[t]
  for(int t=0;t<nt;t++)
    binResampleAverage(out(&t), bin_resampler, [&](const int tK_idx)->rawDataDistributionD const& { return mixdata_alltK({tK_idx,t}); }, ntK);
  return out;
}

template<typename DistributionType, typename BinResampler>
NumericTensor<DistributionType,1> computeQamplitude(const int q, const int tsep_k_pi, const RawKtoPiPiData &raw, const NumericTensor<DistributionType,1> &resampled_bubble, const int Lt, const std::string &descr, const BinResampler &bin_resampler, const computeQamplitudeOpts &opts = computeQamplitudeOpts()){
  std::cout << "Starting generation of resampled amplitude for K->pipi Q" << q+1 << " and t_sep(K->pi)=" << tsep_k_pi << std::endl;

  const int nt = tsep_k_pi + 1; //only compute on  0<=t<=tsep_k_pi

  timer time;

  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions" << std::endl;
  time.start();
  NumericTensor<DistributionType,1> alpha_r({nt}), A0_type4_srcavg_vacsub_r({nt}), mix4_srcavg_vacsub_r({nt}); //[t]
  computeAlphaAndVacuumSubtractions(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
				    raw.A0_type4_alltK_nobub, raw.mix4_alltK_nobub, resampled_bubble,q, raw.nonzerotK(4),tsep_k_pi,Lt, bin_resampler);
  time.stop();
  std::cout << "Alpha and vacuum subtractions: " << time.elapsed()/1.e9 << "s" << std::endl;
  

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages of binned, resampled data" << std::endl;
  time.start();
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=1;i<=4;i++) A0_srcavg_r(i) = binResampleAverageTypeData<DistributionType>(raw.A0_alltK(i), q, bin_resampler); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = binResampleAverageMixDiagram<DistributionType>(raw.mix_alltK(i), bin_resampler);
  time.stop();
  std::cout << "tK averages of binned, resampled data: " << time.elapsed()/1.e9 << "s" << std::endl;

  time.start();
  //Subtract the pseudoscalar operators and mix4 vacuum term
  std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under " << descr << std::endl;
  A0_srcavg_r(3) = A0_srcavg_r(3).transform([&](int const* t, const DistributionType &from){ 
      return DistributionType(from - opts.alpha_scale*alpha_r(t)*mix_srcavg_r(3)(t)); }); 
  A0_srcavg_r(4) = A0_srcavg_r(4).transform([&](int const* t, const DistributionType &from){
      DistributionType mix4_t = mix_srcavg_r(4)(t);
      if(opts.do_vacuum_subtraction) mix4_t = mix4_t - mix4_srcavg_vacsub_r(t);
      return DistributionType(from - opts.alpha_scale*alpha_r(t)*mix4_t);

      //return DistributionType(from - opts.alpha_scale*alpha_r(t)*( mix_srcavg_r(4)(t) - mix4_srcavg_vacsub_r(t) ) );
    }); 

  //Perform the type 4 vacuum subtraction
  if(opts.do_vacuum_subtraction){
    std::cout << "Performing type-4 vacuum subtraction" << std::endl;
    A0_srcavg_r(4) = A0_srcavg_r(4) - A0_type4_srcavg_vacsub_r;
  }

  //Get the full resampled amplitude
  std::cout << "Computing full amplitudes" << std::endl;
  NumericTensor<DistributionType,1> A0_full_srcavg_r = A0_srcavg_r(1) + A0_srcavg_r(2) + A0_srcavg_r(3) + A0_srcavg_r(4);
  time.stop();
  std::cout << "Remainder of computation: " << time.elapsed()/1.e9 << "s" << std::endl;

  return A0_full_srcavg_r;
}
template<typename DistributionType, typename BinResampler>
NumericTensor<DistributionType,1> computeQamplitude(const int q, const int tsep_k_pi, const RawKtoPiPiData &raw, const ProjectedBubbleData &bubble_data, 
						    const int Lt, const std::string &descr, const BinResampler &bin_resampler, 
						    const computeQamplitudeOpts &opts = computeQamplitudeOpts()){
  return computeQamplitude(q, tsep_k_pi, raw, bubble_data.binResample<DistributionType>(bin_resampler), Lt, descr, bin_resampler,opts);
}


//Read and prepare the data for a particular tsep_k_pi_idx
template<typename BinResampler>
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
	     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
	     const int tsep_k_pi, const ProjectedBubbleData &bubble_data, 
	     const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
	     const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
	     const int traj_start, const int traj_inc, const int traj_lessthan, 
	     const int Lt, const int tsep_pipi, 
	     const BinResampler &bin_resampler, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

  std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
  printMem("getData called");
  
  RawKtoPiPiData raw(tsep_k_pi, bubble_data, data_dir, data_file_fmt, type1_pimom_proj, traj_start, traj_inc, traj_lessthan, 
		     Lt, tsep_pipi, opt);
  
  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitude<doubleJackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, Lt, 
														      "double jackknife", bin_resampler);
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitude<jackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, Lt, 
													 "single jackknife", bin_resampler);

    //Insert data into output containers    
    for(int t=0;t<=tsep_k_pi;t++){
      A0_all_j[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_j(&t));
      A0_all_dj[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_dj(&t));
    }
  }
}


struct readKtoPiPiAllDataOptions{
  bool load_amplitude_data;
  std::string load_amplitude_data_file;

  bool save_amplitude_data;
  std::string save_amplitude_data_file;

  bool use_scratch;
  std::string use_scratch_stub;
  bool use_existing_scratch_files;

  readKtoPiPiDataOptions read_opts;

  readKtoPiPiAllDataOptions(): load_amplitude_data(false), save_amplitude_data(false), use_scratch(false){}

  template<typename T>
  void import(const T &from){
#define C(X) X = from.X
    C(load_amplitude_data);
    C(load_amplitude_data_file);
    C(save_amplitude_data);
    C(save_amplitude_data_file);
    C(use_scratch);
    C(use_scratch_stub);
    C(use_existing_scratch_files);
#undef C
    read_opts.import(from);
  }
};


template<typename BinResampler>
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
	     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
	     const std::vector<int> &tsep_k_pi,
	     const std::string &data_dir,  
	     const std::vector<std::string> &data_file_fmt, const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
	     const std::string &bubble_file_fmt, const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,
	     const int traj_start, const int traj_inc, const int traj_lessthan,
	     const int Lt, const int tsep_pipi, 
	     const BinResampler &bin_resampler, const readKtoPiPiAllDataOptions &opt = readKtoPiPiAllDataOptions()){

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
    ProjectedBubbleData bubble_data(data_dir, bubble_file_fmt, traj_start, traj_inc, traj_lessthan, Lt, tsep_pipi, bubble_pimom_proj, opt.read_opts);
  
    //Read and prepare the amplitude data for fitting
    scratch scratch_store(opt.use_scratch, opt.use_scratch_stub, opt.use_existing_scratch_files, tsep_k_pi); //Setup scratch space if in use

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<tsep_k_pi.size();tsep_k_pi_idx++){
      if(scratch_store.doSkipLoad(tsep_k_pi_idx)) continue;
      
      getData(A0_all_j,A0_all_dj,tsep_k_pi[tsep_k_pi_idx],bubble_data,
	      data_dir, data_file_fmt, type1_pimom_proj, traj_start, traj_inc, traj_lessthan,
	      Lt, tsep_pipi, bin_resampler, opt.read_opts);

      scratch_store.writeScratch(A0_all_j, A0_all_dj, tsep_k_pi_idx);
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

void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
	     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
	     const std::vector<int> &tsep_k_pi,
	     const std::string &data_dir,  
	     const std::vector<std::string> &data_file_fmt, const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
	     const std::string &bubble_file_fmt, const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,
	     const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size,
	     const int Lt, const int tsep_pipi, 
	     const readKtoPiPiAllDataOptions &opt = readKtoPiPiAllDataOptions()){
  basicBinResampler bin_resampler(bin_size);
  getData(A0_all_j, A0_all_dj, tsep_k_pi,
	  data_dir, 
	  data_file_fmt, type1_pimom_proj, 
	  bubble_file_fmt, bubble_pimom_proj, 
	  traj_start, traj_inc, traj_lessthan, 
	  Lt, tsep_pipi, bin_resampler, opt);
}

void checkpointRawOnly(const std::vector<int> &tsep_k_pi,
		       const std::string &data_dir,  
		       const std::vector<std::string> &data_file_fmt,
		       const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
		       const std::string &bubble_file_fmt,
		       const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,
		       const int traj_start, const int traj_inc, const int traj_lessthan, 
		       const int Lt, const int tsep_pipi, 
		       const readKtoPiPiDataOptions &opt){
  if(!opt.save_data_checkpoint) error_exit(std::cout << "checkpointRawOnly expect opt.save_data_checkpoint to be true\n");
  ProjectedBubbleData bubble_data(data_dir, bubble_file_fmt, traj_start, traj_inc, traj_lessthan, Lt, tsep_pipi, bubble_pimom_proj, opt);

  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<tsep_k_pi.size();tsep_k_pi_idx++){
    int tsep_k_pi_ = tsep_k_pi[tsep_k_pi_idx];
    std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi_ << std::endl;
    RawKtoPiPiData raw(tsep_k_pi_, bubble_data, data_dir, data_file_fmt, type1_pimom_proj, traj_start, traj_inc, traj_lessthan, 
		       Lt, tsep_pipi, opt);
  }
}

CPSFIT_END_NAMESPACE

#endif

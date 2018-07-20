#ifndef _FIT_KTOPIPI_GETDATA_H
#define _FIT_KTOPIPI_GETDATA_H

#include<config.h>
#include<utils/macros.h>

#include "scratch.h"
#include "amplitude_data.h"

CPSFIT_START_NAMESPACE

template<typename resampledDistributionType, typename Resampler>
void computeAlphaAndVacuumSubtractions(NumericTensor<resampledDistributionType,1> &alpha,
				       NumericTensor<resampledDistributionType,1> &A0_type4_srcavg_vacsub,
				       NumericTensor<resampledDistributionType,1> &mix4_srcavg_vacsub,
				       const NumericTensor<rawDataDistributionD,3> &A0_type4_nobub_alltK,
				       const NumericTensor<rawDataDistributionD,2> &mix4_nobub_alltK,
				       const NumericTensor<resampledDistributionType,1> &bubble_rs,
				       const int q,
				       const std::vector<int> &type4_nonzerotK,
				       const int tsep_k_pi,
				       const int Lt, const Resampler &resampler){
  //Compute mix4 double-jackknife and tK average
  resampledDistributionType zro = bubble_rs({0}); zeroit(zro);
  
#pragma omp parallel for
  for(int t=0;t<Lt;t++){
    mix4_srcavg_vacsub(&t) = zro;
    A0_type4_srcavg_vacsub(&t) = zro;

    resampledDistributionType mix4_nobub_srcavg = zro;
    resampledDistributionType A0_type4_nobub_srcavg = zro;
    
    for(int ii=0;ii<type4_nonzerotK.size();ii++){
      const int tK = type4_nonzerotK[ii];
      const int tB = (tK + tsep_k_pi) % Lt;      

      resampledDistributionType mix4_nobub_rs; resampler.resample(mix4_nobub_rs, mix4_nobub_alltK({tK,t}) );
      resampledDistributionType mix4_vacsub_rs = mix4_nobub_rs * bubble_rs(&tB);

      resampledDistributionType A0_type4_nobub_rs;  resampler.resample(A0_type4_nobub_rs, A0_type4_nobub_alltK({q,tK,t}) );
      resampledDistributionType A0_type4_vacsub_rs = A0_type4_nobub_rs * bubble_rs(&tB);

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


template<typename DistributionType, typename Resampler>
NumericTensor<DistributionType,1> resampleAverageTypeData(const NumericTensor<rawDataDistributionD,3> &typedata_alltK,
							  const int q,
							  const std::vector<int> &typedata_nonzerotK,
							  const int Lt, const Resampler &resampler){
  NumericTensor<DistributionType,1> out({Lt}); //[t]
  for(int t=0;t<Lt;t++)
    resampleAverage(out(&t), resampler, [&](const int i){ return typedata_alltK({q,typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
  return out;
}

template<typename DistributionType, typename Resampler>
NumericTensor<DistributionType,1> resampleAverageMixDiagram(const NumericTensor<rawDataDistributionD,2> &mixdata_alltK,
							    const std::vector<int> &mixdata_nonzerotK,
							    const int Lt, const Resampler &resampler){
  NumericTensor<DistributionType,1> out({Lt}); //[t]
  for(int t=0;t<Lt;t++)
    resampleAverage(out(&t), resampler, [&](const int i){ return mixdata_alltK({mixdata_nonzerotK[i],t}); }, mixdata_nonzerotK.size());
  return out;
}



template<typename DistributionType, typename Resampler>
NumericTensor<DistributionType,1> computeQamplitude(const int q, const int tsep_k_pi, const RawKtoPiPiData &raw, const BubbleData &bubble_data, const int Lt, const std::string &descr, const Resampler &resampler){
  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
  computeAlphaAndVacuumSubtractions(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
				    raw.A0_type4_alltK_nobub, raw.mix4_alltK_nobub, getResampledBubble<DistributionType>::get(bubble_data),q, raw.nonzerotK(4),tsep_k_pi,Lt,resampler);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeData<DistributionType>(raw.A0_alltK(i), q, raw.nonzerotK(i), Lt, resampler); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = resampleAverageMixDiagram<DistributionType>(raw.mix_alltK(i), raw.nonzerotK(i), Lt, resampler);

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
  NumericTensor<DistributionType,1> A0_full_srcavg_r = A0_srcavg_r(1) + A0_srcavg_r(2) + A0_srcavg_r(3) + A0_srcavg_r(4);

  return A0_full_srcavg_r;
}

//Read and prepare the data for a particular tsep_k_pi_idx
template<typename Resampler>
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
	     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
	     const int tsep_k_pi, const BubbleData &bubble_data, 
	     const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
	     const int Lt, const int tsep_pipi, 
	     const Resampler &resampler, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

  std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
  printMem("getData called");
  
  RawKtoPiPiData raw(tsep_k_pi, bubble_data, data_dir, traj_start, traj_inc, traj_lessthan, bin_size, 
		     Lt, tsep_pipi, opt);
  
  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitude<doubleJackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, Lt, "double jackknife", resampler);
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitude<jackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, Lt, "single jackknife", resampler);

    //Insert data into output containers    
    for(int t=0;t<Lt;t++){
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


template<typename Resampler>
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
	     std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
	     const std::vector<int> &tsep_k_pi,
	     const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
	     const int Lt, const int tsep_pipi, 
	     const Resampler &resampler, const readKtoPiPiAllDataOptions &opt = readKtoPiPiAllDataOptions()){

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
    BubbleData bubble_data(data_dir, traj_start, traj_inc, traj_lessthan, bin_size,
			   Lt, tsep_pipi, resampler, opt.read_opts);
    
    //Read and prepare the amplitude data for fitting
    scratch scratch_store(opt.use_scratch, opt.use_scratch_stub, opt.use_existing_scratch_files, tsep_k_pi); //Setup scratch space if in use

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<tsep_k_pi.size();tsep_k_pi_idx++){
      if(scratch_store.doSkipLoad(tsep_k_pi_idx)) continue;
      
      getData(A0_all_j,A0_all_dj,tsep_k_pi[tsep_k_pi_idx],bubble_data,
	      data_dir, traj_start, traj_inc, traj_lessthan, bin_size,
	      Lt, tsep_pipi, resampler, opt.read_opts);

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
	     const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
	     const int Lt, const int tsep_pipi, 
	     const readKtoPiPiAllDataOptions &opt = readKtoPiPiAllDataOptions()){
  basic_resampler resampler;
  getData(A0_all_j, A0_all_dj, tsep_k_pi,
	  data_dir, traj_start, traj_inc, traj_lessthan, bin_size, 
	  Lt, tsep_pipi, resampler, opt);
}


void checkpointRawOnly(const std::vector<int> &tsep_k_pi,
		       const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		       const int Lt, const int tsep_pipi, 
		       const readKtoPiPiDataOptions &opt){
  if(!opt.save_data_checkpoint) error_exit(std::cout << "checkpointRawOnly expect opt.save_data_checkpoint to be true\n");
  basic_resampler resampler;
  BubbleData bubble_data(data_dir, traj_start, traj_inc, traj_lessthan, bin_size,
			 Lt, tsep_pipi, resampler, opt);

  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<tsep_k_pi.size();tsep_k_pi_idx++){
    int tsep_k_pi_ = tsep_k_pi[tsep_k_pi_idx];
    std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi_ << std::endl;
    RawKtoPiPiData raw(tsep_k_pi_, bubble_data, data_dir, traj_start, traj_inc, traj_lessthan, bin_size, 
		       Lt, tsep_pipi, opt);
  }
}

CPSFIT_END_NAMESPACE

#endif

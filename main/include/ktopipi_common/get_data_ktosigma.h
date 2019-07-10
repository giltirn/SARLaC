#ifndef _FIT_KTOSIGMA_GETDATA_H
#define _FIT_KTOSIGMA_GETDATA_H

#include<config.h>
#include<utils/macros.h>

#include "scratch.h"
#include "amplitude_data_ktosigma.h"

CPSFIT_START_NAMESPACE

template<typename DistributionType, typename Resampler>
NumericTensor<DistributionType,1> computeQamplitude(const int q, const int tsep_k_sigma, const RawKtoSigmaData &raw, const ProjectedSigmaBubbleData &bubble_data, const int Lt, const std::string &descr, const int bin_size, const Resampler &resampler){
  int nt = tsep_k_sigma + 1; //only compute on  0<=t<=tsep_k_sigma

  std::vector<int> nonzero_tK(Lt); for(int t=0;t<Lt;t++) nonzero_tK[t] = t;

  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions" << std::endl;
  NumericTensor<DistributionType,1> alpha_r({nt}), A0_type4_srcavg_vacsub_r({nt}), mix4_srcavg_vacsub_r({nt}); //[t]
  computeAlphaAndVacuumSubtractions(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
				    raw.A0_type4_alltK_nobub, raw.mix4_alltK_nobub, bubble_data.binResample<DistributionType>(bin_size, resampler),q, 
				    nonzero_tK,tsep_k_sigma,Lt,bin_size,resampler);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams" << std::endl;
  IndexedContainer<NumericTensor<DistributionType,1>, 3, 2> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=2;i<=4;i++) A0_srcavg_r(i) = binResampleAverageTypeData<DistributionType>(raw.A0_alltK(i), q, bin_size, resampler); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = binResampleAverageMixDiagram<DistributionType>(raw.mix_alltK(i), bin_size, resampler);

  //Subtract the pseudoscalar operators and mix4 vacuum term
  std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under " << descr << std::endl;
  //NOTE: <sigma|P|K> = -|C|*mix3    where the |C| cancels with a corresponding factor absorbed into the definition of alpha here 
  //      (cf https://rbc.phys.columbia.edu/rbc_ukqcd/individual_postings/ckelly/Gparity/operator_subtraction.pdf pg 4)
  A0_srcavg_r(3) = A0_srcavg_r(3).transform([&](int const* t, const DistributionType &from){ return DistributionType(from + alpha_r(t)*mix_srcavg_r(3)(t)); }); 
  A0_srcavg_r(4) = A0_srcavg_r(4).transform([&](int const* t, const DistributionType &from){
      return DistributionType(from - alpha_r(t)*( mix_srcavg_r(4)(t) - mix4_srcavg_vacsub_r(t) ) );
    }); 

  //Perform the type 4 vacuum subtraction
  std::cout << "Performing type-4 vacuum subtraction" << std::endl;
  A0_srcavg_r(4) = A0_srcavg_r(4) - A0_type4_srcavg_vacsub_r;

  //Get the full double-jackknife amplitude
  std::cout << "Computing full amplitudes" << std::endl;
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
  
  RawKtoSigmaData raw(tsep_k_sigma, bubble_data, data_dir, data_file_fmt, traj_start, traj_inc, traj_lessthan, Lt, opt);

  for(int q=0;q<10;q++){
    std::cout << "Starting K->sigma Q" << q+1 << std::endl;

    printMem("Starting new K->sigma Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitude<doubleJackknifeDistributionD>(q, tsep_k_sigma, raw, bubble_data, Lt, "double jackknife", bin_size, resampler);
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitude<jackknifeDistributionD>(q, tsep_k_sigma, raw, bubble_data, Lt, "single jackknife", bin_size, resampler);

    //Insert data into output containers    
    for(int t=0;t<=tsep_k_sigma;t++){
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
    ProjectedSigmaBubbleData bubble_data(data_dir, bubble_file_fmt, traj_start, traj_inc, traj_lessthan, Lt, bubble_quarkmom_proj, opt.read_opts);
  
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

#endif

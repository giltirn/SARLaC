#ifndef FIT_KTOPIPI_GPARITY_SAMPLEAMA_GETDATA_H__
#define FIT_KTOPIPI_GPARITY_SAMPLEAMA_GETDATA_H__

#include<config.h>
#include<utils/macros.h>

#include <ktopipi_common/fitfunc.h>
#include <ktopipi_common/scratch.h>
#include <ktopipi_sampleAMA_common/alpha_vac_sub.h>
#include <ktopipi_sampleAMA_common/resample_average_typedata.h>

CPSFIT_START_NAMESPACE


#define ALPHA_AND_VAC_SUB_SEPARATE
//#define VAC_SUB_SEPARATE

#ifdef VAC_SUB_SEPARATE

template<typename DistributionType>
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const int Lt, const allRawData &raw, const allBubbleData &bubble_data, 
							     const sampleAMA_resamplers &resamplers, const std::string &descr){
  const int sloppy_S = 0;
  const int sloppy_C = 1;
  const int exact_C = 2;

  const int Lt = inputs.args.Lt;
  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";

  NumericTensor<DistributionType,1> dummy({Lt});

  //Get alpha from corrected Green's functions
  NumericTensor<DistributionType,1> alpha_r({Lt}); //[t]
  computeAlphaAndVacuumSubtractionsSampleAMACorrected(alpha_r, dummy, dummy,
						      raw, getResampledBubbleSampleAMA<DistributionType>::get(bubble_data),q,tsep_k_pi,Lt,resamplers);

  //Get vaccum subtractions separately
  std::array<NumericTensor<DistributionType,1>, 3> A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r; //[ens-status][t]
  
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(dummy, A0_type4_srcavg_vacsub_r[sloppy_S], mix4_srcavg_vacsub_r[sloppy_S],
						       raw, bubble_data,q,tsep_k_pi, Lt, Sloppy, resamplers.resampler_S);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(dummy, A0_type4_srcavg_vacsub_r[sloppy_C], mix4_srcavg_vacsub_r[sloppy_C],
						       raw, bubble_data,q,tsep_k_pi, Lt, Sloppy, resamplers.resample_C);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(dummy, A0_type4_srcavg_vacsub_r[exact_C], mix4_srcavg_vacsub_r[exact_C],
						       raw, bubble_data,q,tsep_k_pi, Lt, Exact, resamplers.resample_C);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  std::array< IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> , 3> mix_srcavg_r; //[ens-status][t]

  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,Lt,resamplers); //[t]  (does full correction)
  for(int i=3;i<=4;i++){
    mix_srcavg_r[sloppy_S](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,Sloppy,raw,Lt,resamplers.resampler_S);
    mix_srcavg_r[sloppy_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,Sloppy,raw,Lt,resamplers.resampler_C);
    mix_srcavg_r[exact_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,Exact,raw,Lt,resamplers.resampler_C);
  }

  //Subtract the pseudoscalar operators and mix4 vacuum term
  std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under " << descr << "\n";
  int status_sgn[3] = {1,-1,1}; //  (sloppy_S - sub_sloppy_S) - (sloppy_C - sub_sloppy_C) + (exact_C - sub_exact_C) 
  for(int status = 0 ; status < 3; status++){
    double sgn = status_sgn[status];
    A0_srcavg_r(3) = A0_srcavg_r(3).transform([&](int const* t, const DistributionType &from){ return DistributionType(from - sgn * alpha_r(t)*mix_srcavg_r[status](3)(t)); }); 
    
    A0_srcavg_r(4) = A0_srcavg_r(4).transform([&](int const* t, const DistributionType &from){
	return DistributionType(from - sgn * alpha_r(t)*( mix_srcavg_r[status](4)(t) - mix4_srcavg_vacsub_r[status](t) ) );
      }); 
  }

  //Perform the type 4 vacuum subtraction
  std::cout << "Performing type-4 vacuum subtraction\n";
  for(int status = 0 ; status < 3; status++){
    double sgn = status_sgn[status];
    A0_srcavg_r(4) = A0_srcavg_r(4) - sgn * A0_type4_srcavg_vacsub_r[status];
  }

  //Get the full double-jackknife amplitude
  std::cout << "Computing full amplitudes\n";
  NumericTensor<DistributionType,1> A0_full_srcavg_r = A0_srcavg_r(1) + A0_srcavg_r(2) + A0_srcavg_r(3) + A0_srcavg_r(4);

  return A0_full_srcavg_r;
}



#elif defined(ALPHA_AND_VAC_SUB_SEPARATE)

template<typename DistributionType>
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const int Lt, const allRawData &raw, const allBubbleData &bubble_data, 
							     const sampleAMA_resamplers &resamplers, const std::string &descr){
  const int sloppy_S = 0;
  const int sloppy_C = 1;
  const int exact_C = 2;

  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  std::array<NumericTensor<DistributionType,1>, 3> alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r; //[ens-status][t]
  
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(alpha_r[sloppy_S], A0_type4_srcavg_vacsub_r[sloppy_S], mix4_srcavg_vacsub_r[sloppy_S],
					     raw, bubble_data,q,tsep_k_pi, Lt,Sloppy, resamplers.resampler_S);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(alpha_r[sloppy_C], A0_type4_srcavg_vacsub_r[sloppy_C], mix4_srcavg_vacsub_r[sloppy_C],
					     raw, bubble_data,q,tsep_k_pi, Lt,Sloppy, resamplers.resampler_C);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(alpha_r[exact_C], A0_type4_srcavg_vacsub_r[exact_C], mix4_srcavg_vacsub_r[exact_C],
					     raw, bubble_data,q,tsep_k_pi, Lt,Exact, resamplers.resampler_C);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  std::array< IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> , 3> mix_srcavg_r; //[ens-status][t]

  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,Lt,resamplers); //[t]  (does full correction)
  for(int i=3;i<=4;i++){
    mix_srcavg_r[sloppy_S](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,Sloppy,raw,Lt,resamplers.resampler_S);
    mix_srcavg_r[sloppy_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,Sloppy,raw,Lt,resamplers.resampler_C);
    mix_srcavg_r[exact_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,Exact,raw,Lt,resamplers.resampler_C);
  }

  //Subtract the pseudoscalar operators and mix4 vacuum term
  std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under " << descr << "\n";
  int status_sgn[3] = {1,-1,1}; //  (sloppy_S - sub_sloppy_S) - (sloppy_C - sub_sloppy_C) + (exact_C - sub_exact_C) 
  for(int status = 0 ; status < 3; status++){
    double sgn = status_sgn[status];
    A0_srcavg_r(3) = A0_srcavg_r(3).transform([&](int const* t, const DistributionType &from){ return DistributionType(from - sgn * alpha_r[status](t)*mix_srcavg_r[status](3)(t)); }); 
    
    A0_srcavg_r(4) = A0_srcavg_r(4).transform([&](int const* t, const DistributionType &from){
	return DistributionType(from - sgn * alpha_r[status](t)*( mix_srcavg_r[status](4)(t) - mix4_srcavg_vacsub_r[status](t) ) );
      }); 
  }

  //Perform the type 4 vacuum subtraction
  std::cout << "Performing type-4 vacuum subtraction\n";
  for(int status = 0 ; status < 3; status++){
    double sgn = status_sgn[status];
    A0_srcavg_r(4) = A0_srcavg_r(4) - sgn * A0_type4_srcavg_vacsub_r[status];
  }

  //Get the full double-jackknife amplitude
  std::cout << "Computing full amplitudes\n";
  NumericTensor<DistributionType,1> A0_full_srcavg_r = A0_srcavg_r(1) + A0_srcavg_r(2) + A0_srcavg_r(3) + A0_srcavg_r(4);

  return A0_full_srcavg_r;
}


#else //Correct Green's functions prior to computing alpha and vacuum subtractions


template<typename DistributionType>
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const int Lt, const allRawData &raw, const allBubbleData &bubble_data, 
							     const sampleAMA_resamplers &resamplers, const std::string &descr){
  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
  computeAlphaAndVacuumSubtractionsSampleAMACorrected(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
						      raw, getResampledBubbleSampleAMA<DistributionType>::get(bubble_data),q,tsep_k_pi,Lt,resamplers);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,Lt,resamplers); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,raw,Lt,resamplers);

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

#endif


  
//Read and prepare the data for a particular tsep_k_pi_idx
void getDataSampleAMA(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		      std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		      const int tsep_k_pi, const allBubbleData &bubble_data,
		      const std::vector<std::string> &data_file_fmt_sloppy,  const std::vector<std::string> &data_file_fmt_exact,
		      const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
		      const std::string &data_dir_S, const int traj_start_S, const int traj_lessthan_S,
		      const std::string &data_dir_C, const int traj_start_C, const int traj_lessthan_C,
		      const int traj_inc, const int bin_size,
		      const int Lt, const int tsep_pipi, const sampleAMA_resamplers &resamplers, const readKtoPiPiDataSampleAMAoptions &opt = readKtoPiPiDataSampleAMAoptions()){
  std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
  printMem("getData called");
  
  allRawData raw(bubble_data, tsep_k_pi, 
		 data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		 data_dir_S, traj_start_S, traj_lessthan_S, 
		 data_dir_C, traj_start_C, traj_lessthan_C, 
		 traj_inc, bin_size,
		 Lt,tsep_pipi,opt);
  
  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitudeSampleAMA<doubleJackknifeDistributionD>(q, tsep_k_pi, Lt, raw, bubble_data, resamplers, "double jackknife");
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitudeSampleAMA<jackknifeDistributionD>(q, tsep_k_pi, Lt, raw, bubble_data, resamplers, "single jackknife");

    //Insert data into output containers    
    for(int t=0;t<Lt;t++){
      A0_all_j[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_j(&t));
      A0_all_dj[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_dj(&t));
    }
  }
}


struct readKtoPiPiAllDataSampleAMAoptions{
  bool load_amplitude_data;
  std::string load_amplitude_data_file;

  bool save_amplitude_data;
  std::string save_amplitude_data_file;

  bool use_scratch;
  std::string use_scratch_stub;
  bool use_existing_scratch_files;

  readKtoPiPiDataSampleAMAoptions read_opts;

  readKtoPiPiAllDataSampleAMAoptions(): load_amplitude_data(false), save_amplitude_data(false), use_scratch(false){}

  template<typename T>
  void importGlobalOptions(const T &from){
#define C(X) X = from.X
    C(load_amplitude_data);
    C(load_amplitude_data_file);
    C(save_amplitude_data);
    C(save_amplitude_data_file);
    C(use_scratch);
    C(use_scratch_stub);
    C(use_existing_scratch_files);
#undef C
  }
};


void getDataSampleAMA(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		      std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		      const std::vector<int> &tsep_k_pi,

		      const std::string &bubble_file_fmt_sloppy, const std::string &bubble_file_fmt_exact,
		      const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,

		      const std::vector<std::string> &data_file_fmt_sloppy,  const std::vector<std::string> &data_file_fmt_exact,
		      const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,

		      const std::string &data_dir_S, const int traj_start_S, const int traj_lessthan_S,
		      const std::string &data_dir_C, const int traj_start_C, const int traj_lessthan_C,
		      const int traj_inc, const int bin_size, 
		      const int Lt, const int tsep_pipi, 
		      const readKtoPiPiAllDataSampleAMAoptions &opt = readKtoPiPiAllDataSampleAMAoptions()){

  if(opt.load_amplitude_data){
#ifdef HAVE_HDF5
    HDF5reader reader(opt.load_amplitude_data_file);
    read(reader, A0_all_j, "A0_all_j");
    read(reader, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    sampleAMA_resamplers resamplers(traj_start_S, traj_lessthan_S,
				    traj_start_C, traj_lessthan_C,
				    traj_inc, bin_size);
  
    //Read the bubble data
    allBubbleData bubble_data(bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
			      data_dir_S, traj_start_S, traj_lessthan_S,
			      data_dir_C, traj_start_C, traj_lessthan_C,
			      traj_inc, bin_size, Lt, tsep_pipi, resamplers, opt.read_opts);

    scratch scratch_store(opt.use_scratch, opt.use_scratch_stub, opt.use_existing_scratch_files, tsep_k_pi);

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<tsep_k_pi.size();tsep_k_pi_idx++){
      if(scratch_store.doSkipLoad(tsep_k_pi_idx)) continue;
      getDataSampleAMA(A0_all_j, A0_all_dj, tsep_k_pi[tsep_k_pi_idx], bubble_data, 
		       data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		       data_dir_S, traj_start_S, traj_lessthan_S,
		       data_dir_C, traj_start_C, traj_lessthan_C,
		       traj_inc, bin_size,
		       Lt, tsep_pipi, resamplers, opt.read_opts);

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

void checkpointRawOnly(const std::vector<int> &tsep_k_pi,

		       const std::string &bubble_file_fmt_sloppy, const std::string &bubble_file_fmt_exact,
		       const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,
		       
		       const std::vector<std::string> &data_file_fmt_sloppy,  const std::vector<std::string> &data_file_fmt_exact,
		       const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,

		       const std::string &data_dir_S, const int traj_start_S, const int traj_lessthan_S,
		       const std::string &data_dir_C, const int traj_start_C, const int traj_lessthan_C,
		       const int traj_inc, const int bin_size, 
		       const int Lt, const int tsep_pipi, const readKtoPiPiDataSampleAMAoptions &opt){
  assert(opt.read_opts_sloppy_S.save_data_checkpoint &&
	 opt.read_opts_sloppy_C.save_data_checkpoint &&
	 opt.read_opts_exact_C.save_data_checkpoint);

  sampleAMA_resamplers resamplers(traj_start_S, traj_lessthan_S,
				  traj_start_C, traj_lessthan_C,
				  traj_inc, bin_size);
  
  allBubbleData bubble_data(bubble_file_fmt_sloppy, bubble_file_fmt_exact, bubble_pimom_proj,
			    data_dir_S, traj_start_S, traj_lessthan_S,
			    data_dir_C, traj_start_C, traj_lessthan_C,
			    traj_inc, bin_size, Lt, tsep_pipi, resamplers, opt);

  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<tsep_k_pi.size();tsep_k_pi_idx++){
    std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
    allRawData raw(bubble_data, tsep_k_pi[tsep_k_pi_idx], 
		   data_file_fmt_sloppy, data_file_fmt_exact, type1_pimom_proj,
		   data_dir_S, traj_start_S, traj_lessthan_S, 
		   data_dir_C, traj_start_C, traj_lessthan_C, 
		   traj_inc, bin_size,
		   Lt,tsep_pipi,opt);
  }
}

CPSFIT_END_NAMESPACE

#endif

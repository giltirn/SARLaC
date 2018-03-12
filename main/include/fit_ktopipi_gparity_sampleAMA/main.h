#ifndef FIT_KTOPIPI_GPARITY_SAMPLEAMA_MAIN_H__
#define FIT_KTOPIPI_GPARITY_SAMPLEAMA_MAIN_H__


//#define ALPHA_AND_VAC_SUB_SEPARATE
#define VAC_SUB_SEPARATE

#ifdef VAC_SUB_SEPARATE

template<typename DistributionType>
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const allRawData &raw, const allBubbleData &bubble_data, const allInputs &inputs, const std::string &descr){
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
						      raw, getResampledBubbleSampleAMA<DistributionType>::get(bubble_data),q,tsep_k_pi,inputs);


  //Get vaccum subtractions separately
  std::array<NumericTensor<DistributionType,1>, 3> A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r; //[ens-status][t]
  
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(dummy, A0_type4_srcavg_vacsub_r[sloppy_S], mix4_srcavg_vacsub_r[sloppy_S],
					     raw, bubble_data,q,tsep_k_pi, 'S',Sloppy, inputs);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(dummy, A0_type4_srcavg_vacsub_r[sloppy_C], mix4_srcavg_vacsub_r[sloppy_C],
					     raw, bubble_data,q,tsep_k_pi, 'C',Sloppy, inputs);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(dummy, A0_type4_srcavg_vacsub_r[exact_C], mix4_srcavg_vacsub_r[exact_C],
					     raw, bubble_data,q,tsep_k_pi, 'C',Exact, inputs);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  std::array< IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> , 3> mix_srcavg_r; //[ens-status][t]

  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,inputs); //[t]  (does full correction)
  for(int i=3;i<=4;i++){
    mix_srcavg_r[sloppy_S](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,'S',Sloppy,raw,inputs);
    mix_srcavg_r[sloppy_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,'C',Sloppy,raw,inputs);
    mix_srcavg_r[exact_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,'C',Exact,raw,inputs);
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
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const allRawData &raw, const allBubbleData &bubble_data, const allInputs &inputs, const std::string &descr){
  const int sloppy_S = 0;
  const int sloppy_C = 1;
  const int exact_C = 2;

  const int Lt = inputs.args.Lt;
  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  std::array<NumericTensor<DistributionType,1>, 3> alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r; //[ens-status][t]
  
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(alpha_r[sloppy_S], A0_type4_srcavg_vacsub_r[sloppy_S], mix4_srcavg_vacsub_r[sloppy_S],
					     raw, bubble_data,q,tsep_k_pi, 'S',Sloppy, inputs);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(alpha_r[sloppy_C], A0_type4_srcavg_vacsub_r[sloppy_C], mix4_srcavg_vacsub_r[sloppy_C],
					     raw, bubble_data,q,tsep_k_pi, 'C',Sloppy, inputs);
  computeAlphaAndVacuumSubtractionsSampleAMASeparately(alpha_r[exact_C], A0_type4_srcavg_vacsub_r[exact_C], mix4_srcavg_vacsub_r[exact_C],
					     raw, bubble_data,q,tsep_k_pi, 'C',Exact, inputs);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  std::array< IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> , 3> mix_srcavg_r; //[ens-status][t]

  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,inputs); //[t]  (does full correction)
  for(int i=3;i<=4;i++){
    mix_srcavg_r[sloppy_S](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,'S',Sloppy,raw,inputs);
    mix_srcavg_r[sloppy_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,'C',Sloppy,raw,inputs);
    mix_srcavg_r[exact_C](i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,'C',Exact,raw,inputs);
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
NumericTensor<DistributionType,1> computeQamplitudeSampleAMA(const int q, const int tsep_k_pi, const allRawData &raw, const allBubbleData &bubble_data, const allInputs &inputs, const std::string &descr){
  const int Lt = inputs.args.Lt;
  //Compute alpha and type4/mix4 vacuum subtractions
  std::cout << "Computing " << descr << " alpha and vacuum subtractions\n";
  NumericTensor<DistributionType,1> alpha_r({Lt}), A0_type4_srcavg_vacsub_r({Lt}), mix4_srcavg_vacsub_r({Lt}); //[t]
  computeAlphaAndVacuumSubtractionsSampleAMACorrected(alpha_r, A0_type4_srcavg_vacsub_r, mix4_srcavg_vacsub_r,
						      raw, getResampledBubbleSampleAMA<DistributionType>::get(bubble_data),q,tsep_k_pi,inputs);

  //Compute tK-averages type4 and mix4 diagrams from data including bubble-------------//
  std::cout << "Computing " << descr << " tK averages and mix diagrams\n";
  IndexedContainer<NumericTensor<DistributionType,1>, 4, 1> A0_srcavg_r; //[t]
  IndexedContainer<NumericTensor<DistributionType,1>, 2, 3> mix_srcavg_r; //[t]
  for(int i=1;i<=4;i++) A0_srcavg_r(i) = resampleAverageTypeDataSampleAMA<DistributionType>(i,q,raw,inputs); //[t]
  for(int i=3;i<=4;i++) mix_srcavg_r(i) = resampleAverageMixDiagramSampleAMA<DistributionType>(i,raw,inputs);

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
		      const allBubbleData &bubble_data, const int tsep_k_pi_idx, const allInputs &inputs){
  int tsep_k_pi = inputs.args.tsep_k_pi[tsep_k_pi_idx];
  std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
  printMem("getData called");
  
  allRawData raw(bubble_data, tsep_k_pi, inputs);
  
  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
 
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = computeQamplitudeSampleAMA<doubleJackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, inputs, "double jackknife");
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = computeQamplitudeSampleAMA<jackknifeDistributionD>(q, tsep_k_pi, raw, bubble_data, inputs, "single jackknife");

    //Insert data into output containers    
    for(int t=0;t<inputs.args.Lt;t++){
      A0_all_j[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_j(&t));
      A0_all_dj[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_dj(&t));
    }
  }
}



void getDataSampleAMA(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, 
		      std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		      const allInputs &inputs){
  if(inputs.cmdline.load_amplitude_data){
#ifdef HAVE_HDF5
    HDF5reader reader(inputs.cmdline.load_amplitude_data_file);
    read(reader, A0_all_j, "A0_all_j");
    read(reader, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    //Read the bubble data
    allBubbleData bubble_data(inputs);
    scratch scratch_store(inputs.cmdline.use_scratch, inputs.cmdline.use_scratch_stub, inputs.cmdline.use_existing_scratch_files, inputs.args.tsep_k_pi);

    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<inputs.args.tsep_k_pi.size();tsep_k_pi_idx++){
      if(scratch_store.doSkipLoad(tsep_k_pi_idx)) continue;
      getDataSampleAMA(A0_all_j,A0_all_dj,bubble_data,tsep_k_pi_idx,inputs);
      scratch_store.writeScratch(A0_all_j, A0_all_dj, tsep_k_pi_idx);
    }	

    scratch_store.reloadScratch(A0_all_j, A0_all_dj); 
  }

  if(inputs.cmdline.save_amplitude_data){
#ifdef HAVE_HDF5
    HDF5writer writer(inputs.cmdline.save_amplitude_data_file);
    write(writer, A0_all_j, "A0_all_j");
    write(writer, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }
}

void checkpointRawOnly(const allInputs &inputs){
  allBubbleData bubble_data(inputs);
  for(int tsep_k_pi_idx=0;tsep_k_pi_idx<inputs.args.tsep_k_pi.size();tsep_k_pi_idx++){
    int tsep_k_pi = inputs.args.tsep_k_pi[tsep_k_pi_idx];
    std::cout << "Getting data for tsep_k_pi = " <<  tsep_k_pi << std::endl;
    allRawData raw(bubble_data, tsep_k_pi, inputs);
  }
}
#endif

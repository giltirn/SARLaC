#ifndef _FIT_KTOPIPI_MAIN_H
#define _FIT_KTOPIPI_MAIN_H

#include<distribution_iterate.h>

template<typename T, int N>
struct iterate<NumericTensor<T,N> >{
  static inline int size(const NumericTensor<T,N> &from){ 
    int sz = 1;
    for(int i=0;i<N;i++) sz *= from.size(i);
    return sz;
  }
  static inline std::vector<int> unmap(int i, const NumericTensor<T,N> &from){ 
    std::vector<int> coord(N);
    for(int d=N-1;d>=0;d--){
      coord[d] = i % from.size(d);
      i /= from.size(d);
    }
    return coord;
  }    
  static inline const T& at(const int i, const NumericTensor<T,N> &from){
    std::vector<int> coord = unmap(i,from);
    return from(coord.data());
  }
  static inline T & at(const int i, NumericTensor<T,N> &from){
    std::vector<int> coord = unmap(i,from);
    return from(coord.data());
  }
};  

template<typename DistributionType, typename Accessor>
void average(DistributionType & into, const Accessor &data, const int size){
  //boost::timer::auto_cpu_timer total_time("average(DistributionType & into, const Accessor &data, const int size) %w s\n");
  assert(size > 0);
  into = data(0);
  for(int i=1;i<size;i++) into = into+data(i);
  into = into/double(size);
}

template<typename DistributionType, typename Accessor>
void resampleAverage(DistributionType & into, const Accessor &data, const int size){
  assert(size > 0);
  into.resample(data(0));
  DistributionType tmp;
  for(int i=1;i<size;i++){
    tmp.resample(data(i));
    into = into+tmp;
  }
  into = into/double(size);
}

template<typename resampledDistributionType>
void computeAlphaAndVacuumSubtractions(NumericTensor<resampledDistributionType,1> &alpha,
				       NumericTensor<resampledDistributionType,1> &A0_type4_srcavg_vacsub,
				       NumericTensor<resampledDistributionType,1> &mix4_srcavg_vacsub,
				       const NumericTensor<rawDataDistributionD,3> &A0_type4_nobub_alltK,
				       const NumericTensor<rawDataDistributionD,2> &mix4_nobub_alltK,
				       const NumericTensor<resampledDistributionType,1> &bubble_rs,
				       const int q,
				       const std::vector<int> &type4_nonzerotK,
				       const int tsep_k_pi,
				       const Args &args){
  //Compute mix4 double-jackknife and tK average
  resampledDistributionType zro = bubble_rs({0}); zeroit(zro);
  
#pragma omp parallel for
  for(int t=0;t<args.Lt;t++){
    mix4_srcavg_vacsub(&t) = zro;
    A0_type4_srcavg_vacsub(&t) = zro;

    resampledDistributionType mix4_nobub_srcavg = zro;
    resampledDistributionType A0_type4_nobub_srcavg = zro;
    
    for(int ii=0;ii<type4_nonzerotK.size();ii++){
      const int tK = type4_nonzerotK[ii];
      const int tB = (tK + tsep_k_pi + args.tsep_pipi) % args.Lt;
      
      resampledDistributionType mix4_nobub_rs;  mix4_nobub_rs.resample(mix4_nobub_alltK({tK,t}) );
      resampledDistributionType mix4_vacsub_rs = mix4_nobub_rs * bubble_rs(&tB);

      resampledDistributionType A0_type4_nobub_rs;  A0_type4_nobub_rs.resample(A0_type4_nobub_alltK({q,tK,t}) );
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


template<typename DistributionType>
NumericTensor<DistributionType,1> resampleAverageTypeData(const NumericTensor<rawDataDistributionD,3> &typedata_alltK,
							  const int q,
							  const std::vector<int> &typedata_nonzerotK,
							  const Args &args){
  NumericTensor<DistributionType,1> out({args.Lt}); //[t]
  for(int t=0;t<args.Lt;t++)
    resampleAverage(out(&t), [&](const int i){ return typedata_alltK({q,typedata_nonzerotK[i],t}); }, typedata_nonzerotK.size());
  return out;
}

template<typename DistributionType>
NumericTensor<DistributionType,1> resampleAverageMixDiagram(const NumericTensor<rawDataDistributionD,2> &mixdata_alltK,
							    const std::vector<int> &mixdata_nonzerotK,
							    const Args &args){
  NumericTensor<DistributionType,1> out({args.Lt}); //[t]
  for(int t=0;t<args.Lt;t++)
    resampleAverage(out(&t), [&](const int i){ return mixdata_alltK({mixdata_nonzerotK[i],t}); }, mixdata_nonzerotK.size());
  return out;
}


//Write results at various stages for debugging against other codes

template<typename DistributionType,int N>
void writeToTextFile(const NumericTensor<DistributionType,N> &m, const std::string &filename){
  std::ofstream of(filename.c_str());
  int NN = iterate<NumericTensor<DistributionType,N> >::size(m);
  for(int i=0;i<NN;i++){
    std::vector<int> coord = iterate<NumericTensor<DistributionType,N> >::unmap(i,m);
    const DistributionType &dist = m(coord.data());
    int SS = iterate<DistributionType>::size(dist);
    for(int s=0;s<SS;s++){
      std::vector<int> sample = iterate<DistributionType>::unmap(s,dist);
      for(int cc=0;cc<coord.size();cc++) of << coord[cc] << " ";
      for(int ss=0;ss<sample.size();ss++) of << sample[ss] << " ";
      of << iterate<DistributionType>::at(s,dist) << std::endl;
    }
  }
  of.close();
}



//Read and prepare the data for a particular tsep_k_pi_idx
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_all_dj,
	     const NumericTensor<rawDataDistributionD,1> &bubble, const NumericTensor<doubleJackknifeDistributionD,1> &bubble_dj, const NumericTensor<jackknifeDistributionD,1> &bubble_j,
	     const int tsep_k_pi_idx, const Args &args, const CMDline &cmdline){
  std::cout << "Getting data for tsep_k_pi = " <<  args.tsep_k_pi[tsep_k_pi_idx] << std::endl;
  printMem("getData called");
  int tsep_k_pi = args.tsep_k_pi[tsep_k_pi_idx];

  type1234Data type1, type2, type3, type4;

  if(cmdline.load_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << cmdline.load_data_checkpoint_stub << "_tsepkpi" << args.tsep_k_pi[tsep_k_pi_idx] << ".hdf5";
    std::cout << "Loading checkpoint data for tsep_k_pi = " << args.tsep_k_pi[tsep_k_pi_idx] << " from " << file.str() << std::endl;
    HDF5reader rd(file.str());
    read(rd,type1,"type1");
    read(rd,type2,"type2");
    read(rd,type3,"type3");
    read(rd,type4,"type4");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }else{
    type1 = readType(1, args.traj_start, args.traj_inc, args.traj_lessthan, tsep_k_pi, args.tsep_pipi, args.Lt, args.data_dir);
    type2 = readType(2, args.traj_start, args.traj_inc, args.traj_lessthan, tsep_k_pi, args.tsep_pipi, args.Lt, args.data_dir);
    type3 = readType(3, args.traj_start, args.traj_inc, args.traj_lessthan, tsep_k_pi, args.tsep_pipi, args.Lt, args.data_dir);
    type4 = readType(4, args.traj_start, args.traj_inc, args.traj_lessthan, tsep_k_pi, args.tsep_pipi, args.Lt, args.data_dir);
  }

  if(cmdline.save_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << cmdline.save_data_checkpoint_stub << "_tsepkpi" << args.tsep_k_pi[tsep_k_pi_idx] << ".hdf5";
    std::cout << "Saving checkpoint data for tsep_k_pi = " << tsep_k_pi << " to " << file.str() << std::endl;
    HDF5writer wr(file.str());
    write(wr,type1,"type1");
    write(wr,type2,"type2");
    write(wr,type3,"type3");
    write(wr,type4,"type4");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }
  
  const int nsample = type1.getNsample();
  
  std::vector<int> type1_nonzerotK = type1.getNonZeroKaonTimeslices();
  std::cout << "Type1 data non-zero for tK=" << type1_nonzerotK << std::endl;
  std::vector<int> type2_nonzerotK = type2.getNonZeroKaonTimeslices();
  std::cout << "Type2 data non-zero for tK=" << type2_nonzerotK << std::endl;
  std::vector<int> type3_nonzerotK = type3.getNonZeroKaonTimeslices();
  std::cout << "Type3 data non-zero for tK=" << type3_nonzerotK << std::endl;
  std::vector<int> type4_nonzerotK = type4.getNonZeroKaonTimeslices();
  std::cout << "Type4 data non-zero for tK=" << type4_nonzerotK << std::endl;

  //Get the type1-4 and mix3/mix4 components of the data.
  //Note that we have not yet multiplied the type4/mix4 data by the pipi bubble, hence these  are just the K->vac component which are needed for computing alpha
  std::cout << "Computing raw data diagrams prior to multiplying by bubble\n";
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK = computeAmplitudeType1<computeAmplitudeAlltKtensorControls>(type1); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK = computeAmplitudeType2<computeAmplitudeAlltKtensorControls>(type2); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type3_alltK = computeAmplitudeType3<computeAmplitudeAlltKtensorControls>(type3); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_nobub = computeAmplitudeType4<computeAmplitudeAlltKtensorControls>(type4); //[Qidx][tK][t]
  
  NumericTensor<rawDataDistributionD,2> mix3_alltK({args.Lt,args.Lt}); //[tK][t]
  NumericTensor<rawDataDistributionD,2> mix4_alltK_nobub({args.Lt,args.Lt}); //[tK][t]
  for(int tK=0;tK<args.Lt;tK++)
    for(int t=0;t<args.Lt;t++){
      mix3_alltK({tK,t}) = type3(tK,t).mix();
      mix4_alltK_nobub({tK,t}) = type4(tK,t).mix();
    }

  printMem("Pre typedata free");

  //Data is no longer needed, so free it
  type1.freeData(); type2.freeData(); type3.freeData(); type4.freeData();

  printMem("Post typedata free");
  
  //Compute the type4/mix4 data with the bubble included
  std::cout << "Computing raw type4/mix4 data with bubble included" << std::endl;
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK({10,args.Lt,args.Lt});
  NumericTensor<rawDataDistributionD,2> mix4_alltK({args.Lt,args.Lt});
  for(int tK=0;tK<args.Lt;tK++)
    for(int t=0;t<args.Lt;t++){
      int tB = (tK + tsep_k_pi + args.tsep_pipi) % args.Lt;
      mix4_alltK({tK,t}) = mix4_alltK_nobub({tK,t})*bubble(&tB);
      for(int q=0;q<10;q++)
	A0_type4_alltK({q,tK,t}) = A0_type4_alltK_nobub({q,tK,t})*bubble(&tB);
    }
  
  for(int q=0;q<10;q++){
    std::cout << "Starting Q" << q+1 << std::endl;

    printMem("Starting new Q");
    
    //Compute alpha and type4/mix4 vacuum subtractions
    std::cout << "Computing double-jackknife alpha and vacuum subtractions\n";
    NumericTensor<doubleJackknifeDistributionD,1> alpha_dj({args.Lt}); //[t]
    NumericTensor<doubleJackknifeDistributionD,1> A0_type4_srcavg_vacsub_dj({args.Lt}); //[t]
    NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_vacsub_dj({args.Lt}); //[t]
    computeAlphaAndVacuumSubtractions(alpha_dj,A0_type4_srcavg_vacsub_dj,mix4_srcavg_vacsub_dj,
				      A0_type4_alltK_nobub,mix4_alltK_nobub,bubble_dj,q,type4_nonzerotK,tsep_k_pi,args);

    std::cout << "Computing single-jackknife alpha and vacuum subtractions\n";
    NumericTensor<jackknifeDistributionD,1> alpha_j({args.Lt}); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type4_srcavg_vacsub_j({args.Lt}); //[t]
    NumericTensor<jackknifeDistributionD,1> mix4_srcavg_vacsub_j({args.Lt}); //[t]
    computeAlphaAndVacuumSubtractions(alpha_j,A0_type4_srcavg_vacsub_j,mix4_srcavg_vacsub_j,
				      A0_type4_alltK_nobub,mix4_alltK_nobub,bubble_j,q,type4_nonzerotK,tsep_k_pi,args);

    //Compute double-jackknife, tK-averages type4 and mix4 diagrams from data including bubble-------------//
    std::cout << "Computing double-jackknife tK averages and mix diagrams\n";
    NumericTensor<doubleJackknifeDistributionD,1> A0_type1_srcavg_dj = resampleAverageTypeData<doubleJackknifeDistributionD>(A0_type1_alltK, q, type1_nonzerotK, args); //[t]
    NumericTensor<doubleJackknifeDistributionD,1> A0_type2_srcavg_dj = resampleAverageTypeData<doubleJackknifeDistributionD>(A0_type2_alltK, q, type2_nonzerotK, args); //[t]
    NumericTensor<doubleJackknifeDistributionD,1> A0_type3_srcavg_dj = resampleAverageTypeData<doubleJackknifeDistributionD>(A0_type3_alltK, q, type3_nonzerotK, args); //[t]
    NumericTensor<doubleJackknifeDistributionD,1> A0_type4_srcavg_dj = resampleAverageTypeData<doubleJackknifeDistributionD>(A0_type4_alltK, q, type4_nonzerotK, args); //[t]
    NumericTensor<doubleJackknifeDistributionD,1> mix3_srcavg_dj = resampleAverageMixDiagram<doubleJackknifeDistributionD>(mix3_alltK, type3_nonzerotK, args);
    NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_dj = resampleAverageMixDiagram<doubleJackknifeDistributionD>(mix4_alltK, type4_nonzerotK, args);

    std::cout << "Computing single-jackknife tK averages and mix diagrams\n";
    NumericTensor<jackknifeDistributionD,1> A0_type1_srcavg_j = resampleAverageTypeData<jackknifeDistributionD>(A0_type1_alltK, q, type1_nonzerotK, args); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type2_srcavg_j = resampleAverageTypeData<jackknifeDistributionD>(A0_type2_alltK, q, type2_nonzerotK, args); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type3_srcavg_j = resampleAverageTypeData<jackknifeDistributionD>(A0_type3_alltK, q, type3_nonzerotK, args); //[t]
    NumericTensor<jackknifeDistributionD,1> A0_type4_srcavg_j = resampleAverageTypeData<jackknifeDistributionD>(A0_type4_alltK, q, type4_nonzerotK, args); //[t]
    NumericTensor<jackknifeDistributionD,1> mix3_srcavg_j = resampleAverageMixDiagram<jackknifeDistributionD>(mix3_alltK, type3_nonzerotK, args);
    NumericTensor<jackknifeDistributionD,1> mix4_srcavg_j = resampleAverageMixDiagram<jackknifeDistributionD>(mix4_alltK, type4_nonzerotK, args);

    //Subtract the pseudoscalar operators and mix4 vacuum term
    std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under double-jackknife\n";
    A0_type3_srcavg_dj = A0_type3_srcavg_dj.transform([&](int const* t, const doubleJackknifeDistributionD &from){ return doubleJackknifeDistributionD(from - alpha_dj(t)*mix3_srcavg_dj(t)); }); 
    A0_type4_srcavg_dj = A0_type4_srcavg_dj.transform(
						      [&](int const* t, const doubleJackknifeDistributionD &from){
							return doubleJackknifeDistributionD(from - alpha_dj(t)*( mix4_srcavg_dj(t) - mix4_srcavg_vacsub_dj(t) ) );
						      }
						      ); 
    std::cout << "Subtracting pseudoscalar operators and mix4 vacuum term under single-jackknife\n";
    A0_type3_srcavg_j = A0_type3_srcavg_j.transform([&](int const* t, const jackknifeDistributionD &from){ return jackknifeDistributionD(from - alpha_j(t)*mix3_srcavg_j(t)); }); 
    A0_type4_srcavg_j = A0_type4_srcavg_j.transform(
						    [&](int const* t, const jackknifeDistributionD &from){
						      return jackknifeDistributionD(from - alpha_j(t)*( mix4_srcavg_j(t) - mix4_srcavg_vacsub_j(t) ) );
						    }
						    ); 
    //Perform the type 4 vacuum subtraction
    std::cout << "Performing type-4 vacuum subtraction\n";
    A0_type4_srcavg_dj = A0_type4_srcavg_dj - A0_type4_srcavg_vacsub_dj;
    A0_type4_srcavg_j = A0_type4_srcavg_j - A0_type4_srcavg_vacsub_j;
  
    //Get the full double-jackknife amplitude
    std::cout << "Computing full amplitudes\n";
    NumericTensor<doubleJackknifeDistributionD,1> A0_full_srcavg_dj = A0_type1_srcavg_dj + A0_type2_srcavg_dj + A0_type3_srcavg_dj + A0_type4_srcavg_dj;
    NumericTensor<jackknifeDistributionD,1> A0_full_srcavg_j = A0_type1_srcavg_j + A0_type2_srcavg_j + A0_type3_srcavg_j + A0_type4_srcavg_j;

    //Insert data into output containers    
    for(int t=0;t<args.Lt;t++){
      A0_all_j[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_j(&t));
      A0_all_dj[q].push_back(amplitudeDataCoord(t,tsep_k_pi), A0_full_srcavg_dj(&t));
    }
  }
}

void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j, std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > &A0_all_dj,
	     const Args &args, const CMDline &cmdline){
  if(cmdline.load_amplitude_data){
#ifdef HAVE_HDF5
    HDF5reader reader(cmdline.load_amplitude_data_file);
    read(reader, A0_all_j, "A0_all_j");
    read(reader, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Reading amplitude data requires HDF5\n");
#endif
  }else{
    //Read the bubble data
    NumericTensor<rawDataDistributionD,1> bubble = getA2projectedBubble(args,cmdline);
    NumericTensor<jackknifeDistributionD,1> bubble_j = bubble.transform(resampleFunctor<jackknifeDistributionD,rawDataDistributionD>());
    NumericTensor<doubleJackknifeDistributionD,1> bubble_dj = bubble.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
    
    //Read and prepare the amplitude data for fitting

    //Setup scratch space if in use
    std::vector<std::string> scratch_files(args.tsep_k_pi.size());
    if(cmdline.use_scratch){
#ifndef HAVE_HDF5
      error_exit("getData: scratch usage requires HDF5\n");
#endif
      for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
	std::ostringstream f; f << "scratch_" << tsep_k_pi_idx;
	scratch_files[tsep_k_pi_idx] = f.str();
      }
    }

    //Load data
    for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
      if(cmdline.use_scratch && cmdline.use_existing_scratch_files && fileExists(scratch_files[tsep_k_pi_idx])) continue;
      
      getData(A0_all_j,A0_all_dj,bubble,bubble_dj,bubble_j,tsep_k_pi_idx,args,cmdline);
      if(cmdline.use_scratch){
	printMem("Pre-scratch write");
#ifdef HAVE_HDF5
	HDF5writer writer(scratch_files[tsep_k_pi_idx]);
	write(writer, A0_all_j, "A0_all_j");
	write(writer, A0_all_dj, "A0_all_dj");
	//Clear the data
	for(int q=0;q<10;q++){
	  A0_all_j[q].clear(); A0_all_dj[q].clear();
	}
#endif
	printMem("Post-scratch write");
      }
    }	

    //Reload scratch data
    if(cmdline.use_scratch){
      for(int tsep_k_pi_idx=0;tsep_k_pi_idx<args.tsep_k_pi.size();tsep_k_pi_idx++){
#ifdef HAVE_HDF5
	HDF5reader reader(scratch_files[tsep_k_pi_idx]);
	std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > tmp_A0_all_j;
	std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> > tmp_A0_all_dj;
	read(reader, tmp_A0_all_j, "A0_all_j");
	read(reader, tmp_A0_all_dj, "A0_all_dj");
	for(int q=0;q<10;q++){
	  for(int i=0;i<tmp_A0_all_j[q].size();i++){
	    A0_all_j[q].push_back(tmp_A0_all_j[q][i]);
	    A0_all_dj[q].push_back(tmp_A0_all_dj[q][i]);
	  }
	}
#endif
      }
    }	
      
  }

  if(cmdline.save_amplitude_data){
#ifdef HAVE_HDF5
    HDF5writer writer(cmdline.save_amplitude_data_file);
    write(writer, A0_all_j, "A0_all_j");
    write(writer, A0_all_dj, "A0_all_dj");
#else
    error_exit("getData: Saving amplitude data requires HDF5\n");
#endif
  }

  
}



#endif

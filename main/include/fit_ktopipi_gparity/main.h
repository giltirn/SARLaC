#ifndef _FIT_KTOPIPI_MAIN_H
#define _FIT_KTOPIPI_MAIN_H

template<typename distributionType>
struct iterate;

template<typename T>
struct iterate<doubleJackknifeDistribution<T> >{
  static inline int size(const doubleJackknifeDistribution<T> &from){ return from.size() * (from.size()-1); } //j + (from.size()-1)*i
  static inline const T& at(const int i, const doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline T & at(const int i, doubleJackknifeDistribution<T> &from){
    const int nn = from.size()-1;
    return from.sample(i/nn).sample(i%nn);
  }
  static inline std::vector<int> unmap(const int i, const doubleJackknifeDistribution<T> &from){ 
    const int nn = from.size()-1;
    return std::vector<int>({i/nn, i%nn});
  }
};
template<typename T>
struct iterate<rawDataDistribution<T> >{
  static inline int size(const rawDataDistribution<T> &from){ return from.size(); } 
  static inline const T& at(const int i, const rawDataDistribution<T> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, rawDataDistribution<T> &from){
    return from.sample(i);
  }
  static inline std::vector<int> unmap(const int i, const rawDataDistribution<T> &from){ 
    return std::vector<int>({i});
  }
};
template<typename T>
struct iterate<jackknifeDistribution<T> >{
  static inline int size(const jackknifeDistribution<T> &from){ return from.size(); } 
  static inline const T& at(const int i, const jackknifeDistribution<T> &from){
    return from.sample(i);
  }
  static inline T & at(const int i, jackknifeDistribution<T> &from){
    return from.sample(i);
  }
  static inline std::vector<int> unmap(const int i, const jackknifeDistribution<T> &from){ 
    return std::vector<int>({i});
  }
};
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

template<typename distributionType>
struct multiplyBubbleFunctor{
  const NumericTensor<distributionType,1> &bubble;
  int tsep_k_pi;
  int tsep_pipi;
  int Lt;
  int tK_dim;
  multiplyBubbleFunctor(const NumericTensor<distributionType,1> &_bubble, const int _tsep_k_pi, const int _tsep_pipi, const int _Lt, const int _tK_dim): bubble(_bubble), tsep_k_pi(_tsep_k_pi), tsep_pipi(_tsep_pipi), Lt(_Lt), tK_dim(_tK_dim){}
  inline distributionType operator()(int const* coord, const distributionType &from) const{
    typedef iterate<distributionType> iter;
    int tK = coord[tK_dim]; int tB = (tK + tsep_k_pi + tsep_pipi) % Lt;
    distributionType out(from);
    for(int i=0;i<iter::size(out);i++)
      iter::at(i,out) = iter::at(i,out) * iter::at(i,bubble({tB}));
    return out;
  }
};


template<typename DistributionType, typename Accessor>
void average(DistributionType & into, const Accessor &data, const int size){
  assert(size > 0);
  into = data(0);
  for(int i=1;i<size;i++) into = into+data(i);
  into = into/double(size);
}

void computeAlphaAndVacuumSubtractions(NumericTensor<doubleJackknifeDistributionD,2> &alpha,
				       NumericTensor<doubleJackknifeDistributionD,2> &A0_type4_srcavg_vacsub,
				       NumericTensor<doubleJackknifeDistributionD,1> &mix4_srcavg_vacsub,
				       const NumericTensor<rawDataDistributionD,3> &A0_type4_nobub_alltK,
				       const NumericTensor<rawDataDistributionD,2> &mix4_nobub_alltK,
				       const NumericTensor<doubleJackknifeDistributionD,1> &bubble_dj,
				       const std::vector<int> &type4_nonzerotK,
				       const int tsep_k_pi,
				       const Args &args){
  for(int t=0;t<args.Lt;t++){
    //Compute mix4 double-jackknife and tK average
    NumericTensor<doubleJackknifeDistributionD,1> mix4_nobub_alltK_dj({args.Lt}, [&](const int* coord){ doubleJackknifeDistributionD out; out.resample(mix4_nobub_alltK({coord[0],t})); return out; } ); //[tK]
    doubleJackknifeDistributionD mix4_nobub_srcavg_dj; average(mix4_nobub_srcavg_dj, [&](const int i){ return mix4_nobub_alltK_dj(&type4_nonzerotK[i]); },  type4_nonzerotK.size());

    //Compute mix4 vacuum subtraction
    NumericTensor<doubleJackknifeDistributionD,1> mix4_alltK_vacsub = mix4_nobub_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,tsep_k_pi,args.tsep_pipi,args.Lt,0)); //[tK]
    average(mix4_srcavg_vacsub({t}), [&](const int i){ return mix4_alltK_vacsub(&type4_nonzerotK[i]); },  type4_nonzerotK.size());
    
    for(int q=0;q<10;q++){
      //Compute type4 double-jackknife and tK average
      NumericTensor<doubleJackknifeDistributionD,1> A0_type4_nobub_alltK_dj({args.Lt}, [&](const int* coord){ doubleJackknifeDistributionD out; out.resample(A0_type4_nobub_alltK({q,coord[0],t})); return out; } ); //[tK]
      doubleJackknifeDistributionD A0_type4_nobub_srcavg_dj; average(A0_type4_nobub_srcavg_dj, [&](const int i){ return A0_type4_nobub_alltK_dj(&type4_nonzerotK[i]); },  type4_nonzerotK.size());

      //Compute alpha
      alpha({q,t}) = A0_type4_nobub_srcavg_dj/mix4_nobub_srcavg_dj;

      //Compute type4 vacuum subtraction
      NumericTensor<doubleJackknifeDistributionD,1> A0_type4_alltK_vacsub = A0_type4_nobub_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,tsep_k_pi,args.tsep_pipi,args.Lt,0));//[tK]
      average(A0_type4_srcavg_vacsub({q,t}), [&](const int i){ return A0_type4_alltK_vacsub(&type4_nonzerotK[i]); },  type4_nonzerotK.size());
    }
  }
}

NumericTensor<doubleJackknifeDistributionD,2> resampleAverageTypeData(const NumericTensor<rawDataDistributionD,3> &typedata_alltK,
								      const std::vector<int> &typedata_nonzerotK,
								      const Args &args){
  NumericTensor<doubleJackknifeDistributionD,2> out({10,args.Lt}); //[Qidx][t]
  for(int t=0;t<args.Lt;t++){     
    for(int q=0;q<10;q++){
      NumericTensor<doubleJackknifeDistributionD,1> typedata_alltK_dj({args.Lt}, [&](const int* coord){ doubleJackknifeDistributionD out; out.resample(typedata_alltK({q,coord[0],t})); return out; } ); //[tK]
      average(out({q,t}), [&](const int i){ return typedata_alltK_dj(&typedata_nonzerotK[i]); },  typedata_nonzerotK.size());
    }
  }
  return out;
}

NumericTensor<doubleJackknifeDistributionD,1> resampleAverageMixDiagram(const NumericTensor<rawDataDistributionD,2> &mixdata_alltK,
									const std::vector<int> &mixdata_nonzerotK,
									const Args &args){
  NumericTensor<doubleJackknifeDistributionD,1> out({args.Lt}); //[t]
  for(int t=0;t<args.Lt;t++){    
    NumericTensor<doubleJackknifeDistributionD,1> mixdata_alltK_dj({args.Lt}, [&](const int* coord){ doubleJackknifeDistributionD out; out.resample(mixdata_alltK({coord[0],t})); return out; } ); //[tK]
    average(out({t}), [&](const int i){ return mixdata_alltK_dj(&mixdata_nonzerotK[i]); },  mixdata_nonzerotK.size());    
  }
  return out;
}

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
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit, std::vector<NumericVector<jackknifeDistributionD> > &sigma_fit,
	     const NumericTensor<rawDataDistributionD,1> &bubble, const NumericTensor<doubleJackknifeDistributionD,1> &bubble_dj,
	     const int tsep_k_pi_idx, const Args &args, const CMDline &cmdline){
  std::cout << "Getting data for tsep_k_pi = " <<  args.tsep_k_pi[tsep_k_pi_idx] << std::endl;
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
  
  //Compute alpha and type4/mix4 vacuum subtractions
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK = computeAmplitudeType4<computeAmplitudeAlltKtensorControls>(type4); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,2> mix4_alltK({args.Lt,args.Lt}); //[tK][t]
  for(int tK=0;tK<args.Lt;tK++) for(int t=0;t<args.Lt;t++) mix4_alltK({tK,t}) = type4(tK,t).mix();
  
  NumericTensor<doubleJackknifeDistributionD,2> alpha({10,args.Lt}); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_vacsub({10,args.Lt}); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_vacsub({args.Lt}); //[t]
  computeAlphaAndVacuumSubtractions(alpha,A0_type4_srcavg_vacsub,mix4_srcavg_vacsub,A0_type4_alltK,mix4_alltK,bubble_dj,type4_nonzerotK,tsep_k_pi,args);
        
  //Save some quantities to disk for external viewing pleasure
  {
    std::ostringstream os; os << "alpha_tsep_k_pi" << tsep_k_pi << ".dat";
    writeToTextFile(alpha,os.str());
  }
  {
    std::ostringstream os; os << "type4_vacsub_tsep_k_pi" << tsep_k_pi << ".dat";
    writeToTextFile(A0_type4_srcavg_vacsub,os.str());
  }
  {
    std::ostringstream os; os << "mix4_vacsub_tsep_k_pi" << tsep_k_pi << ".dat";
    writeToTextFile(mix4_srcavg_vacsub,os.str());
  }

  //Apply the bubble diagram to the type4 and mix4 data
  A0_type4_alltK = A0_type4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,tsep_k_pi,args.tsep_pipi,args.Lt,1));
  mix4_alltK = mix4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,tsep_k_pi,args.tsep_pipi,args.Lt,0));

  //Get type1, type2, type3 and mix3 data
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK = computeAmplitudeType1<computeAmplitudeAlltKtensorControls>(type1); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK = computeAmplitudeType2<computeAmplitudeAlltKtensorControls>(type2); //[Qidx][tK][t]
  NumericTensor<rawDataDistributionD,3> A0_type3_alltK = computeAmplitudeType3<computeAmplitudeAlltKtensorControls>(type3); //[Qidx][tK][t]

  NumericTensor<rawDataDistributionD,2> mix3_alltK({args.Lt,args.Lt}); //[tK][t]
  for(int tK=0;tK<args.Lt;tK++) for(int t=0;t<args.Lt;t++) mix3_alltK({tK,t}) = type3(tK,t).mix();

  //Compute double-jackknife, tK-averages type4 and mix4 diagrams from data including bubble-------------//
  NumericTensor<doubleJackknifeDistributionD,2> A0_type1_srcavg_dj = resampleAverageTypeData(A0_type1_alltK, type1_nonzerotK, args); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,2> A0_type2_srcavg_dj = resampleAverageTypeData(A0_type2_alltK, type2_nonzerotK, args);
  NumericTensor<doubleJackknifeDistributionD,2> A0_type3_srcavg_dj = resampleAverageTypeData(A0_type3_alltK, type3_nonzerotK, args);
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_dj = resampleAverageTypeData(A0_type4_alltK, type4_nonzerotK, args);
  NumericTensor<doubleJackknifeDistributionD,1> mix3_srcavg_dj = resampleAverageMixDiagram(mix3_alltK, type3_nonzerotK, args);
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_dj = resampleAverageMixDiagram(mix4_alltK, type4_nonzerotK, args);

  //Save to ASCII
  {
    NumericTensor<doubleJackknifeDistributionD,2> const* type_data[4] = {&A0_type1_srcavg_dj, &A0_type2_srcavg_dj, &A0_type3_srcavg_dj, &A0_type4_srcavg_dj};
    for(int type=0;type<4;type++){
      std::ostringstream os; os << "A0_type" << type+1 << "_dj_tsep_k_pi" << tsep_k_pi << ".dat";
      writeToTextFile(*type_data[type], os.str());
    }
  }

  //Subtract the pseudoscalar operators and mix4 vacuum term
  A0_type3_srcavg_dj = A0_type3_srcavg_dj.transform([&](int const* coord, const doubleJackknifeDistributionD &from){ return doubleJackknifeDistributionD(from - alpha(coord)*mix3_srcavg_dj(coord+1)); }); 
  A0_type4_srcavg_dj = A0_type4_srcavg_dj.transform(
						    [&](int const* coord, const doubleJackknifeDistributionD &from){
						      return doubleJackknifeDistributionD(from - alpha(coord)*( mix4_srcavg_dj(coord+1) - mix4_srcavg_vacsub(coord+1) ) );
						    }
						    ); 
  
  //Perform the type 4 vacuum subtraction
  A0_type4_srcavg_dj = A0_type4_srcavg_dj - A0_type4_srcavg_vacsub;

  //Get the full double-jackknife amplitude
  NumericTensor<doubleJackknifeDistributionD,2> A0_full_srcavg_dj = A0_type1_srcavg_dj + A0_type2_srcavg_dj + A0_type3_srcavg_dj + A0_type4_srcavg_dj;

  {
    std::ostringstream os; os << "A0_dj_tsep_k_pi" << tsep_k_pi << ".dat";
    writeToTextFile(A0_full_srcavg_dj,os.str());
  }

  
  //Get the single-elimination jackknife distributions
  NumericTensor<jackknifeDistributionD,2> A0_type1_srcavg_j = A0_type1_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type2_srcavg_j = A0_type2_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type3_srcavg_j = A0_type3_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type4_srcavg_j = A0_type4_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });

  NumericTensor<jackknifeDistributionD,2> A0_full_srcavg_j = A0_type1_srcavg_j + A0_type2_srcavg_j + A0_type3_srcavg_j + A0_type4_srcavg_j;

  NumericTensor<jackknifeDistributionD,2> sigma_j = A0_full_srcavg_dj.transform([](int const *c, const doubleJackknifeDistributionD &d){ return jackknifeDistributionD(sqrt(doubleJackknifeDistributionD::covariance(d,d))); });

  
  {
    std::ostringstream os; os << "sigma_dj_tsep_k_pi" << tsep_k_pi << ".dat";
    writeToTextFile(sigma_j,os.str());
  }
  {
    std::ostringstream os; os << "A0_jack_tsep_k_pi" << tsep_k_pi << ".dat";
    writeToTextFile(A0_full_srcavg_j,os.str());
  }

  //Separate out the data in the desired fit range
  for(int q=0;q<10;q++)
    for(int t=0;t<args.Lt;t++)
      if(t >= args.tmin_k_op && t <= tsep_k_pi-args.tmin_op_pi){
	amplitudeDataCoord coord(t,tsep_k_pi);
	A0_fit[q].push_back(coord, A0_full_srcavg_j({q,t}));
	sigma_fit[q].push_back(sigma_j({q,t}));
      }
}

#endif

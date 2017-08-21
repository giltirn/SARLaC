#ifndef _FIT_KTOPIPI_MAIN_H
#define _FIT_KTOPIPI_MAIN_H


//Functor for tensor reduction by average
template<typename T, int Rank>
struct averageDimensionFunctor{
  const int dim;
  std::vector<int> const* use;
  
  averageDimensionFunctor(const int _dim, std::vector<int> const*_use = NULL): dim(_dim), use(_use){}
  
  void operator()(T &o, int const *coord, const NumericTensor<T,Rank> &from) const{
    int full_coord[Rank];
    int i=0; for(int ii=0;ii<Rank;ii++) if(ii!=dim) full_coord[ii] = coord[i++];    
    zeroit(o);
    if(use != NULL){
      assert(use->size()> 0);
      full_coord[dim] = use->at(0);
      o = from(full_coord);      
      for(int i=1;i<use->size();i++){
	full_coord[dim] = use->at(i);
	o = o + from(full_coord);
      }
      o = o/double(use->size());
    }else{
      assert(from.size(dim)>0);
      full_coord[dim] = 0;
      o = from(full_coord);
      for(int i=1;i<from.size(dim);i++){
	full_coord[dim] = i;
	o = o + from(full_coord);
      }
      o = o/double(from.size(dim));
    }
  }
};
    
template<typename Resampled, typename Raw>
struct resampleFunctor{
  inline Resampled operator()(int const* coord,const Raw &from) const{
    Resampled o(from.size());
    o.resample(from);
    return o;
  }
};

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



//Read and prepare the data for a particular tsep_k_pi_idx
void getData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit, std::vector<NumericVector<jackknifeDistributionD> > &sigma_fit,
	     const NumericTensor<rawDataDistributionD,1> &bubble, const NumericTensor<doubleJackknifeDistributionD,1> &bubble_dj,
	     const int tsep_k_pi_idx, const Args &args, const CMDline &cmdline){
  std::cout << "Getting data for tsep_k_pi = " <<  args.tsep_k_pi[tsep_k_pi_idx] << std::endl;

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
    type1 = readType(1, args.traj_start, args.traj_inc, args.traj_lessthan, args.tsep_k_pi[tsep_k_pi_idx], args.tsep_pipi, args.Lt, args.data_dir);
    type2 = readType(2, args.traj_start, args.traj_inc, args.traj_lessthan, args.tsep_k_pi[tsep_k_pi_idx], args.tsep_pipi, args.Lt, args.data_dir);
    type3 = readType(3, args.traj_start, args.traj_inc, args.traj_lessthan, args.tsep_k_pi[tsep_k_pi_idx], args.tsep_pipi, args.Lt, args.data_dir);
    type4 = readType(4, args.traj_start, args.traj_inc, args.traj_lessthan, args.tsep_k_pi[tsep_k_pi_idx], args.tsep_pipi, args.Lt, args.data_dir);
  }

  if(cmdline.save_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << cmdline.save_data_checkpoint_stub << "_tsepkpi" << args.tsep_k_pi[tsep_k_pi_idx] << ".hdf5";
    std::cout << "Saving checkpoint data for tsep_k_pi = " << args.tsep_k_pi[tsep_k_pi_idx] << " to " << file.str() << std::endl;
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
  std::vector<int> type2_nonzerotK = type2.getNonZeroKaonTimeslices();
  std::vector<int> type3_nonzerotK = type3.getNonZeroKaonTimeslices();
  std::vector<int> type4_nonzerotK = type4.getNonZeroKaonTimeslices();

  //Get the type 4 contribution without the bubble and double-jackknife resample
  NumericTensor<rawDataDistributionD,3> A0_type4_alltK = computeAmplitudeType4<computeAmplitudeAlltKtensorControls>(type4); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type4_alltK_dj = A0_type4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Get and double-jackknife resample the mix4 (no bubble) diagram
  NumericTensor<rawDataDistributionD,2> mix4_alltK({args.Lt,args.Lt}); //[tK][t]
  for(int tK=0;tK<args.Lt;tK++) for(int t=0;t<args.Lt;t++) mix4_alltK({tK,t}) = type4(tK,t).mix();
  NumericTensor<doubleJackknifeDistributionD,2> mix4_alltK_dj = mix4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  
  //tK-average the type4 and mix4 (no bubble) double-jackknife data then double-jackknife resample  
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_dj = A0_type4_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_dj = mix4_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]
   
  //Compute alpha from the above
  NumericTensor<doubleJackknifeDistributionD,2> alpha({10,args.Lt}, [&](int const *c){ return A0_type4_srcavg_dj({c[0],c[1]})/mix4_srcavg_dj({c[1]}); }); //[Qidx][t]

  //Compute vacuum subtractions
  NumericTensor<doubleJackknifeDistributionD,3> A0_type4_alltK_vacsub = A0_type4_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,args.tsep_k_pi[tsep_k_pi_idx],args.tsep_pipi,args.Lt,1));//[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,2> mix4_alltK_vacsub = mix4_alltK_dj.transform(multiplyBubbleFunctor<doubleJackknifeDistributionD>(bubble_dj,args.tsep_k_pi[tsep_k_pi_idx],args.tsep_pipi,args.Lt,0)); //[tK][t]

  //Source average the vacuum subtractions
  NumericTensor<doubleJackknifeDistributionD,2> A0_type4_srcavg_vacsub = A0_type4_alltK_vacsub.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  NumericTensor<doubleJackknifeDistributionD,1> mix4_srcavg_vacsub = mix4_alltK_vacsub.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]
  
  //Apply the bubble diagram to the type4 and mix4 data
  A0_type4_alltK = A0_type4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,args.tsep_k_pi[tsep_k_pi_idx],args.tsep_pipi,args.Lt,1));
  mix4_alltK = mix4_alltK.transform(multiplyBubbleFunctor<rawDataDistributionD>(bubble,args.tsep_k_pi[tsep_k_pi_idx],args.tsep_pipi,args.Lt,0));

  //Recompute double-jackknife resamplings of full type4 and mix4
  A0_type4_alltK_dj = A0_type4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  mix4_alltK_dj = mix4_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());

  //Recompute double-jackknife resamplings of full, source-averaged type4 and mix4
  A0_type4_srcavg_dj = A0_type4_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type4_nonzerotK)); //[Qidx][t]
  mix4_srcavg_dj = mix4_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type4_nonzerotK)); //[t]

  //Get, double-jackknife resample and source-average the type1, type2, type3 and mix3 contributions
  NumericTensor<rawDataDistributionD,3> A0_type1_alltK = computeAmplitudeType1<computeAmplitudeAlltKtensorControls>(type1); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type1_alltK_dj = A0_type1_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type1_srcavg_dj = A0_type1_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type1_nonzerotK)); //[Qidx][t]
  
  NumericTensor<rawDataDistributionD,3> A0_type2_alltK = computeAmplitudeType2<computeAmplitudeAlltKtensorControls>(type2); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type2_alltK_dj = A0_type2_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type2_srcavg_dj = A0_type2_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type2_nonzerotK)); //[Qidx][t]

  NumericTensor<rawDataDistributionD,3> A0_type3_alltK = computeAmplitudeType3<computeAmplitudeAlltKtensorControls>(type3); //[Qidx][tK][t]
  NumericTensor<doubleJackknifeDistributionD,3> A0_type3_alltK_dj = A0_type3_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,2> A0_type3_srcavg_dj = A0_type3_alltK_dj.reduce(1, averageDimensionFunctor<doubleJackknifeDistributionD,3>(1, &type3_nonzerotK)); //[Qidx][t]

  NumericTensor<rawDataDistributionD,2> mix3_alltK({args.Lt,args.Lt}); //[tK][t]
  for(int tK=0;tK<args.Lt;tK++) for(int t=0;t<args.Lt;t++) mix3_alltK({tK,t}) = type3(tK,t).mix();
  NumericTensor<doubleJackknifeDistributionD,2> mix3_alltK_dj = mix3_alltK.transform(resampleFunctor<doubleJackknifeDistributionD,rawDataDistributionD>());
  NumericTensor<doubleJackknifeDistributionD,1> mix3_srcavg_dj = mix3_alltK_dj.reduce(0, averageDimensionFunctor<doubleJackknifeDistributionD,2>(0, &type3_nonzerotK)); //[t]
  
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
  
  //Get the single-elimination jackknife distributions
  NumericTensor<jackknifeDistributionD,2> A0_type1_srcavg_j = A0_type1_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type2_srcavg_j = A0_type2_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type3_srcavg_j = A0_type3_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });
  NumericTensor<jackknifeDistributionD,2> A0_type4_srcavg_j = A0_type4_srcavg_dj.transform([](int const* coord, const doubleJackknifeDistributionD &from){ return from.toJackknife(); });

  NumericTensor<jackknifeDistributionD,2> A0_full_srcavg_j = A0_type1_srcavg_j + A0_type2_srcavg_j + A0_type3_srcavg_j + A0_type4_srcavg_j;

  NumericTensor<jackknifeDistributionD,2> sigma_j = A0_full_srcavg_dj.transform([](int const *c, const doubleJackknifeDistributionD &d){ return jackknifeDistributionD(sqrt(doubleJackknifeDistributionD::covariance(d,d))); });

  //Separate out the data in the desired fit range
  for(int q=0;q<10;q++)
    for(int t=0;t<args.Lt;t++)
      if(t >= args.tmin_k_op && t <= args.tsep_k_pi[tsep_k_pi_idx]-args.tmin_op_pi){
	amplitudeDataCoord coord(t,args.tsep_k_pi[tsep_k_pi_idx]);
	A0_fit[q].push_back(coord, A0_full_srcavg_j({q,t}));
	sigma_fit[q].push_back(sigma_j({q,t}));
      }
}

#endif

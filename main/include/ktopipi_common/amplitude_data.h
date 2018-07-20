#ifndef _KTOPIPI_FIT_AMPLITUDE_DATA_H_
#define _KTOPIPI_FIT_AMPLITUDE_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include "read_data.h"
#include "utils.h"
#include "compute_amplitude.h"

CPSFIT_START_NAMESPACE

struct readKtoPiPiDataOptions{
  bool load_data_checkpoint;
  std::string load_data_checkpoint_stub;
  
  bool save_data_checkpoint;
  std::string save_data_checkpoint_stub;

  bool use_symmetric_quark_momenta;
  std::string symmetric_quark_momenta_figure_file_extension;

  readKtoPiPiDataOptions(): load_data_checkpoint(false), save_data_checkpoint(false), use_symmetric_quark_momenta(false){}

  template<typename T>
  void import(const T &from){
#define C(X) X = from.X
    C(load_data_checkpoint);
    C(load_data_checkpoint_stub);
    C(save_data_checkpoint);
    C(save_data_checkpoint_stub);
    C(use_symmetric_quark_momenta);
    C(symmetric_quark_momenta_figure_file_extension);
#undef C
  }
};


template<typename Resampled, typename Raw, typename Resampler>
struct resampleFunctorGeneral{
  const Resampler &resampler;
  
  resampleFunctorGeneral(const Resampler &_resampler): resampler(_resampler){}

  inline Resampled operator()(int const* coord,const Raw &from) const{
    Resampled o(from.size());
    resampler.resample(o,from);
    return o;
  }
};

NumericTensor<rawDataDistributionD,1> getA2projectedBubble(const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan,
							   const int Lt, const int tsep_pipi, 
							   const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
  NumericTensor<rawDataDistributionD,1> bubble;
  if(opt.load_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << opt.load_data_checkpoint_stub << "_bubble.hdf5";
    std::cout << "Loading checkpoint data for bubble from " << file.str() << std::endl;
    HDF5reader rd(file.str());
    readBubble(rd,bubble,"bubble");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }else{
    bubble = readA2projectedBubble(traj_start,traj_inc,traj_lessthan,tsep_pipi,Lt,data_dir,opt.use_symmetric_quark_momenta);
  }
  if(opt.save_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << opt.save_data_checkpoint_stub << "_bubble.hdf5";
    std::cout << "Saving checkpoint data for bubble to " << file.str() << std::endl;
    HDF5writer wr(file.str());
    writeBubble(wr,bubble,"bubble");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }
  return bubble;
}



//Structure for containing raw and resampled A2-projected bubble data
struct BubbleData{
  NumericTensor<rawDataDistributionD,1> bubble;
  NumericTensor<rawDataDistributionD,1> bubble_binned;
  NumericTensor<jackknifeDistributionD,1> bubble_j;
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj;

  BubbleData(){}

  template<typename Resampler>
  BubbleData(const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size,
	     const int Lt, const int tsep_pipi, 
	     const Resampler &resampler, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    bubble = getA2projectedBubble(data_dir,traj_start,traj_inc,traj_lessthan,
				  Lt,tsep_pipi,opt);
    bubble_binned = bin(bubble,bin_size);
    bubble_j = bubble_binned.transform(resampleFunctorGeneral<jackknifeDistributionD,rawDataDistributionD,Resampler>(resampler));
    bubble_dj = bubble_binned.transform(resampleFunctorGeneral<doubleJackknifeDistributionD,rawDataDistributionD,Resampler>(resampler));
  }
};
template<typename DistributionType> struct getResampledBubble{};
template<> struct getResampledBubble<jackknifeDistributionD>{ static inline const NumericTensor<jackknifeDistributionD,1> &get(const BubbleData &bubble_data){ return bubble_data.bubble_j; }  };
template<> struct getResampledBubble<doubleJackknifeDistributionD>{ static inline const NumericTensor<doubleJackknifeDistributionD,1> &get(const BubbleData &bubble_data){ return bubble_data.bubble_dj; }  };


//Compute and store the raw amplitude data and data necessary for mix and vacuum subtractions
class RawKtoPiPiData{
public:
  //Type data
  IndexedContainer<NumericTensor<rawDataDistributionD,3>, 4, 1> A0_alltK; //(type idx){q,tK,t}

  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_nobub; //{q,tK,t}
  
  //Mix data
  IndexedContainer<NumericTensor<rawDataDistributionD,2>, 2, 3> mix_alltK; //(mix idx){tK,t}

  NumericTensor<rawDataDistributionD,2> mix4_alltK_nobub;  

  //Info on non-zero timeslices
  IndexedContainer<std::vector<int>, 4, 1> nonzerotK; //(type idx)

private:
  void getTypeData(IndexedContainer<type1234Data, 4, 1> &type_data, const int tsep_k_pi, 
		   const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan,
		   const int Lt, const int tsep_pipi, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    //Read the data
    if(opt.load_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.load_data_checkpoint_stub << "_tsepkpi" << tsep_k_pi << ".hdf5";
      std::cout << "Loading checkpoint data for tsep_k_pi = " << tsep_k_pi << " from " << file.str() << std::endl;
      HDF5reader rd(file.str());
      for(int i=1;i<=4;i++) read(rd,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }else{
      for(int i=1;i<=4;i++) type_data(i) = readType(i, traj_start, traj_inc, traj_lessthan, tsep_k_pi, tsep_pipi, Lt, 
						    data_dir, opt.use_symmetric_quark_momenta, opt.symmetric_quark_momenta_figure_file_extension);
    }

    //Write checkpoint if necessary
    if(opt.save_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.save_data_checkpoint_stub << "_tsepkpi" << tsep_k_pi << ".hdf5";
      std::cout << "Saving checkpoint data for tsep_k_pi = " << tsep_k_pi << " to " << file.str() << std::endl;
      HDF5writer wr(file.str());
      for(int i=1;i<=4;i++) write(wr,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }
  }
 
public:
  RawKtoPiPiData(){}

  RawKtoPiPiData(const int tsep_k_pi, const BubbleData &bubble_data, 
		 const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		 const int Lt, const int tsep_pipi, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    std::cout << "Reading data with tsep_k_pi=" << tsep_k_pi << std::endl;

    IndexedContainer<type1234Data, 4, 1> type_data;
    getTypeData(type_data, tsep_k_pi, data_dir, traj_start, traj_inc, traj_lessthan,
		Lt, tsep_pipi, opt);
    
    for(int i=1;i<=4;i++){
      nonzerotK(i) = type_data(i).getNonZeroKaonTimeslices();
      std::cout << "Type" << i << " data non-zero for tK=" << nonzerotK(i) << std::endl;
    }
    
    //Get the type1-4 and mix3/mix4 components of the data.
    //Note that we have not yet multiplied the type4/mix4 data by the pipi bubble, hence these  are just the K->vac component which are needed for computing alpha
    std::cout << "Computing raw data diagrams prior to multiplying by bubble\n";
    for(int i=1;i<=3;i++) A0_alltK(i) = computeAmplitudeType<computeAmplitudeAlltKtensorControls>(i,type_data(i)); //[Qidx][tK][t]
    A0_type4_alltK_nobub = computeAmplitudeType<computeAmplitudeAlltKtensorControls>(4,type_data(4)); //[Qidx][tK][t]
    
    mix_alltK(3) = NumericTensor<rawDataDistributionD,2>({Lt,Lt}); //[tK][t]
    mix4_alltK_nobub = NumericTensor<rawDataDistributionD,2>({Lt,Lt}); //[tK][t]
    for(int tK=0;tK<Lt;tK++)
      for(int t=0;t<Lt;t++){
	mix_alltK(3)({tK,t}) = type_data(3)(tK,t).mix();
	mix4_alltK_nobub({tK,t}) = type_data(4)(tK,t).mix();
      }
    
    printMem("Pre typedata free");
    
    //Data is no longer needed, so free it
    for(int i=1;i<=4;i++) type_data(i).freeData();
    
    printMem("Post typedata free");
    
    //Compute the type4/mix4 data with the bubble included
    std::cout << "Computing raw type4/mix4 data with bubble included" << std::endl;
    A0_alltK(4) = NumericTensor<rawDataDistributionD,3>({10,Lt,Lt});
    mix_alltK(4) = NumericTensor<rawDataDistributionD,2>({Lt,Lt});
    for(int tK=0;tK<Lt;tK++) for(int t=0;t<Lt;t++){
	int tB = (tK + tsep_k_pi) % Lt;
	mix_alltK(4)({tK,t}) = mix4_alltK_nobub({tK,t})*bubble_data.bubble(&tB);
	for(int q=0;q<10;q++)
	  A0_alltK(4)({q,tK,t}) = A0_type4_alltK_nobub({q,tK,t})*bubble_data.bubble(&tB);
      }
    
    //Bin everything we are going to use henceforth
    std::cout << "Binning data\n";
    for(int i=1;i<=4;i++) A0_alltK(i) = bin(A0_alltK(i), bin_size);
    for(int i=3;i<=4;i++) mix_alltK(i) = bin(mix_alltK(i), bin_size);
    A0_type4_alltK_nobub = bin(A0_type4_alltK_nobub, bin_size);
    mix4_alltK_nobub = bin(mix4_alltK_nobub, bin_size);  
    
    std::cout << "Finished reading raw data with tsep_k_pi=" << tsep_k_pi << std::endl;
  }
};

CPSFIT_END_NAMESPACE

#endif

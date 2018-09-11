#ifndef _KTOSIGMA_FIT_AMPLITUDE_DATA_H_
#define _KTOSIGMA_FIT_AMPLITUDE_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include "read_data_ktosigma.h"
#include "utils.h"
#include "compute_amplitude_ktosigma.h"

CPSFIT_START_NAMESPACE

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
  template<typename ReadPolicy>
  void getTypeData(IndexedContainer<type1234Data, 3, 2> &type_data, const int tsep_k_sigma, 
		   const int Lt, const ReadPolicy &rp, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
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
      for(int i=2;i<=4;i++) type_data(i) = readKtoSigmaType(i, tsep_k_sigma, Lt, rp);
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

    template<typename ReadPolicy>
    void getAllData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		    const int bin_size, const int Lt, const ReadPolicy &rp,
		    const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    std::cout << "Reading K->sigma data with tsep_k_sigma=" << tsep_k_sigma << std::endl;

    IndexedContainer<type1234Data, 3, 2> type_data;
    getTypeData(type_data, tsep_k_sigma, Lt, rp, opt);
    
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
 
public:
  RawKtoSigmaData(){}

  template<typename ReadPolicy>
  RawKtoSigmaData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		  const int bin_size, const int Lt, const ReadPolicy &rp,
		  const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    getAllData(tsep_k_sigma, bubble_data, bin_size, Lt, rp, opt);
  }

  RawKtoSigmaData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		  const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
		  const int traj_start, const int traj_inc, const int traj_lessthan, const int bin_size, 
		  const int Lt, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    KtoSigmaFilenamePolicyGen fp(data_file_fmt[0], data_file_fmt[1], data_file_fmt[2]);
    BasicKtoSigmaReadPolicy<KtoSigmaFilenamePolicyGen> rp(data_dir, traj_start, traj_inc, traj_lessthan, fp);
    getAllData(tsep_k_sigma, bubble_data, bin_size, Lt, rp, opt);
  }

};


CPSFIT_END_NAMESPACE

#endif
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

  ProjectedSigmaBubbleData(){}

  ProjectedSigmaBubbleData(const std::string &data_dir, const std::string &file_fmt, \
		  const int traj_start, const int traj_inc, const int traj_lessthan,
		  const int Lt, const std::vector<std::pair<threeMomentum, double> > &bubble_quarkmom_proj,
		  const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    bubble = getProjectedSigmaBubble(data_dir,file_fmt,
				     traj_start,traj_inc,traj_lessthan,
				     Lt, bubble_quarkmom_proj, opt);
  }

  template<typename DistributionType, typename Resampler>
  inline NumericTensor<DistributionType,1> binResample(const int bin_size, const Resampler &resampler) const{
    NumericTensor<DistributionType,1> out({bubble.size(0)}, 
					  [&](const int* t){ 
					    DistributionType r; resampler.resample(r, bubble(t).bin(bin_size)); return r;
					  }
					  );
    return out;
  }
};
void write(HDF5writer &writer, const ProjectedSigmaBubbleData &value, const std::string &tag){
  writeProjectedBubble(writer,value.bubble,tag);
}
void read(HDF5reader &reader, ProjectedSigmaBubbleData &value, const std::string &tag){
  readProjectedBubble(reader,value.bubble,tag);
}

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

  GENERATE_HDF5_SERIALIZE_METHOD( (A0_alltK)(A0_type4_alltK_nobub)(mix_alltK)(mix4_alltK_nobub) )

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
      for(int i=2;i<=4;i++) CPSfit::read(rd,type_data(i),stringize("type%d",i));
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
      for(int i=2;i<=4;i++) CPSfit::write(wr,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }
  }

  template<typename ReadPolicy>
  void getAllData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		  const int Lt, const ReadPolicy &rp,
		  const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

    std::cout << "Reading K->sigma data with tsep_k_sigma=" << tsep_k_sigma << std::endl;
    
    std::vector<int> nonzero_tK(Lt); for(int t=0;t<Lt;t++) nonzero_tK[t] = t;

    IndexedContainer<type1234Data, 3, 2> type_data;
    getTypeData(type_data, tsep_k_sigma, Lt, rp, opt);
    
    //Get the type1-4 and mix3/mix4 components of the data.
    //Note that we have not yet multiplied the type4/mix4 data by the pipi bubble, hence these  are just the K->vac component which are needed for computing alpha
    std::cout << "Computing raw data diagrams prior to multiplying by bubble\n";
    typedef computeAmplitudeAlltKtensorControls::inputType input;
    for(int i=2;i<=3;i++) A0_alltK(i) = computeKtoSigmaAmplitudeType<computeAmplitudeAlltKtensorControls>(i, input(type_data(i), nonzero_tK) ); //[Qidx][tK][t]
    A0_type4_alltK_nobub = computeKtoSigmaAmplitudeType<computeAmplitudeAlltKtensorControls>(4, input(type_data(4), nonzero_tK) ); //[Qidx][tK][t]

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
    
    std::cout << "Finished reading raw K->sigma data with tsep_k_sigma=" << tsep_k_sigma << std::endl;
  }
 
public:
  RawKtoSigmaData(){}

  template<typename ReadPolicy>
  RawKtoSigmaData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		  const int Lt, const ReadPolicy &rp,
		  const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    getAllData(tsep_k_sigma, bubble_data, Lt, rp, opt);
  }

  RawKtoSigmaData(const int tsep_k_sigma, const ProjectedSigmaBubbleData &bubble_data, 
		  const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
		  const int traj_start, const int traj_inc, const int traj_lessthan, 
		  const int Lt, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    KtoSigmaFilenamePolicyGen fp(data_file_fmt[0], data_file_fmt[1], data_file_fmt[2]);
    BasicKtoSigmaReadPolicy<KtoSigmaFilenamePolicyGen> rp(data_dir, traj_start, traj_inc, traj_lessthan, fp);
    getAllData(tsep_k_sigma, bubble_data, Lt, rp, opt);
  }

};

GENERATE_HDF5_SERIALIZE_FUNC(RawKtoSigmaData);


CPSFIT_END_NAMESPACE

#endif

#ifndef _KTOPIPI_FIT_AMPLITUDE_DATA_H_
#define _KTOPIPI_FIT_AMPLITUDE_DATA_H_

#include<config.h>
#include<utils/macros.h>

#include "read_data_ktopipi.h"
#include "utils.h"
#include "compute_amplitude_ktopipi.h"

CPSFIT_START_NAMESPACE

struct readKtoPiPiDataOptions{
  bool load_data_checkpoint;
  std::string load_data_checkpoint_stub;
  
  bool save_data_checkpoint;
  std::string save_data_checkpoint_stub;

  readKtoPiPiDataOptions(): load_data_checkpoint(false), save_data_checkpoint(false){}

  template<typename T>
  void import(const T &from){
#define C(X) X = from.X
    C(load_data_checkpoint);
    C(load_data_checkpoint_stub);
    C(save_data_checkpoint);
    C(save_data_checkpoint_stub);
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

NumericTensor<rawDataDistributionD,1> getProjectedBubble(const std::string &data_dir, const std::string &file_fmt,
							 const int traj_start, const int traj_inc, const int traj_lessthan, 
							 const int Lt, const int tsep_pipi, 
							 const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj, 
							 const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){

  NumericTensor<rawDataDistributionD,1> bubble;
  if(opt.load_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << opt.load_data_checkpoint_stub << "_bubble.hdf5";
    std::cout << "Loading checkpoint data for bubble from " << file.str() << std::endl;
    HDF5reader rd(file.str());
    readProjectedBubble(rd,bubble,"bubble");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }else{
    bubble = readProjectedPiPiBubble(data_dir, file_fmt, traj_start,traj_inc,traj_lessthan, Lt, tsep_pipi, bubble_pimom_proj);
  }
  if(opt.save_data_checkpoint){
#ifdef HAVE_HDF5
    std::ostringstream file; file << opt.save_data_checkpoint_stub << "_bubble.hdf5";
    std::cout << "Saving checkpoint data for bubble to " << file.str() << std::endl;
    HDF5writer wr(file.str());
    writeProjectedBubble(wr,bubble,"bubble");
#else
    error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
  }
  return bubble;
}



//Structure for containing raw and resampled A2-projected bubble data
struct ProjectedBubbleData{
  NumericTensor<rawDataDistributionD,1> bubble;

  ProjectedBubbleData(){}

  ProjectedBubbleData(const std::string &data_dir, const std::string &file_fmt, \
	     const int traj_start, const int traj_inc, const int traj_lessthan,
	     const int Lt, const int tsep_pipi, 
	     const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,
	     const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
  
    bubble = getProjectedBubble(data_dir,file_fmt,
				traj_start,traj_inc,traj_lessthan,
				Lt, tsep_pipi, bubble_pimom_proj, opt);
  }

  template<typename DistributionType, typename BinResampler>
  inline NumericTensor<DistributionType,1> binResample(const BinResampler &bin_resampler) const{
    NumericTensor<DistributionType,1> out({bubble.size(0)}, 
					  [&](const int* t){ 
					    DistributionType r; bin_resampler.binResample(r, bubble(t)); return r;
					  }
					  );
    return out;
  }
};
void write(HDF5writer &writer, const ProjectedBubbleData &value, const std::string &tag){
  writeProjectedBubble(writer,value.bubble,tag);
}
void read(HDF5reader &reader, ProjectedBubbleData &value, const std::string &tag){
  readProjectedBubble(reader,value.bubble,tag);
}



//Compute and store the raw amplitude data and data necessary for mix and vacuum subtractions
class RawKtoPiPiData{
public:
  //Type data
  IndexedContainer<NumericTensor<rawDataDistributionD,3>, 4, 1> A0_alltK; //(type idx){q,tK_idx,t}   tK_idx is the index of the non-zero tK for this type

  NumericTensor<rawDataDistributionD,3> A0_type4_alltK_nobub; //{q,tK_idx,t}
  
  //Mix data
  IndexedContainer<NumericTensor<rawDataDistributionD,2>, 2, 3> mix_alltK; //(mix idx){tK_idx,t}

  NumericTensor<rawDataDistributionD,2> mix4_alltK_nobub;   //{tK_idx, t}

  //Info on non-zero timeslices
  IndexedContainer<std::vector<int>, 4, 1> nonzerotK; //(type idx)

  GENERATE_HDF5_SERIALIZE_METHOD( (A0_alltK)(A0_type4_alltK_nobub)(mix_alltK)(mix4_alltK_nobub)(nonzerotK) );

private:
  template<typename ReadPolicy>
  static void getTypeData(IndexedContainer<type1234Data, 4, 1> &type_data, const int tsep_k_pi, 
			  const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
			  const int Lt, const int tsep_pipi, const ReadPolicy &rp, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    //Read the data
    if(opt.load_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.load_data_checkpoint_stub << "_tsepkpi" << tsep_k_pi << ".hdf5";
      std::cout << "Loading checkpoint data for tsep_k_pi = " << tsep_k_pi << " from " << file.str() << std::endl;
      HDF5reader rd(file.str());
      for(int i=1;i<=4;i++) CPSfit::read(rd,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }else{
      for(int i=1;i<=4;i++) type_data(i) = readKtoPiPiType(i, tsep_k_pi, tsep_pipi, Lt, type1_pimom_proj, rp);
    }

    //Write checkpoint if necessary
    if(opt.save_data_checkpoint){
#ifdef HAVE_HDF5
      std::ostringstream file; file << opt.save_data_checkpoint_stub << "_tsepkpi" << tsep_k_pi << ".hdf5";
      std::cout << "Saving checkpoint data for tsep_k_pi = " << tsep_k_pi << " to " << file.str() << std::endl;
      HDF5writer wr(file.str());
      for(int i=1;i<=4;i++) CPSfit::write(wr,type_data(i),stringize("type%d",i));
#else
      error_exit(std::cout << "Checkpointing of data requires HDF5\n");
#endif
    }
  }

  template<typename ReadPolicy>
  void getAllData(const int tsep_k_pi, const ProjectedBubbleData &bubble_data, 
	     const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
	     const int Lt, const int tsep_pipi, const ReadPolicy &rd,
	     const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    
    std::cout << "Reading data with tsep_k_pi=" << tsep_k_pi << std::endl;

    IndexedContainer<type1234Data, 4, 1> type_data;
    getTypeData(type_data, tsep_k_pi, type1_pimom_proj, Lt, tsep_pipi, rd, opt);
    
    int nt = tsep_k_pi + 1; //only compute on  0<=t<=tsep_k_pi

    int ntK[5];

    for(int i=1;i<=4;i++){
      nonzerotK(i) = type_data(i).getNonZeroKaonTimeslices();
      ntK[i] = nonzerotK(i).size();
      std::cout << "Type" << i << " data non-zero for tK=" << nonzerotK(i) << std::endl;
    }
    
    //Get the type1-4 and mix3/mix4 components of the data.
    //Note that we have not yet multiplied the type4/mix4 data by the pipi bubble, hence these  are just the K->vac component which are needed for computing alpha
    std::cout << "Computing raw data diagrams prior to multiplying by bubble" << std::endl;
    typedef computeAmplitudeAlltKtensorControls::inputType input;
    for(int i=1;i<=3;i++) A0_alltK(i) = computeAmplitudeType<computeAmplitudeAlltKtensorControls>(i, input(type_data(i), nonzerotK(i), tsep_k_pi ) ); //[Qidx][tK_idx][t]
    A0_type4_alltK_nobub = computeAmplitudeType<computeAmplitudeAlltKtensorControls>(4, input(type_data(4), nonzerotK(4), tsep_k_pi) ); //[Qidx][tK_idx][t]
    
    mix_alltK(3) = NumericTensor<rawDataDistributionD,2>({ntK[3],nt}); //[tK_idx][t]
    for(int tK_idx=0;tK_idx<ntK[3];tK_idx++)
      for(int t=0;t<nt;t++)
	mix_alltK(3)({tK_idx,t}) = type_data(3)(nonzerotK(3)[tK_idx],t).mix();


    mix4_alltK_nobub = NumericTensor<rawDataDistributionD,2>({ntK[4],nt}); //[tK_idx][t]
    for(int tK_idx=0;tK_idx<ntK[4];tK_idx++)
      for(int t=0;t<nt;t++)
	mix4_alltK_nobub({tK_idx,t}) = type_data(4)(nonzerotK(4)[tK_idx],t).mix();
  
    printMem("Pre typedata free");
    
    //Data is no longer needed, so free it
    for(int i=1;i<=4;i++) type_data(i).freeData();
    
    printMem("Post typedata free");
    
    //Compute the type4/mix4 data with the bubble included
    std::cout << "Computing raw type4/mix4 data with bubble included" << std::endl;
    A0_alltK(4) = NumericTensor<rawDataDistributionD,3>({10,ntK[4],nt});
    mix_alltK(4) = NumericTensor<rawDataDistributionD,2>({ntK[4],nt});
    for(int tK_idx=0;tK_idx<ntK[4];tK_idx++){
      int tK = nonzerotK(4)[tK_idx];
      for(int t=0;t<nt;t++){
	int tB = (tK + tsep_k_pi) % Lt;
	mix_alltK(4)({tK_idx,t}) = mix4_alltK_nobub({tK_idx,t})*bubble_data.bubble(&tB);
	for(int q=0;q<10;q++)
	  A0_alltK(4)({q,tK_idx,t}) = A0_type4_alltK_nobub({q,tK_idx,t})*bubble_data.bubble(&tB);
      }
    }
    
    std::cout << "Finished reading raw data with tsep_k_pi=" << tsep_k_pi << std::endl;
  }
 
public:
  RawKtoPiPiData(){}

  template<typename ReadPolicy>
  RawKtoPiPiData(const int tsep_k_pi, const ProjectedBubbleData &bubble_data, 
		 const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
		 const int Lt, const int tsep_pipi, const ReadPolicy &rd,
		 const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    getAllData(tsep_k_pi, bubble_data, type1_pimom_proj, Lt, tsep_pipi, rd, opt);
  }

  RawKtoPiPiData(const int tsep_k_pi, const ProjectedBubbleData &bubble_data, 
		 const std::string &data_dir, const std::vector<std::string> &data_file_fmt,
		 const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
		 const int traj_start, const int traj_inc, const int traj_lessthan,
		 const int Lt, const int tsep_pipi, const readKtoPiPiDataOptions &opt = readKtoPiPiDataOptions()){
    KtoPiPiFilenamePolicyGen fp(data_file_fmt[0], data_file_fmt[1], data_file_fmt[2], data_file_fmt[3]);
    BasicKtoPiPiReadPolicy<KtoPiPiFilenamePolicyGen> rp(data_dir, traj_start, traj_inc, traj_lessthan, fp);
    getAllData(tsep_k_pi, bubble_data, type1_pimom_proj, Lt, tsep_pipi, rp, opt);
  }
};

GENERATE_HDF5_SERIALIZE_FUNC(RawKtoPiPiData);




CPSFIT_END_NAMESPACE

#endif

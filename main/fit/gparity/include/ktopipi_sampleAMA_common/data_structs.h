#ifndef _KTOPIPI_SAMPLE_AMA_DATA_STRUCTS_H_
#define _KTOPIPI_SAMPLE_AMA_DATA_STRUCTS_H_

#include<config.h>
#include<utils/macros.h>

#include<ktopipi_common/amplitude_data_ktopipi.h>

SARLAC_START_NAMESPACE

struct readKtoPiPiDataSampleAMAoptions{
  readKtoPiPiDataOptions read_opts_sloppy_S;
  readKtoPiPiDataOptions read_opts_sloppy_C;
  readKtoPiPiDataOptions read_opts_exact_C;
};

enum SloppyExact {Sloppy, Exact};

struct allBubbleData{
  NumericTensor<rawDataDistributionD,1> bubble_sloppy_S;  
  NumericTensor<rawDataDistributionD,1> bubble_sloppy_S_binned;

  NumericTensor<rawDataDistributionD,1> bubble_sloppy_C;
  NumericTensor<rawDataDistributionD,1> bubble_sloppy_C_binned;

  NumericTensor<rawDataDistributionD,1> bubble_exact_C;
  NumericTensor<rawDataDistributionD,1> bubble_exact_C_binned;

  //full superjackknife
  NumericTensor<jackknifeDistributionD,1> bubble_j;
  NumericTensor<doubleJackknifeDistributionD,1> bubble_dj;

  allBubbleData(const std::string &bubble_file_fmt_sloppy, const std::string &bubble_file_fmt_exact,
		const std::vector<std::pair<threeMomentum, double> > &bubble_pimom_proj,
		const std::string &data_dir_S, const int traj_start_S, const int traj_lessthan_S,
		const std::string &data_dir_C, const int traj_start_C, const int traj_lessthan_C,
		const int traj_inc, const int bin_size, const int Lt, const int tsep_pipi,
		const sampleAMA_resamplers &resamplers, const readKtoPiPiDataSampleAMAoptions &opt = readKtoPiPiDataSampleAMAoptions())
    :bubble_j({Lt}), bubble_dj({Lt}){

    bubble_sloppy_S = getProjectedBubble(data_dir_S, bubble_file_fmt_sloppy, traj_start_S, traj_inc, traj_lessthan_S,
					 Lt, tsep_pipi, bubble_pimom_proj, opt.read_opts_sloppy_S);
      
    bubble_sloppy_C = getProjectedBubble(data_dir_C, bubble_file_fmt_sloppy, traj_start_C, traj_inc, traj_lessthan_C,
					 Lt, tsep_pipi, bubble_pimom_proj, opt.read_opts_sloppy_C);
      
    bubble_exact_C = getProjectedBubble(data_dir_C, bubble_file_fmt_exact, traj_start_C, traj_inc, traj_lessthan_C,
					Lt, tsep_pipi, bubble_pimom_proj, opt.read_opts_exact_C);

    bubble_sloppy_S_binned = bin(bubble_sloppy_S,bin_size);
    bubble_sloppy_C_binned = bin(bubble_sloppy_C,bin_size);
    bubble_exact_C_binned = bin(bubble_exact_C,bin_size);

    for(int t=0;t<Lt;t++){
      bubble_j(&t) = sampleAMAresampleCorrect<jackknifeDistributionD>(bubble_sloppy_S_binned(&t), bubble_sloppy_C_binned(&t), bubble_exact_C_binned(&t), 
								      resamplers.resampler_S, resamplers.resampler_C, stringize("Bubble(%d)",t));
      bubble_dj(&t) = sampleAMAresampleCorrect<doubleJackknifeDistributionD>(bubble_sloppy_S_binned(&t), bubble_sloppy_C_binned(&t), bubble_exact_C_binned(&t), 
									     resamplers.resampler_S, resamplers.resampler_C);
    }
  }
  
  ProjectedBubbleData extractUnbinnedBubble(const char ens, const SloppyExact se) const{
    ProjectedBubbleData out;
    if(ens == 'S'){
      assert(se == Sloppy);
      out.bubble = bubble_sloppy_S;
    }else{
      out.bubble = se == Sloppy ? bubble_sloppy_C : bubble_exact_C;
    }
    return out;
  }

  const NumericTensor<rawDataDistributionD,1> & rawBubble(const char ens, const SloppyExact se) const{
    if(ens == 'S'){
      assert(se == Sloppy);
      return bubble_sloppy_S_binned;
    }else{
      return se == Sloppy ? bubble_sloppy_C_binned : bubble_exact_C_binned;
    }
  }

};
template<typename DistributionType> struct getResampledBubbleSampleAMA{};
template<> struct getResampledBubbleSampleAMA<jackknifeDistributionD>{ 
  static inline const NumericTensor<jackknifeDistributionD,1> &get(const allBubbleData &bubble_data){ return bubble_data.bubble_j; }  
};
template<> struct getResampledBubbleSampleAMA<doubleJackknifeDistributionD>{
  static inline const NumericTensor<doubleJackknifeDistributionD,1> &get(const allBubbleData &bubble_data){ return bubble_data.bubble_dj; }  
};


struct allRawData{
  RawKtoPiPiData raw_sloppy_S;
  RawKtoPiPiData raw_sloppy_C;
  RawKtoPiPiData raw_exact_C;

  allRawData(const allBubbleData &bubble_data, const int tsep_k_pi,
	     const std::vector<std::string> &data_file_fmt_sloppy,  const std::vector<std::string> &data_file_fmt_exact,
	     const std::vector<std::pair<threeMomentum, double> > &type1_pimom_proj,
	     const std::string &data_dir_S, const int traj_start_S, const int traj_lessthan_S,
	     const std::string &data_dir_C, const int traj_start_C, const int traj_lessthan_C,
	     const int traj_inc, const int bin_size,
	     const int Lt, const int tsep_pipi, const readKtoPiPiDataSampleAMAoptions &opt = readKtoPiPiDataSampleAMAoptions()){

    ProjectedBubbleData bubble_sloppy_S = bubble_data.extractUnbinnedBubble('S',Sloppy);
    ProjectedBubbleData bubble_sloppy_C = bubble_data.extractUnbinnedBubble('C',Sloppy);
    ProjectedBubbleData bubble_exact_C = bubble_data.extractUnbinnedBubble('C',Exact);

    std::cout << "allRawData loading sloppy_S" << std::endl;
    raw_sloppy_S = RawKtoPiPiData(tsep_k_pi, bubble_sloppy_S, data_dir_S, 
				  data_file_fmt_sloppy, type1_pimom_proj, traj_start_S, traj_inc, traj_lessthan_S,
				  Lt, tsep_pipi, opt.read_opts_sloppy_S);

    std::cout << "allRawData loading sloppy_C" << std::endl;
    raw_sloppy_C = RawKtoPiPiData(tsep_k_pi, bubble_sloppy_C, data_dir_C,
				  data_file_fmt_sloppy, type1_pimom_proj, traj_start_C, traj_inc, traj_lessthan_C,
				  Lt, tsep_pipi, opt.read_opts_sloppy_C);

    std::cout << "allRawData loading exact_C" << std::endl;
    raw_exact_C = RawKtoPiPiData(tsep_k_pi, bubble_exact_C, data_dir_C,
				 data_file_fmt_exact, type1_pimom_proj, traj_start_C, traj_inc, traj_lessthan_C,
				 Lt, tsep_pipi, opt.read_opts_exact_C);
  }

  const RawKtoPiPiData &getRaw(const char ens, const SloppyExact se) const{
    if(ens == 'S'){
      assert(se == Sloppy);
      return raw_sloppy_S;
    }else{
      return se == Sloppy ? raw_sloppy_C : raw_exact_C;
    }
  }
};

SARLAC_END_NAMESPACE


#endif

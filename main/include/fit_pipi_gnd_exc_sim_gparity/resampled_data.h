#ifndef _PIPI_GND_EXC_SIM_FIT_RESAMPLED_DATA_H
#define _PIPI_GND_EXC_SIM_FIT_RESAMPLED_DATA_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

struct generateResampledDataOptions{
  bool load_combined_data;
  std::string load_combined_data_file;

  bool save_combined_data;
  std::string save_combined_data_file;

  generateResampledDataOptions(): load_combined_data(false), save_combined_data(false){}
};

void generateResampledData(jackknifeCorrelationFunction &j_data_gnd_gnd, doubleJackCorrelationFunction &dj_data_gnd_gnd, 
			   jackknifeCorrelationFunction &j_data_exc_exc, doubleJackCorrelationFunction &dj_data_exc_exc, 
			   jackknifeCorrelationFunction &j_data_gnd_exc, doubleJackCorrelationFunction &dj_data_gnd_exc,
			   const bubbleDataAllMomenta &raw_bubble_gnd_gnd, const bubbleDataAllMomenta &raw_bubble_exc_exc, const bubbleDataAllMomenta &raw_bubble_gnd_exc,
			   const rawCorrelationFunction &raw_data_gnd_gnd, const rawCorrelationFunction &raw_data_exc_exc, const rawCorrelationFunction &raw_data_gnd_exc,
			   const int tsep_pipi, const int bin_size, 
			   const bool do_vacuum_subtraction, const generateResampledDataOptions &opt = generateResampledDataOptions()){

  if(opt.load_combined_data){
    HDF5reader rd(opt.load_combined_data_file);
    read(rd, j_data_gnd_gnd, "j_data_gnd_gnd");
    read(rd, j_data_exc_exc, "j_data_exc_exc");
    read(rd, j_data_gnd_exc, "j_data_gnd_exc");
    read(rd, dj_data_gnd_gnd, "dj_data_gnd_gnd");
    read(rd, dj_data_exc_exc, "dj_data_exc_exc");
    read(rd, dj_data_gnd_exc, "dj_data_gnd_exc");
  }else{
    rawCorrelationFunction const* raw[3] = { &raw_data_gnd_gnd, &raw_data_exc_exc, &raw_data_gnd_exc };
    jackknifeCorrelationFunction *dsets_j[3] = {  &j_data_gnd_gnd, &j_data_exc_exc, &j_data_gnd_exc};
    doubleJackCorrelationFunction *dsets_dj[3] = {  &dj_data_gnd_gnd, &dj_data_exc_exc, &dj_data_gnd_exc};

    //Bin/resample
    for(int i=0;i<3;i++){
      *dsets_j[i] = binResample<jackknifeCorrelationFunction>(*raw[i], bin_size);
      *dsets_dj[i] = binResample<doubleJackCorrelationFunction>(*raw[i], bin_size);
    }

    //Compute vacuum subtractions
    if(do_vacuum_subtraction){
      std::cout << "Computing vacuum subtractions" << std::endl;
      bubbleDataAllMomenta const* bub[3] = { &raw_bubble_gnd_gnd, &raw_bubble_exc_exc, &raw_bubble_gnd_exc };
      static const PiPiProjector proj_src[3] = { PiPiProjector::A1momSet111, PiPiProjector::A1momSet311, PiPiProjector::A1momSet111 };
      static const PiPiProjector proj_snk[3] = { PiPiProjector::A1momSet111, PiPiProjector::A1momSet311, PiPiProjector::A1momSet311 };

      for(int i=0;i<3;i++){
	*dsets_j[i] = *dsets_j[i] - computePiPi2ptVacSub<jackknifeCorrelationFunction>(*bub[i], bin_size, tsep_pipi, proj_src[i], proj_snk[i]);
	*dsets_dj[i] = *dsets_dj[i] - computePiPi2ptVacSub<doubleJackCorrelationFunction>(*bub[i], bin_size, tsep_pipi, proj_src[i], proj_snk[i]);
      }
    }

    std::cout << "Folding data" << std::endl;
    for(int i=0;i<3;i++){
      *dsets_j[i] = fold(*dsets_j[i], 2*tsep_pipi);
      *dsets_dj[i] = fold(*dsets_dj[i], 2*tsep_pipi);
    }
  }

  if(opt.save_combined_data){
    HDF5writer wr(opt.save_combined_data_file);
    write(wr, j_data_gnd_gnd, "j_data_gnd_gnd");
    write(wr, j_data_exc_exc, "j_data_exc_exc");
    write(wr, j_data_gnd_exc, "j_data_gnd_exc");
    write(wr, dj_data_gnd_gnd, "dj_data_gnd_gnd");
    write(wr, dj_data_exc_exc, "dj_data_exc_exc");
    write(wr, dj_data_gnd_exc, "dj_data_gnd_exc");
  }

}

CPSFIT_END_NAMESPACE

#endif

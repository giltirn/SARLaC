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

void generateResampledData(doubleJackCorrelationFunction &dj_data_gnd_gnd, doubleJackCorrelationFunction &dj_data_exc_exc, doubleJackCorrelationFunction &dj_data_gnd_exc,
			   const bubbleDataAllMomenta &raw_bubble_gnd_gnd, const bubbleDataAllMomenta &raw_bubble_exc_exc, const bubbleDataAllMomenta &raw_bubble_gnd_exc,
			   const rawCorrelationFunction &raw_data_gnd_gnd, const rawCorrelationFunction &raw_data_exc_exc, const rawCorrelationFunction &raw_data_gnd_exc,
			   const int tsep_pipi, const int bin_size, 
			   const bool do_vacuum_subtraction, const generateResampledDataOptions &opt = generateResampledDataOptions()){

  static std::vector<threeMomentum> pion_mom = { {1,1,1}, {-1,-1,-1},
					  {-1,1,1}, {1,-1,-1},
					  {1,-1,1}, {-1,1,-1},
					  {1,1,-1}, {-1,-1,1},
  
					  {3,1,1}, {-3,-1,-1},
					  {-3,1,1}, {3,-1,-1},
					  {3,-1,1}, {-3,1,-1},
					  {3,1,-1}, {-3,-1,1},
					  
					  {1,3,1}, {-1,-3,-1},
					  {-1,3,1}, {1,-3,-1},
					  {1,-3,1}, {-1,3,-1},
					  {1,3,-1}, {-1,-3,1},
					  
					  {1,1,3}, {-1,-1,-3},
					  {-1,1,3}, {1,-1,-3},
					  {1,-1,3}, {-1,1,-3},
					  {1,1,-3}, {-1,-1,3}
  };

  if(opt.load_combined_data){
    HDF5reader rd(opt.load_combined_data_file);
    read(rd, dj_data_gnd_gnd, "dj_data_gnd_gnd");
    read(rd, dj_data_exc_exc, "dj_data_exc_exc");
    read(rd, dj_data_gnd_exc, "dj_data_gnd_exc");
  }else{
    dj_data_gnd_gnd = binDoubleJackknifeResample(raw_data_gnd_gnd, bin_size);
    dj_data_exc_exc = binDoubleJackknifeResample(raw_data_exc_exc, bin_size);
    dj_data_gnd_exc = binDoubleJackknifeResample(raw_data_gnd_exc, bin_size);
  
    //Compute vacuum subtractions
    if(do_vacuum_subtraction){
      std::cout << "Computing vacuum subtractions" << std::endl;
      doubleJackCorrelationFunction vac_sub_dj = computePiPi2ptVacSub(raw_bubble_gnd_gnd, bin_size, tsep_pipi, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
      dj_data_gnd_gnd = dj_data_gnd_gnd - vac_sub_dj;

      vac_sub_dj = computePiPi2ptVacSub(raw_bubble_exc_exc, bin_size, tsep_pipi, PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);
      dj_data_exc_exc = dj_data_exc_exc - vac_sub_dj;

      vac_sub_dj = computePiPi2ptVacSub(raw_bubble_gnd_exc, bin_size, tsep_pipi, PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);
      dj_data_gnd_exc = dj_data_gnd_exc - vac_sub_dj;
    }

    std::cout << "Folding data" << std::endl;
    //Fold data
    dj_data_gnd_gnd = foldPiPi2pt(dj_data_gnd_gnd, tsep_pipi);
    dj_data_exc_exc = foldPiPi2pt(dj_data_exc_exc, tsep_pipi);
    dj_data_gnd_exc = foldPiPi2pt(dj_data_gnd_exc, tsep_pipi);
  }

  if(opt.save_combined_data){
    HDF5writer wr(opt.save_combined_data_file);
    write(wr, dj_data_gnd_gnd, "dj_data_gnd_gnd");
    write(wr, dj_data_exc_exc, "dj_data_exc_exc");
    write(wr, dj_data_gnd_exc, "dj_data_gnd_exc");
  }

}

CPSFIT_END_NAMESPACE

#endif

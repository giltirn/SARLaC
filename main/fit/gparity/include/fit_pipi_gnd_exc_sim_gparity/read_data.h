#ifndef _PIPI_GND_EXC_SIM_FIT_READ_DATA_H
#define _PIPI_GND_EXC_SIM_FIT_READ_DATA_H

#include<config.h>
#include<utils/macros.h>

SARLAC_START_NAMESPACE

void readRawPiPiGndExcData(bubbleDataAllMomenta &raw_bubble_gnd_gnd, bubbleDataAllMomenta &raw_bubble_exc_exc, bubbleDataAllMomenta &raw_bubble_gnd_exc,
		rawDataCorrelationFunctionD &raw_data_gnd_gnd, rawDataCorrelationFunctionD &raw_data_exc_exc, rawDataCorrelationFunctionD &raw_data_gnd_exc,
		const std::string &data_dir, const std::string &figure_file_format, const std::string &bubble_file_format,
		const int Lt, const int tsep_pipi, const int tstep_pipi,
		const int traj_start, const int traj_inc, const int traj_lessthan){
  
  readPiPi2pt(raw_data_gnd_gnd, raw_bubble_gnd_gnd, data_dir, figure_file_format, bubble_file_format, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
	      PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
  
  readPiPi2pt(raw_data_exc_exc, raw_bubble_exc_exc, data_dir, figure_file_format, bubble_file_format, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
	      PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);
  
  readPiPi2pt(raw_data_gnd_exc, raw_bubble_gnd_exc, data_dir, figure_file_format, bubble_file_format, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
	      PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);
  
}

void loadPiPiGndExcCheckPoint(bubbleDataAllMomenta &raw_bubble_gnd_gnd, bubbleDataAllMomenta &raw_bubble_exc_exc, bubbleDataAllMomenta &raw_bubble_gnd_exc,
		    rawDataCorrelationFunctionD &raw_data_gnd_gnd, rawDataCorrelationFunctionD &raw_data_exc_exc, rawDataCorrelationFunctionD &raw_data_gnd_exc,
		    const std::string &file){
  HDF5reader rd(file);
  read(rd, raw_bubble_gnd_gnd, "raw_bubble_gnd_gnd");
  read(rd, raw_bubble_exc_exc, "raw_bubble_exc_exc");
  read(rd, raw_bubble_gnd_exc, "raw_bubble_gnd_exc");
  read(rd, raw_data_gnd_gnd, "raw_data_gnd_gnd");
  read(rd, raw_data_exc_exc, "raw_data_exc_exc");
  read(rd, raw_data_gnd_exc, "raw_data_gnd_exc");
}

void savePiPiGndExcCheckPoint(const bubbleDataAllMomenta &raw_bubble_gnd_gnd, const bubbleDataAllMomenta &raw_bubble_exc_exc, const bubbleDataAllMomenta &raw_bubble_gnd_exc,
		    const rawDataCorrelationFunctionD &raw_data_gnd_gnd, const rawDataCorrelationFunctionD &raw_data_exc_exc, const rawDataCorrelationFunctionD &raw_data_gnd_exc,
		    const std::string &file){
  HDF5writer wr(file);
  write(wr, raw_bubble_gnd_gnd, "raw_bubble_gnd_gnd");
  write(wr, raw_bubble_exc_exc, "raw_bubble_exc_exc");
  write(wr, raw_bubble_gnd_exc, "raw_bubble_gnd_exc");
  write(wr, raw_data_gnd_gnd, "raw_data_gnd_gnd");
  write(wr, raw_data_exc_exc, "raw_data_exc_exc");
  write(wr, raw_data_gnd_exc, "raw_data_gnd_exc");
}


struct readRawDataOptions{
  bool load_hdf5_data_checkpoint;
  std::string load_hdf5_data_checkpoint_stub;
  bool save_hdf5_data_checkpoint;
  std::string save_hdf5_data_checkpoint_stub;
  readRawDataOptions(): load_hdf5_data_checkpoint(false), save_hdf5_data_checkpoint(false){}
};

void getRawPiPiGndExcData(bubbleDataAllMomenta &raw_bubble_gnd_gnd, bubbleDataAllMomenta &raw_bubble_exc_exc, bubbleDataAllMomenta &raw_bubble_gnd_exc,
		rawDataCorrelationFunctionD &raw_data_gnd_gnd, rawDataCorrelationFunctionD &raw_data_exc_exc, rawDataCorrelationFunctionD &raw_data_gnd_exc,
		const std::string &data_dir, const std::string &figure_file_format, const std::string &bubble_file_format,
		const int Lt, const int tsep_pipi, const int tstep_pipi,
		const int traj_start, const int traj_inc, const int traj_lessthan,
		const readRawDataOptions &opt = readRawDataOptions()){

  if(opt.load_hdf5_data_checkpoint){
    loadPiPiGndExcCheckPoint(raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
		   raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc, opt.load_hdf5_data_checkpoint_stub);
  }else{
    readRawPiPiGndExcData(raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
		raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc,
		data_dir, figure_file_format, bubble_file_format,
		Lt, tsep_pipi, tstep_pipi, traj_start, traj_inc, traj_lessthan);
  }

  if(opt.save_hdf5_data_checkpoint){
    savePiPiGndExcCheckPoint(raw_bubble_gnd_gnd, raw_bubble_exc_exc, raw_bubble_gnd_exc,
		   raw_data_gnd_gnd, raw_data_exc_exc, raw_data_gnd_exc, opt.save_hdf5_data_checkpoint_stub);
  }
}

SARLAC_END_NAMESPACE

#endif

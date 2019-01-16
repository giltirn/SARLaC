#ifndef _FIT_PIPI_COMOVING_GEVP_RESAMPLED_DATA_H
#define _FIT_PIPI_COMOVING_GEVP_RESAMPLED_DATA_H

void saveCheckpoint(const std::map<threeMomentum, ResampledData<jackknifeCorrelationFunction> > &data_j,
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  write(wr, data_j, "j_data");
}


void loadCheckpoint(std::map<threeMomentum, ResampledData<jackknifeCorrelationFunction> > &data_j,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  read(rd, data_j, "j_data");
}

#endif

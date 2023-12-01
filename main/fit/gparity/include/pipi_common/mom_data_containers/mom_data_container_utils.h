#pragma once

#include<config.h>
#include<utils/macros.h>

#include "pipi_figure_mom_data_container.h"
#include "pipi_bubble_mom_data_container.h"
SARLAC_START_NAMESPACE

inline std::string checkpointFilename(const std::string &stub, const std::string &extra_descr){
  std::ostringstream filename;
  filename << stub;
  if(extra_descr.size() > 0) filename << "_" << extra_descr;
  filename << ".hdf5";
  return filename.str();
}

void saveRawDataCheckpoint(const figureDataAllMomenta &raw_data, const bubbleDataAllMomenta &raw_bubble_data, const std::string &filename_stub, const std::string &extra_descr){
  saveHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(filename_stub, extra_descr) );
}
void loadRawDataCheckpoint(figureDataAllMomenta &raw_data, bubbleDataAllMomenta &raw_bubble_data, const std::string &filename_stub, const std::string &extra_descr){
  loadHDF5checkpoint(raw_data, raw_bubble_data, checkpointFilename(filename_stub, extra_descr) );
}


SARLAC_END_NAMESPACE

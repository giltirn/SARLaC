#ifndef _PIPI_SIGMA_SIM_FIT_READ_DATA_H_
#define _PIPI_SIGMA_SIM_FIT_READ_DATA_H_


#include<config.h>
#include<utils/macros.h>

#include <common/common_defs.h>
#include <pipi_common/pipi_common.h>

SARLAC_START_NAMESPACE

void writeCheckpoint(const std::string &file, 
		     const rawDataCorrelationFunctionD &pipi_raw,
		     const rawDataCorrelationFunctionD &pipi_to_sigma_raw,
		     const rawDataCorrelationFunctionD &sigma2pt_raw,
		     const bubbleDataAllMomenta &pipi_self_data,
		     const sigmaSelfContraction &sigma_self_data){
  HDF5writer wr(file);
  write(wr, pipi_raw, "pipi_raw");
  write(wr, pipi_to_sigma_raw, "pipi_to_sigma_raw");
  write(wr, sigma2pt_raw, "sigma2pt_raw");
  write(wr, pipi_self_data,"pipi_self_data");
  write(wr, sigma_self_data, "sigma_self_data");
}
void readCheckpoint(rawDataCorrelationFunctionD &pipi_raw,
		    rawDataCorrelationFunctionD &pipi_to_sigma_raw,
		    rawDataCorrelationFunctionD &sigma2pt_raw,
		    bubbleDataAllMomenta &pipi_self_data,
		    sigmaSelfContraction &sigma_self_data,
		    const std::string &file){
  HDF5reader rd(file);
  read(rd, pipi_raw, "pipi_raw");
  read(rd, pipi_to_sigma_raw, "pipi_to_sigma_raw");
  read(rd, sigma2pt_raw, "sigma2pt_raw");
  read(rd, pipi_self_data, "pipi_self_data");
  read(rd, sigma_self_data, "sigma_self_data");
}

//compute_pipitosigma_disconn_ReRe = false was original strategy; setting to true gives a slight stat error improvement
void readData(rawDataCorrelationFunctionD &pipi_raw,
	      rawDataCorrelationFunctionD &pipi_to_sigma_raw,
	      rawDataCorrelationFunctionD &sigma2pt_raw,
	      bubbleDataAllMomenta &pipi_self_data,
	      sigmaSelfContraction &sigma_self_data,
	      const std::string &data_dir, 
	      const std::string &pipi2pt_figure_file_fmt, 
	      const std::string &sigma2pt_file_fmt, 
	      const std::string &pipitosigma_file_fmt, 
	      const std::string &pipi_bubble_file_fmt, 
	      const std::string &sigma_bubble_file_fmt,
	      const int tsep_pipi,
	      const int tstep_pipi2pt, int tstep_pipitosigma,
	      const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan,
	      bool compute_pipitosigma_disconn_ReRe){

  //pipi 2pt + pipi self (complex)
  bubbleDataAllMomentaZ pipi_self_data_Z;
  readPiPi2pt(pipi_raw, pipi_self_data_Z, data_dir, pipi2pt_figure_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi2pt, Lt, traj_start, traj_inc, traj_lessthan, 
	      PiPiProjector::A1momSet111, PiPiProjector::A1momSet111, 0);;
  pipi_self_data = reIm(pipi_self_data_Z, 0);

  //sigma 2pt + sigma self (complex)
  figureData sigma2pt_data;
  readSigmaSigma(sigma2pt_data, sigma2pt_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
  sigma2pt_raw = sourceAverage(sigma2pt_data);

  sigmaSelfContractionZ sigma_self_data_Z;
  readSigmaSelf(sigma_self_data_Z, sigma_bubble_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
  sigma_self_data = reIm(sigma_self_data_Z, 0);
  
  //pipi->sigma
  readReconstructPiPiToSigmaWithDisconnAllTsrcOptions opt;
  opt.compute_disconn_ReRe = compute_pipitosigma_disconn_ReRe;

  PiPiProjectorA1Basis111 proj_pipi;
  pipi_to_sigma_raw = readReconstructPiPiToSigmaWithDisconnAllTsrc(pipitosigma_file_fmt, data_dir, Lt, tstep_pipitosigma, proj_pipi, traj_start, traj_inc, traj_lessthan,
								   pipi_self_data_Z, sigma_self_data_Z,
								   opt);
}


struct rawData{
  rawDataCorrelationFunctionD pipi_raw, pipi_to_sigma_raw, sigma2pt_raw;
  bubbleDataAllMomenta pipi_self_data;
  sigmaSelfContraction sigma_self_data;

  void readDataFromOrigFiles(const std::string &data_dir, 
		const std::string &pipi2pt_figure_file_fmt, 
		const std::string &sigma2pt_file_fmt, 
		const std::string &pipitosigma_file_fmt, 
		const std::string &pipi_bubble_file_fmt, 
		const std::string &sigma_bubble_file_fmt,
		const int tsep_pipi,
		const int tstep_pipi2pt, int tstep_pipitosigma,
		const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan,
		bool compute_pipitosigma_disconn_ReRe){
        readData(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,
		 data_dir, 
		 pipi2pt_figure_file_fmt, sigma2pt_file_fmt, pipitosigma_file_fmt,
		 pipi_bubble_file_fmt, sigma_bubble_file_fmt,
		 tsep_pipi, tstep_pipi2pt, tstep_pipitosigma,
		 Lt, traj_start, traj_inc, traj_lessthan, compute_pipitosigma_disconn_ReRe);
  }
  void readDataFromCheckpoint(const std::string &file){
    readCheckpoint(pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data,file);
  }
  void saveCheckpoint(const std::string &file) const{
    writeCheckpoint(file,pipi_raw,pipi_to_sigma_raw,sigma2pt_raw,pipi_self_data,sigma_self_data);
  }
};



SARLAC_END_NAMESPACE

#endif

#ifndef _PIPI_SIGMA_SIM_FIT_READ_DATA_H_
#define _PIPI_SIGMA_SIM_FIT_READ_DATA_H_


#include<config.h>
#include<utils/macros.h>

#include <fit_sigmasigma_gparity/read_data.h>
#include <fit_pipitosigma_gparity/read_data.h>

CPSFIT_START_NAMESPACE

void readPiPi2pt(rawCorrelationFunction &pipi_raw, bubbleDataAllMomentaZ &raw_bubble_data,
		 const std::string &data_dir, 
		 const std::string &figure_file_fmt, const std::string &bubble_file_fmt, 
		 const int tsep_pipi, const int tstep_pipi, const int Lt,
		 const int traj_start, const int traj_inc, const int traj_lessthan, 
		 const std::vector<threeMomentum> &pion_mom,
		 const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  PiPiCorrelatorBasicSelector corr_select(proj_src, proj_snk,PiPiMomAllowed::All,{0,0,0});
  
  //Read C, D, R diagrams
  PiPiSymmetrySubsetFigureFileMapping ffn(data_dir, figure_file_fmt, traj_start, tsep_pipi, pion_mom, {0,0,0});
  figureDataAllMomenta raw_data;
  char figs[3] = {'C','D','R'};

  for(int f=0;f<3;f++){
    readFigure(raw_data, figs[f], data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, ffn, pion_mom, corr_select);
    zeroUnmeasuredSourceTimeslices(raw_data, figs[f], tstep_pipi);
  }

  bubbleFilenamePolicyGeneric bpsrc(bubble_file_fmt, {0,0,0}, Source);
  bubbleFilenamePolicyGeneric bpsnk(bubble_file_fmt, {0,0,0}, Sink);

  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, bpsrc, bpsnk, pion_mom, corr_select);

  //Compute V diagram using the real part of the bubble
  bubbleDataAllMomenta raw_bubble_data_real = reIm(raw_bubble_data, 0);

  computeV(raw_data, raw_bubble_data_real, tsep_pipi, pion_mom, corr_select);

  //Combine diagrams to construct raw correlator
  getRawPiPiCorrFunc(pipi_raw, raw_data, corr_select, 0, pion_mom, 1);
}



inline void readPiPi2pt(rawCorrelationFunction &pipi_raw, bubbleDataAllMomenta &raw_bubble_data,
			const std::string &data_dir, 
			const std::string &figure_file_fmt, const std::string &bubble_file_fmt, 
			const int tsep_pipi, const int tstep_pipi, const int Lt,
			const int traj_start, const int traj_inc, const int traj_lessthan, const std::vector<threeMomentum> &pion_mom,
			const PiPiProjector proj_src = PiPiProjector::A1, const PiPiProjector proj_snk = PiPiProjector::A1){
  bubbleDataAllMomentaZ raw_bubble_data_Z;
  readPiPi2pt(pipi_raw, raw_bubble_data_Z, data_dir, figure_file_fmt,  bubble_file_fmt, tsep_pipi, tstep_pipi, Lt,
	      traj_start, traj_inc, traj_lessthan, pion_mom, proj_src, proj_snk);
  raw_bubble_data = reIm(raw_bubble_data_Z, 0);
}



void writeCheckpoint(const std::string &file, 
		     const rawCorrelationFunction &pipi_raw,
		     const rawCorrelationFunction &pipi_to_sigma_raw,
		     const rawCorrelationFunction &sigma2pt_raw,
		     const bubbleDataAllMomenta &pipi_self_data,
		     const sigmaSelfContraction &sigma_self_data){
  HDF5writer wr(file);
  write(wr, pipi_raw, "pipi_raw");
  write(wr, pipi_to_sigma_raw, "pipi_to_sigma_raw");
  write(wr, sigma2pt_raw, "sigma2pt_raw");
  write(wr, pipi_self_data,"pipi_self_data");
  write(wr, sigma_self_data, "sigma_self_data");
}
void readCheckpoint(rawCorrelationFunction &pipi_raw,
		    rawCorrelationFunction &pipi_to_sigma_raw,
		    rawCorrelationFunction &sigma2pt_raw,
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
void readData(rawCorrelationFunction &pipi_raw,
	      rawCorrelationFunction &pipi_to_sigma_raw,
	      rawCorrelationFunction &sigma2pt_raw,
	      bubbleDataAllMomenta &pipi_self_data,
	      sigmaSelfContraction &sigma_self_data,
	      const std::string &data_dir, 
	      const std::string &pipi2pt_figure_file_fmt, 
	      const std::string &sigma2pt_file_fmt, 
	      const std::string &pipitosigma_file_fmt, 
	      const std::string &pipi_bubble_file_fmt, 
	      const std::string &sigma_bubble_file_fmt,
	      const int tsep_pipi, const std::vector<threeMomentum> &pion_mom,
	      const int tstep_pipi2pt, int tstep_pipitosigma,
	      const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan,
	      bool compute_pipitosigma_disconn_ReRe){

  //pipi 2pt + pipi self (complex)
  bubbleDataAllMomentaZ pipi_self_data_Z;
  readPiPi2pt(pipi_raw, pipi_self_data_Z, data_dir, pipi2pt_figure_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi2pt, Lt, traj_start, traj_inc, traj_lessthan, pion_mom);
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

  pipi_to_sigma_raw = readReconstructPiPiToSigmaWithDisconnAllTsrc(pipitosigma_file_fmt, data_dir, Lt, tstep_pipitosigma, traj_start, traj_inc, traj_lessthan,
								   pion_mom, pipi_self_data_Z, sigma_self_data_Z,
								   opt);
}

CPSFIT_END_NAMESPACE

#endif

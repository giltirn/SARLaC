#ifndef _PIPI_SIGMA_SIM_FIT_READ_DATA_H_
#define _PIPI_SIGMA_SIM_FIT_READ_DATA_H_

void readPiPi2pt(rawCorrelationFunction &pipi_raw, bubbleDataAllMomenta &raw_bubble_data,
		 const std::string &data_dir, const int tsep_pipi, const int tstep_pipi, const int Lt,
		 const int traj_start, const int traj_inc, const int traj_lessthan, const std::vector<threeMomentum> &pion_mom){
  PiPiProjectA1 proj;
  PiPiMomAllowAll allow;
  readFigureStationaryPolicy ffn(false);
  
  figureDataAllMomenta raw_data;
  char figs[3] = {'C','D','R'};

  for(int f=0;f<3;f++){
    readFigure(raw_data, figs[f], data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, ffn, pion_mom, proj, proj, allow);
    zeroUnmeasuredSourceTimeslices(raw_data, figs[f], tstep_pipi);
  }

  readBubbleStationaryPolicy bpsrc(false,Source);
  readBubbleStationaryPolicy bpsnk(false,Sink);

  readBubble(raw_bubble_data, data_dir, tsep_pipi, Lt, traj_start, traj_inc, traj_lessthan, bpsrc, bpsnk, pion_mom, proj, proj, allow);

  computeV(raw_data, raw_bubble_data, tsep_pipi, pion_mom, proj, proj, allow);

  getRawPiPiCorrFunc(pipi_raw, raw_data, raw_bubble_data, proj, proj, allow, 0, pion_mom, 1);
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

void readData(rawCorrelationFunction &pipi_raw,
	      rawCorrelationFunction &pipi_to_sigma_raw,
	      rawCorrelationFunction &sigma2pt_raw,
	      bubbleDataAllMomenta &pipi_self_data,
	      sigmaSelfContraction &sigma_self_data,
	      const std::vector<threeMomentum> &pion_mom,
	      const PiPiSigmaSimArgs &args){
  //pipi 2pt + pipi self
  readPiPi2pt(pipi_raw, pipi_self_data, args.data_dir, args.tsep_pipi, args.tstep_pipi2pt, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan, pion_mom);
  
  //pipi->sigma
  figureData pipi_to_sigma_data;
  readPiPiToSigma(pipi_to_sigma_data, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
  pipi_to_sigma_raw = sourceAverage(pipi_to_sigma_data);

  //sigma 2pt + sigma self
  figureData sigma2pt_data;
  readSigmaSigma(sigma2pt_data, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
  sigma2pt_raw = sourceAverage(sigma2pt_data);
  readSigmaSelf(sigma_self_data, args.data_dir, args.Lt, args.traj_start, args.traj_inc, args.traj_lessthan);
}



#endif

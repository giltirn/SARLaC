#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RAW_DATA_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RAW_DATA_H

class RawData{
  NumericSquareMatrix<bubbleDataAllMomenta> pipi_bubble;
  NumericSquareMatrix<bubbleDataAllMomentaZ> pipi_bubble_Z;
  NumericSquareMatrix<rawCorrelationFunction> correlators;
  sigmaSelfContraction sigma_self;
  sigmaSelfContractionZ sigma_self_Z;
public:
  //Real self-contractions
  inline bubbleDataAllMomenta & PiPiBubble(const Operator srcop, const Operator snkop){  return pipi_bubble((int)srcop, (int)snkop);  }
  inline const bubbleDataAllMomenta & PiPiBubble(const Operator srcop, const Operator snkop) const{  return pipi_bubble((int)srcop, (int)snkop);  }
  inline sigmaSelfContraction & SigmaBubble(){ return sigma_self; }
  inline const sigmaSelfContraction & SigmaBubble() const{ return sigma_self; }

  //Complex self-contractions
  inline bubbleDataAllMomentaZ & PiPiBubbleZ(const Operator srcop, const Operator snkop){ return pipi_bubble_Z((int)srcop, (int)snkop);  }
  inline const bubbleDataAllMomentaZ & PiPiBubbleZ(const Operator srcop, const Operator snkop) const{ return pipi_bubble_Z((int)srcop, (int)snkop);  }
  inline sigmaSelfContractionZ & SigmaBubbleZ(){ return sigma_self_Z; }
  inline const sigmaSelfContractionZ & SigmaBubbleZ() const{ return sigma_self_Z; }

  //Correlator data
  inline rawCorrelationFunction & correlator(const Operator srcop, const Operator snkop){ return correlators((int)srcop, (int)snkop); }
  inline const rawCorrelationFunction & correlator(const Operator srcop, const Operator snkop) const{ return correlators((int)srcop, (int)snkop); }
  
  RawData(): pipi_bubble(2), pipi_bubble_Z(2), correlators(3){}

  void read(const int Lt, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan,
	    const std::string &pipi_fig_file_fmt, const std::string &pipi_bubble_file_fmt, const int tsep_pipi, const int tstep_pipi,
	    const std::string &pipitosigma_file_fmt, const int tstep_pipi_to_sigma,
	    const std::string &sigma2pt_file_fmt, const std::string &sigma_bubble_file_fmt){

    //PiPi 2pt and PiPi bubble
    figureData::useFileCache() = true;
    readPiPi2pt(correlator(PiPiGnd,PiPiGnd), PiPiBubbleZ(PiPiGnd,PiPiGnd), data_dir, pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
		PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
    
    figureData::getFileCache().clear();
    
    readPiPi2pt(correlator(PiPiExc,PiPiExc), PiPiBubbleZ(PiPiExc,PiPiExc), data_dir, pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
		PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);
    
    figureData::getFileCache().clear();
    
    readPiPi2pt(correlator(PiPiGnd,PiPiExc), PiPiBubbleZ(PiPiGnd,PiPiExc), data_dir, pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, traj_start, traj_inc, traj_lessthan, 
		PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);
    
    figureData::getFileCache().clear();

    PiPiBubble(PiPiGnd,PiPiGnd) = reIm(PiPiBubbleZ(PiPiGnd,PiPiGnd),0);
    PiPiBubble(PiPiGnd,PiPiExc) = reIm(PiPiBubbleZ(PiPiGnd,PiPiExc),0);
    PiPiBubble(PiPiExc,PiPiExc) = reIm(PiPiBubbleZ(PiPiExc,PiPiExc),0);

    //Sigma 2pt and sigma bubble
    figureData sigma2pt_data;
    readSigmaSigma(sigma2pt_data, sigma2pt_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
    correlator(Sigma,Sigma) = sourceAverage(sigma2pt_data);

    readSigmaSelf(SigmaBubbleZ(), sigma_bubble_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
    SigmaBubble() = reIm(SigmaBubbleZ(), 0);
    
    
    //Pipi->sigma
    readReconstructPiPiToSigmaWithDisconnAllTsrcOptions opt;
    opt.compute_disconn_ReRe = true; 

    PiPiProjectorA1Basis111 proj_pipi_gnd;
    PiPiProjectorA1Basis311 proj_pipi_exc;
    correlator(PiPiGnd, Sigma) = readReconstructPiPiToSigmaWithDisconnAllTsrc(pipitosigma_file_fmt, data_dir, Lt, tstep_pipi_to_sigma, proj_pipi_gnd, traj_start, traj_inc, traj_lessthan,
									      PiPiBubbleZ(PiPiGnd,PiPiGnd), SigmaBubbleZ(), opt);

    correlator(PiPiExc, Sigma) = readReconstructPiPiToSigmaWithDisconnAllTsrc(pipitosigma_file_fmt, data_dir, Lt, tstep_pipi_to_sigma, proj_pipi_exc, traj_start, traj_inc, traj_lessthan,
									      PiPiBubbleZ(PiPiExc,PiPiExc), SigmaBubbleZ(), opt);


  }
};


#endif

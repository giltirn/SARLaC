#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RAW_DATA_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RAW_DATA_H

class RawData{
  NumericSquareMatrix<bubbleDataAllMomenta> pipi_bubble;
  NumericSquareMatrix<bubbleDataAllMomentaZ> pipi_bubble_Z;
  NumericSquareMatrix<rawCorrelationFunction> correlators;
  sigmaSelfContraction sigma_self;
  sigmaSelfContractionZ sigma_self_Z;
  std::set<std::pair<Operator,Operator> > contains;
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

  inline bool doOp(const Operator op, const std::vector<Operator> &incl_ops) const{ 
    return std::find(incl_ops.begin(), incl_ops.end(), op) != incl_ops.end();
  }

  bool haveData(const Operator opa, const Operator opb) const{ 
    return contains.find({opa,opb}) != contains.end();
  }

  void write(HDF5writer &wr, const std::string &nm) const{
    wr.enter(nm);
    CPSfit::write(wr, pipi_bubble, "pipi_bubble");
    CPSfit::write(wr, pipi_bubble_Z, "pipi_bubble_Z");
    CPSfit::write(wr, correlators, "correlators");
    CPSfit::write(wr, sigma_self, "sigma_self");
    CPSfit::write(wr, sigma_self_Z, "sigma_self_Z");
    CPSfit::write(wr, contains, "contains");
    wr.leave();
  }
  void read(HDF5reader &rd, const std::string &nm){
    rd.enter(nm);
    CPSfit::read(rd, pipi_bubble, "pipi_bubble");
    CPSfit::read(rd, pipi_bubble_Z, "pipi_bubble_Z");
    CPSfit::read(rd, correlators, "correlators");
    CPSfit::read(rd, sigma_self, "sigma_self");
    CPSfit::read(rd, sigma_self_Z, "sigma_self_Z");
    CPSfit::read(rd, contains, "contains");
    rd.leave();
  }

  void read(const int Lt, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan,
	    const std::string &pipi_fig_file_fmt, const std::string &pipi_bubble_file_fmt, const int tsep_pipi, const int tstep_pipi,
	    const std::string &pipitosigma_file_fmt, const int tstep_pipi_to_sigma,
	    const std::string &sigma2pt_file_fmt, const std::string &sigma_bubble_file_fmt,
	    const std::vector<Operator> &incl_ops){

    //PiPi 2pt and PiPi bubble
    figureData::useFileCache() = true;

    if(doOp(Operator::PiPiGnd, incl_ops)){
      readPiPi2pt(correlator(Operator::PiPiGnd,Operator::PiPiGnd), PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiGnd), data_dir, 
		  pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		  traj_start, traj_inc, traj_lessthan, 
		  PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
      figureData::getFileCache().clear();
      PiPiBubble(Operator::PiPiGnd,Operator::PiPiGnd) = reIm(PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiGnd),0);
      contains.insert({Operator::PiPiGnd,Operator::PiPiGnd});
    }    
    if(doOp(Operator::PiPiExc, incl_ops)){ 
      readPiPi2pt(correlator(Operator::PiPiExc,Operator::PiPiExc), PiPiBubbleZ(Operator::PiPiExc,Operator::PiPiExc), data_dir, 
		  pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		  traj_start, traj_inc, traj_lessthan, 
		  PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);
      
      figureData::getFileCache().clear();
      PiPiBubble(Operator::PiPiExc,Operator::PiPiExc) = reIm(PiPiBubbleZ(Operator::PiPiExc,Operator::PiPiExc),0);
      contains.insert({Operator::PiPiExc,Operator::PiPiExc});
    }
    if(doOp(Operator::PiPiGnd, incl_ops) && doOp(Operator::PiPiExc, incl_ops)){ 
      readPiPi2pt(correlator(Operator::PiPiGnd,Operator::PiPiExc), PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiExc), data_dir, 
		  pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		  traj_start, traj_inc, traj_lessthan, 
		  PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);
    
      figureData::getFileCache().clear();
      PiPiBubble(Operator::PiPiGnd,Operator::PiPiExc) = reIm(PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiExc),0);
      contains.insert({Operator::PiPiGnd,Operator::PiPiExc});
    }
    
    if(doOp(Operator::Sigma, incl_ops)){
      //Sigma 2pt and sigma bubble
      figureData sigma2pt_data;
      readSigmaSigma(sigma2pt_data, sigma2pt_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
      correlator(Operator::Sigma,Operator::Sigma) = sourceAverage(sigma2pt_data);
      readSigmaSelf(SigmaBubbleZ(), sigma_bubble_file_fmt, data_dir, Lt, traj_start, traj_inc, traj_lessthan);
      SigmaBubble() = reIm(SigmaBubbleZ(), 0);
      contains.insert({Operator::Sigma,Operator::Sigma});
    }

    PiPiProjectorA1Basis111 proj_pipi_gnd;
    PiPiProjectorA1Basis311 proj_pipi_exc;
    readReconstructPiPiToSigmaWithDisconnAllTsrcOptions opt;
    opt.compute_disconn_ReRe = true; 

    if(doOp(Operator::PiPiGnd, incl_ops) && doOp(Operator::Sigma, incl_ops)){
      correlator(Operator::PiPiGnd, Operator::Sigma) = readReconstructPiPiToSigmaWithDisconnAllTsrc(pipitosigma_file_fmt, data_dir, Lt, tstep_pipi_to_sigma, proj_pipi_gnd, 
										traj_start, traj_inc, traj_lessthan,
										PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiGnd), SigmaBubbleZ(), opt);
      contains.insert({Operator::PiPiGnd,Operator::Sigma});
    }
    if(doOp(Operator::PiPiExc, incl_ops) && doOp(Operator::Sigma, incl_ops)){
      correlator(Operator::PiPiExc, Operator::Sigma) = readReconstructPiPiToSigmaWithDisconnAllTsrc(pipitosigma_file_fmt, data_dir, Lt, tstep_pipi_to_sigma, proj_pipi_exc, 
										traj_start, traj_inc, traj_lessthan,
										PiPiBubbleZ(Operator::PiPiExc,Operator::PiPiExc), SigmaBubbleZ(), opt);
      contains.insert({Operator::PiPiExc,Operator::Sigma});
    }


  }
};


#endif

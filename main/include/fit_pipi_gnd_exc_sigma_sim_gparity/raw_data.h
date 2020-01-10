#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RAW_DATA_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RAW_DATA_H

#include<utils/utils/time.h>

class RawData{
  NumericSquareMatrix<bubbleDataAllMomenta> pipi_bubble;
  NumericSquareMatrix<bubbleDataAllMomentaZ> pipi_bubble_Z;
  NumericSquareMatrix<rawDataCorrelationFunctionD> correlators;
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
  inline rawDataCorrelationFunctionD & correlator(const Operator srcop, const Operator snkop){ return correlators((int)srcop, (int)snkop); }
  inline const rawDataCorrelationFunctionD & correlator(const Operator srcop, const Operator snkop) const{ return correlators((int)srcop, (int)snkop); }
  
  RawData(): pipi_bubble(2), pipi_bubble_Z(2), correlators(3){}

  inline bool doOp(const Operator op, const std::vector<Operator> &incl_ops) const{ 
    return std::find(incl_ops.begin(), incl_ops.end(), op) != incl_ops.end();
  }

  bool haveData(const Operator opa, const Operator opb) const{ 
    return contains.find({opa,opb}) != contains.end();
  }

  template<typename T>
  void iterateOverRawDistributions(const T &action){
    for(auto it=contains.begin(); it != contains.end(); it++){
      {//correlator      
	rawDataCorrelationFunctionD &raw = correlator(it->first, it->second);
	for(int i=0;i<raw.size();i++) action(raw.value(i));
      }
      if( (it->first == Operator::PiPiGnd || it->first == Operator::PiPiExc) &&
	  (it->second == Operator::PiPiGnd || it->second == Operator::PiPiExc) ){
	{//bubble
	  bubbleDataAllMomenta &bub = PiPiBubble(it->first, it->second);
	  for(auto p=bub.begin();p!=bub.end();p++)
	    for(int t=0;t<p->second.getLt();t++)
	      action(p->second.at(t));
	}
	{//bubbleZ
	  bubbleDataAllMomentaZ &bub = PiPiBubbleZ(it->first, it->second);
	  for(auto p=bub.begin();p!=bub.end();p++)
	    for(int t=0;t<p->second.getLt();t++)
	      action(p->second.at(t));
	}
      }
    }
    if(haveData(Operator::PiPiGnd,Operator::Sigma) || 
       haveData(Operator::PiPiExc,Operator::Sigma) ||
       haveData(Operator::Sigma,Operator::Sigma)
       ){
      {//sigma
	sigmaSelfContraction &bub = SigmaBubble();
	for(int t=0;t<bub.getLt();t++)
	  action(bub.at(t));
      }
      {//sigma Z
	sigmaSelfContractionZ &bub = SigmaBubbleZ();
	for(int t=0;t<bub.getLt();t++)
	  action(bub.at(t));
      }
    }
  }

  //For testing of bin size dependence we can randomly scramble the data to remove any real autocorrelations
  void scrambleSamples(){
    auto first_op_pair = *contains.begin();
    int nsample = correlator(first_op_pair.first, first_op_pair.second).value(0).size();
   
    if(!RNG.isInitialized()) RNG.initialize(1234);
    std::vector<int> reord(nsample);
    std::list<int> rem; 
    for(int i=0;i<nsample;i++) rem.push_back(i);
      
    for(int i=0;i<nsample;i++){
      int off = (int)uniformRandom<float>(0,rem.size());
      auto it = std::next(rem.begin(), off);
      reord[i] = *it;
      rem.erase(it);
    }
      
    //Check indices are unique
    std::cout << "Reordered samples: ";
    std::set<int> con;
    for(int i=0;i<nsample;i++){
      std::cout << reord[i] << " ";
      con.insert(reord[i]);
    }
    std::cout << std::endl;
    assert(con.size() == nsample); 
    
    iterateOverRawDistributions(
				[&](auto &v){ 
				  typedef typename std::decay<decltype(v)>::type DistributionType;
				  v = DistributionType(nsample, [&](const int s){ return v.sample(reord[s]); }); 
				}
      );
  }


  template<typename DistributionType>
  inline static DistributionType removeSamplesInRange(const DistributionType &raw, const int start, const int lessthan){
    int remove_num_samples = lessthan - start;
    assert(remove_num_samples < raw.size());
    DistributionType out(raw.size() - remove_num_samples);
    int s=0;
    for(int p=0;p<start;p++)
      out.sample(s++) = raw.sample(p);
    for(int p=lessthan;p<raw.size();p++)
      out.sample(s++) = raw.sample(p);
    return out;
  }

  void removeSamplesInRange(const int start, const int lessthan){
    std::cout << "RawData removing samples in range [" << start << ", " << lessthan << ")" << std::endl;
    iterateOverRawDistributions(
				[&](auto &v){ 
				  v = removeSamplesInRange(v,start,lessthan);				  
				}
				);
    std::cout << "RawData samples removed" << std::endl;
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
    timer time;

    rd.enter(nm);
    std::cout << "Reading pipi bubble" << std::endl;
    time.start();
    CPSfit::read(rd, pipi_bubble, "pipi_bubble");
    time.stop();
    std::cout << time.elapsed()/1e9 << "s" << std::endl;

    std::cout << "Reading pipi bubbleZ" << std::endl;
    time.start();
    CPSfit::read(rd, pipi_bubble_Z, "pipi_bubble_Z");
    time.stop();
    std::cout << time.elapsed()/1e9 << "s" << std::endl;


    std::cout << "Reading correlators" << std::endl;
    time.start();
    CPSfit::read(rd, correlators, "correlators");
    time.stop();
    std::cout << time.elapsed()/1e9 << "s" << std::endl;


    std::cout << "Reading sigma bubble" << std::endl;
    time.start();
    CPSfit::read(rd, sigma_self, "sigma_self");
    time.stop();
    std::cout << time.elapsed()/1e9 << "s" << std::endl;


    std::cout << "Reading sigma bubbleZ" << std::endl;
    time.start();
    CPSfit::read(rd, sigma_self_Z, "sigma_self_Z");
    time.stop();
    std::cout << time.elapsed()/1e9 << "s" << std::endl;

    std::cout << "Reading map" << std::endl;
    time.start();
    CPSfit::read(rd, contains, "contains");
    time.stop();
    std::cout << time.elapsed()/1e9 << "s" << std::endl;

    std::cout << "Reading complete" << std::endl;
    rd.leave();
  }

  void read(const int Lt, const std::string &data_dir, const int traj_start, const int traj_inc, const int traj_lessthan,
	    const std::string &pipi_fig_file_fmt, const std::string &pipi_bubble_file_fmt, const int tsep_pipi, const int tstep_pipi,
	    const std::string &pipitosigma_file_fmt, const int tstep_pipi_to_sigma,
	    const std::string &sigma2pt_file_fmt, const std::string &sigma_bubble_file_fmt,
	    const std::vector<Operator> &incl_ops, const int isospin){

    //PiPi 2pt and PiPi bubble
    if(doOp(Operator::PiPiGnd, incl_ops)){
      readPiPi2pt(correlator(Operator::PiPiGnd,Operator::PiPiGnd), PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiGnd), data_dir, 
		  pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		  traj_start, traj_inc, traj_lessthan, 
		  PiPiProjector::A1momSet111, PiPiProjector::A1momSet111, isospin);
      PiPiBubble(Operator::PiPiGnd,Operator::PiPiGnd) = reIm(PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiGnd),0);
      contains.insert({Operator::PiPiGnd,Operator::PiPiGnd});
    }    
    if(doOp(Operator::PiPiExc, incl_ops)){ 
      readPiPi2pt(correlator(Operator::PiPiExc,Operator::PiPiExc), PiPiBubbleZ(Operator::PiPiExc,Operator::PiPiExc), data_dir, 
		  pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		  traj_start, traj_inc, traj_lessthan, 
		  PiPiProjector::A1momSet311, PiPiProjector::A1momSet311, isospin);
      
      PiPiBubble(Operator::PiPiExc,Operator::PiPiExc) = reIm(PiPiBubbleZ(Operator::PiPiExc,Operator::PiPiExc),0);
      contains.insert({Operator::PiPiExc,Operator::PiPiExc});
    }
    if(doOp(Operator::PiPiGnd, incl_ops) && doOp(Operator::PiPiExc, incl_ops)){ 
      readPiPi2pt(correlator(Operator::PiPiGnd,Operator::PiPiExc), PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiExc), data_dir, 
		  pipi_fig_file_fmt, pipi_bubble_file_fmt, tsep_pipi, tstep_pipi, Lt, 
		  traj_start, traj_inc, traj_lessthan, 
		  PiPiProjector::A1momSet111, PiPiProjector::A1momSet311, isospin);
    
      PiPiBubble(Operator::PiPiGnd,Operator::PiPiExc) = reIm(PiPiBubbleZ(Operator::PiPiGnd,Operator::PiPiExc),0);
      contains.insert({Operator::PiPiGnd,Operator::PiPiExc});
    }
    
    if(doOp(Operator::Sigma, incl_ops)){
      if(isospin != 0) error_exit(std::cout << "Sigma operator is only applicable for I=0\n");

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

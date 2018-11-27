#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RESAMPLED_DATA_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RESAMPLED_DATA_H

//Multiplier of tsep_pipi around which the fold pivots
int foldOffsetMultiplier(const Operator a, const Operator b){
  if(  (a==Operator::PiPiGnd || a==Operator::PiPiExc) &&  (b==Operator::PiPiGnd || b==Operator::PiPiExc) ) return 2;
  else if( (a==Operator::PiPiGnd || a==Operator::PiPiExc) && b == Operator::Sigma ) return 1;
  else if( a == Operator::Sigma && b== Operator::Sigma ) return 0;
  assert(0);
}

template<typename resampledCorrelationFunctionType>
class ResampledData{
  NumericSquareMatrix<resampledCorrelationFunctionType> correlators;
  std::set<std::pair<Operator,Operator> > contains;
public:
  //Correlator data
  inline resampledCorrelationFunctionType & correlator(const Operator srcop, const Operator snkop){ return correlators((int)srcop, (int)snkop); }
  inline const resampledCorrelationFunctionType & correlator(const Operator srcop, const Operator snkop) const{ return correlators((int)srcop, (int)snkop); }
  
  ResampledData(): correlators(3){}

  bool haveData(const Operator opa, const Operator opb) const{ 
    return contains.find({opa,opb}) != contains.end();
  }

  void generatedResampledData(const RawData &raw_data, const int bin_size, const int Lt, const int tsep_pipi, const bool do_vacuum_subtraction = true){
    const static std::vector<std::pair<Operator,Operator> > rp = {  {Operator::PiPiGnd, Operator::PiPiGnd}, {Operator::PiPiGnd,Operator::PiPiExc}, {Operator::PiPiExc,Operator::PiPiExc}, {Operator::PiPiGnd,Operator::Sigma}, {Operator::PiPiExc,Operator::Sigma}, {Operator::Sigma,Operator::Sigma} };

    for(auto it=rp.begin();it!=rp.end();it++)
      if(raw_data.haveData(it->first, it->second))
	contains.insert({it->first, it->second});
    
    //Resample
    for(auto it=rp.begin();it!=rp.end();it++)
      if(raw_data.haveData(it->first, it->second))
	correlator(it->first, it->second) = binResample<resampledCorrelationFunctionType>(raw_data.correlator(it->first,it->second), bin_size);

    //Compute vacuum subtractions
    if(do_vacuum_subtraction){
      if(raw_data.haveData(Operator::PiPiGnd, Operator::PiPiGnd))
	correlator(Operator::PiPiGnd, Operator::PiPiGnd) = correlator(Operator::PiPiGnd, Operator::PiPiGnd) - 
	  computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(Operator::PiPiGnd, Operator::PiPiGnd), bin_size, tsep_pipi, 
								 PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);
      
      if(raw_data.haveData(Operator::PiPiGnd, Operator::PiPiExc))
	correlator(Operator::PiPiGnd, Operator::PiPiExc) = correlator(Operator::PiPiGnd, Operator::PiPiExc) - 
	  computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(Operator::PiPiGnd, Operator::PiPiExc), bin_size, tsep_pipi, PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);

      if(raw_data.haveData(Operator::PiPiExc, Operator::PiPiExc))
	correlator(Operator::PiPiExc, Operator::PiPiExc) = correlator(Operator::PiPiExc, Operator::PiPiExc) - 
	  computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(Operator::PiPiExc, Operator::PiPiExc), bin_size, tsep_pipi, PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);

      if(raw_data.haveData(Operator::PiPiGnd, Operator::Sigma))
	correlator(Operator::PiPiGnd, Operator::Sigma) = correlator(Operator::PiPiGnd, Operator::Sigma) -
	  computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), raw_data.PiPiBubble(Operator::PiPiGnd,Operator::PiPiGnd), PiPiProjector::A1momSet111, bin_size);

      if(raw_data.haveData(Operator::PiPiExc, Operator::Sigma))
	correlator(Operator::PiPiExc, Operator::Sigma) = correlator(Operator::PiPiExc, Operator::Sigma) -
	  computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), raw_data.PiPiBubble(Operator::PiPiExc,Operator::PiPiExc), PiPiProjector::A1momSet311, bin_size);

      if(raw_data.haveData(Operator::Sigma, Operator::Sigma))
	correlator(Operator::Sigma,Operator::Sigma) = correlator(Operator::Sigma,Operator::Sigma) -
	  computeSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), bin_size);
    }
 
    //Fold data
    for(auto it=rp.begin();it!=rp.end();it++)
      if(raw_data.haveData(it->first, it->second))
	correlator(it->first, it->second) = fold( 
						 correlator(it->first, it->second), 
						 foldOffsetMultiplier(it->first,it->second) * tsep_pipi
						  );
  }

  void write(HDF5writer &wr, const std::string &nm) const{
    wr.enter(nm);
    CPSfit::write(wr, correlators, "correlators");
    CPSfit::write(wr, contains, "contains");
    wr.leave();
  }
  void read(HDF5reader &rd, const std::string &nm){
    rd.enter(nm);
    CPSfit::read(rd, correlators, "correlators");
    CPSfit::read(rd, contains, "contains");
    rd.leave();
  }
};

void saveCheckpoint(const ResampledData<jackknifeCorrelationFunction> &data_j,
		    const ResampledData<doubleJackCorrelationFunction> &data_dj, 
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  data_j.write(wr, "j_data");
  data_dj.write(wr, "dj_data");
}


void loadCheckpoint(ResampledData<jackknifeCorrelationFunction> &data_j,
		    ResampledData<doubleJackCorrelationFunction> &data_dj, 		    
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  data_j.read(rd, "j_data");
  data_dj.read(rd, "dj_data");
}

#endif

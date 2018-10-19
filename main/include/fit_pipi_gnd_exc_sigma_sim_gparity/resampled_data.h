#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RESAMPLED_DATA_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RESAMPLED_DATA_H

//Multiplier of tsep_pipi around which the fold pivots
int foldOffsetMultiplier(const Operator a, const Operator b){
  if(  (a==PiPiGnd || a==PiPiExc) &&  (b==PiPiGnd || b==PiPiExc) ) return 2;
  else if( (a==PiPiGnd || a==PiPiExc) && b == Sigma ) return 1;
  else if( a == Sigma && b== Sigma ) return 0;
  assert(0);
}

template<typename resampledCorrelationFunctionType>
class ResampledData{
  NumericSquareMatrix<resampledCorrelationFunctionType> correlators;
public:
  //Correlator data
  inline resampledCorrelationFunctionType & correlator(const Operator srcop, const Operator snkop){ return correlators((int)srcop, (int)snkop); }
  inline const resampledCorrelationFunctionType & correlator(const Operator srcop, const Operator snkop) const{ return correlators((int)srcop, (int)snkop); }
  
  ResampledData(): correlators(3){}

  void generatedResampledData(const RawData &raw_data, const int bin_size, const int Lt, const int tsep_pipi, const bool do_vacuum_subtraction = true){
    const static std::vector<std::pair<Operator,Operator> > rp = {  {PiPiGnd, PiPiGnd}, {PiPiGnd,PiPiExc}, {PiPiExc,PiPiExc}, {PiPiGnd,Sigma}, {PiPiExc,Sigma}, {Sigma,Sigma} };
    
    //Resample
    for(auto it=rp.begin();it!=rp.end();it++)
      correlator(it->first, it->second) = binResample<resampledCorrelationFunctionType>(raw_data.correlator(it->first,it->second), bin_size);

    //Compute vacuum subtractions
    if(do_vacuum_subtraction){
      correlator(PiPiGnd, PiPiGnd) = correlator(PiPiGnd, PiPiGnd) - 
	computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(PiPiGnd, PiPiGnd), bin_size, tsep_pipi, PiPiProjector::A1momSet111, PiPiProjector::A1momSet111);

      correlator(PiPiGnd, PiPiExc) = correlator(PiPiGnd, PiPiExc) - 
	computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(PiPiGnd, PiPiExc), bin_size, tsep_pipi, PiPiProjector::A1momSet111, PiPiProjector::A1momSet311);

      correlator(PiPiExc, PiPiExc) = correlator(PiPiExc, PiPiExc) - 
	computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(PiPiExc, PiPiExc), bin_size, tsep_pipi, PiPiProjector::A1momSet311, PiPiProjector::A1momSet311);

      correlator(PiPiGnd, Sigma) = correlator(PiPiGnd, Sigma) -
	computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), raw_data.PiPiBubble(PiPiGnd,PiPiGnd), PiPiProjector::A1momSet111, bin_size);

      correlator(PiPiExc, Sigma) = correlator(PiPiExc, Sigma) -
	computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), raw_data.PiPiBubble(PiPiExc,PiPiExc), PiPiProjector::A1momSet311, bin_size);

      correlator(Sigma,Sigma) = correlator(Sigma,Sigma) -
	computeSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), bin_size);
    }
 
    //Fold data
    for(auto it=rp.begin();it!=rp.end();it++)
      correlator(it->first, it->second) = fold( 
					       correlator(it->first, it->second), 
					       foldOffsetMultiplier(it->first,it->second) * tsep_pipi
						);
  }




};

void saveCheckpoint(const ResampledData<jackknifeCorrelationFunction> &data_j,
		    const ResampledData<doubleJackCorrelationFunction> &data_dj, 
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  for(int i=0;i<3;i++)
    for(int j=i; j<3; j++){
      std::string post = stringize("_%s_%s", anyToStr( (Operator)i ).c_str(), anyToStr( (Operator)j ).c_str() );
      write(wr, data_j.correlator((Operator)i, (Operator)j), "j_data" + post);
      write(wr, data_dj.correlator((Operator)i, (Operator)j), "dj_data" + post);
    }
}


void loadCheckpoint(ResampledData<jackknifeCorrelationFunction> &data_j,
		    ResampledData<doubleJackCorrelationFunction> &data_dj, 		    
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  for(int i=0;i<3;i++)
    for(int j=i; j<3; j++){
      std::string post = stringize("_%s_%s", anyToStr( (Operator)i ).c_str(), anyToStr( (Operator)j ).c_str() );
      read(rd, data_j.correlator((Operator)i, (Operator)j), "j_data" + post);
      read(rd, data_dj.correlator((Operator)i, (Operator)j), "dj_data" + post);
    }
}



#endif

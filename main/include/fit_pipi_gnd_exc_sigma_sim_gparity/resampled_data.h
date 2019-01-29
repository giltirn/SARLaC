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
  int getNsample() const{
    auto it = contains.begin();
    return correlator(it->first, it->second).value(0).size();
  }

  static inline resampledCorrelationFunctionType timesliceAvg(const resampledCorrelationFunctionType &in){
    resampledCorrelationFunctionType out(in);
    auto v = in.value(0);
    for(int i=1;i<in.size();i++) v = v + in.value(i);
    v = v / double(in.size());
    for(int i=0;i<out.size();i++) out.value(i) = v;
    return out;
  }
  resampledCorrelationFunctionType computeVacSub(const RawData &raw_data, 
						 const Operator op1, const Operator op2,
						 const int bin_size, const int tsep_pipi,
						 const bool timeslice_avg_vac_sub = false){
    resampledCorrelationFunctionType v;
    if( (op1 == Operator::PiPiGnd || op1 == Operator::PiPiExc) &&
	(op2 == Operator::PiPiGnd || op2 == Operator::PiPiExc) ){
      PiPiProjector proj_src = op1 == Operator::PiPiGnd ?  PiPiProjector::A1momSet111 :  PiPiProjector::A1momSet311;
      PiPiProjector proj_snk = op2 == Operator::PiPiGnd ?  PiPiProjector::A1momSet111 :  PiPiProjector::A1momSet311;
      v = computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(op1,op2), bin_size, tsep_pipi, proj_src, proj_snk);
    }else if( (op1 == Operator::PiPiGnd || op1 == Operator::PiPiExc) &&
	      op2 == Operator::Sigma ){
      PiPiProjector proj_src = op1 == Operator::PiPiGnd ?  PiPiProjector::A1momSet111 :  PiPiProjector::A1momSet311;
      v = computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), raw_data.PiPiBubble(op1,op1), proj_src, bin_size);
    }else if( op1 == Operator::Sigma && op2 == Operator::Sigma ){
      v = computeSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), bin_size);
    }else assert(0);
    
    if(timeslice_avg_vac_sub) v = timesliceAvg(v);
    return v;
  }

  void generatedResampledData(const RawData &raw_data, const int bin_size, const int Lt, const int tsep_pipi, const bool do_vacuum_subtraction = true, const bool timeslice_avg_vac_sub = false){
    const static std::vector<std::pair<Operator,Operator> > rp = {  {Operator::PiPiGnd, Operator::PiPiGnd}, {Operator::PiPiGnd,Operator::PiPiExc}, {Operator::PiPiExc,Operator::PiPiExc}, {Operator::PiPiGnd,Operator::Sigma}, {Operator::PiPiExc,Operator::Sigma}, {Operator::Sigma,Operator::Sigma} };

    for(auto it=rp.begin();it!=rp.end();it++)
      if(raw_data.haveData(it->first, it->second)){
	contains.insert({it->first, it->second});

	auto &corr = correlator(it->first, it->second);

	//Resample	
	corr = binResample<resampledCorrelationFunctionType>(raw_data.correlator(it->first,it->second), bin_size);

	//Vacuum subtract
	if(do_vacuum_subtraction) corr = corr - 
				    computeVacSub(raw_data, it->first, it->second, bin_size, tsep_pipi, timeslice_avg_vac_sub);

	//Fold
	corr = fold( 
		    corr, 
		    foldOffsetMultiplier(it->first,it->second) * tsep_pipi
		     );

      } 
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


void saveCheckpoint(const ResampledData<jackknifeCorrelationFunction> &data_j,
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  data_j.write(wr, "j_data");
}


void loadCheckpoint(ResampledData<jackknifeCorrelationFunction> &data_j,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  data_j.read(rd, "j_data");
}


#endif

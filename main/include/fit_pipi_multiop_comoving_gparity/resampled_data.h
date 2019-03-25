#ifndef _FIT_PIPI_COMOVING_RESAMPLED_DATA_H
#define _FIT_PIPI_COMOVING_RESAMPLED_DATA_H

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

  void generatedResampledData(const std::vector<Operator> &ops, const RawData &raw_data, const int bin_size, const int isospin, const int Lt, const int tsep_pipi, const threeMomentum &p_tot, const bool do_vacuum_subtraction){
    if(isospin == 2) assert(do_vacuum_subtraction == false);

    static const std::vector< std::pair<Operator,Operator> > rp = { {Operator::PiPiComoveGnd, Operator::PiPiComoveGnd},
								    {Operator::PiPiComoveGnd, Operator::PiPiComoveExc1},
								    {Operator::PiPiComoveGnd, Operator::PiPiComoveExc2},
								    {Operator::PiPiComoveExc1, Operator::PiPiComoveExc1},
								    {Operator::PiPiComoveExc1, Operator::PiPiComoveExc2},
								    {Operator::PiPiComoveExc2, Operator::PiPiComoveExc2} };

    static const std::map<Operator, PiPiProjector> op_proj = { {Operator::PiPiComoveGnd, PiPiProjector::MovingSwaveGround},
							       {Operator::PiPiComoveExc1, PiPiProjector::MovingSwaveExc1},
							       {Operator::PiPiComoveExc2, PiPiProjector::MovingSwaveExc2} };
    

    for(auto it=rp.begin();it!=rp.end();it++){
      auto o1 = it->first; auto o2 = it->second;
      
      if(std::find(ops.begin(),ops.end(),o1) == ops.end() || std::find(ops.begin(),ops.end(),o2) == ops.end() ) continue;

      if(raw_data.haveData(o1,o2)){
	contains.insert({o1,o2});
	auto &corr = correlator(o1,o2);

	corr = binResample<resampledCorrelationFunctionType>(raw_data.correlator(o1,o2), bin_size);

	if(do_vacuum_subtraction) //technically this should not be needed!
	  corr = corr - computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(o1,o2), bin_size, tsep_pipi, p_tot,
									       op_proj.find(o1)->second , op_proj.find(o2)->second );
	corr = fold(corr, 2*tsep_pipi); 
      }
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

template<typename T>
void write(HDF5writer &wr, const ResampledData<T> &d, const std::string &nm){
  d.write(wr, nm);
}
template<typename T>
void read(HDF5reader &rd, ResampledData<T> &d, const std::string &nm){
  d.read(rd, nm);
}



void saveCheckpoint(const std::map<threeMomentum, ResampledData<jackknifeCorrelationFunction> > &data_j,
		    const std::map<threeMomentum, ResampledData<doubleJackCorrelationFunction> > &data_dj, 
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  write(wr, data_j, "j_data");
  write(wr, data_dj, "dj_data");
}


void loadCheckpoint(std::map<threeMomentum, ResampledData<jackknifeCorrelationFunction> > &data_j,
		    std::map<threeMomentum, ResampledData<doubleJackCorrelationFunction> > &data_dj,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  read(rd, data_j, "j_data");
  read(rd, data_dj, "dj_data");
}

#endif

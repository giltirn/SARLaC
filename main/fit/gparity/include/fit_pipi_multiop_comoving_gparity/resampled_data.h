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
  
  template<typename Resampler>
  void generatedResampledData(const std::vector<Operator> &ops, const RawData &raw_data, const Resampler &resampler, const int isospin, const int Lt, const int tsep_pipi, const threeMomentum &p_tot, const bool do_vacuum_subtraction){
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

	corr = binResample<resampledCorrelationFunctionType>(raw_data.correlator(o1,o2), resampler);

	if(do_vacuum_subtraction) //technically this should not be needed!
	  corr = corr - computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(o1,o2), resampler, tsep_pipi, p_tot,
									       op_proj.find(o1)->second , op_proj.find(o2)->second );
	corr = fold(corr, 2*tsep_pipi); 
      }
    }
  }
  inline void generatedResampledData(const std::vector<Operator> &ops, const RawData &raw_data, const int bin_size, const int isospin, const int Lt, const int tsep_pipi, const threeMomentum &p_tot, const bool do_vacuum_subtraction){
    basicBinResampler resampler(bin_size);
    generatedResampledData(ops, raw_data, resampler, isospin, Lt, tsep_pipi, p_tot, do_vacuum_subtraction);
  }


  void write(HDF5writer &wr, const std::string &nm) const{
    wr.enter(nm);
    SARLaC::write(wr, correlators, "correlators");
    SARLaC::write(wr, contains, "contains");
    wr.leave();
  }
  void read(HDF5reader &rd, const std::string &nm){
    rd.enter(nm);
    SARLaC::read(rd, correlators, "correlators");
    SARLaC::read(rd, contains, "contains");
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



void saveCheckpoint(const std::map<threeMomentum, ResampledData<jackknifeCorrelationFunctionD> > &data_j,
		    const std::map<threeMomentum, ResampledData<doubleJackknifeCorrelationFunctionD> > &data_dj, 
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  write(wr, data_j, "j_data");
  write(wr, data_dj, "dj_data");
}


void loadCheckpoint(std::map<threeMomentum, ResampledData<jackknifeCorrelationFunctionD> > &data_j,
		    std::map<threeMomentum, ResampledData<doubleJackknifeCorrelationFunctionD> > &data_dj,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  read(rd, data_j, "j_data");
  read(rd, data_dj, "dj_data");
}

void saveCheckpoint(const std::map<threeMomentum, ResampledData<bootstrapCorrelationFunctionD> > &data_b,
		    const std::map<threeMomentum, ResampledData<bootJackknifeCorrelationFunctionD> > &data_bj, 
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  write(wr, data_b, "b_data");
  write(wr, data_bj, "bj_data");
}


void loadCheckpoint(std::map<threeMomentum, ResampledData<bootstrapCorrelationFunctionD> > &data_b,
		    std::map<threeMomentum, ResampledData<bootJackknifeCorrelationFunctionD> > &data_bj,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  read(rd, data_b, "b_data");
  read(rd, data_bj, "bj_data");
}



#endif

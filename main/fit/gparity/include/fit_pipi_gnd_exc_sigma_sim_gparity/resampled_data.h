#ifndef _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RESAMPLED_DATA_H
#define _FIT_PIPI_GND_EXC_SIGMA_GPARITY_RESAMPLED_DATA_H

#include "filters.h"
#include "fitfunc.h"
#include <pipi_common/resampled_data.h>

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
  template<typename binResampler>
  static resampledCorrelationFunctionType computeVacSub(const RawData &raw_data, 
							const Operator op1, const Operator op2,
							const binResampler &resampler, const int tsep_pipi,
							const bool timeslice_avg_vac_sub = false){
    resampledCorrelationFunctionType v;
    if( (op1 == Operator::PiPiGnd || op1 == Operator::PiPiExc) &&
	(op2 == Operator::PiPiGnd || op2 == Operator::PiPiExc) ){
      PiPiProjector proj_src = op1 == Operator::PiPiGnd ?  PiPiProjector::A1momSet111 :  PiPiProjector::A1momSet311;
      PiPiProjector proj_snk = op2 == Operator::PiPiGnd ?  PiPiProjector::A1momSet111 :  PiPiProjector::A1momSet311;
      v = computePiPi2ptVacSub<resampledCorrelationFunctionType>(raw_data.PiPiBubble(op1,op2), resampler, tsep_pipi, proj_src, proj_snk);
    }else if( (op1 == Operator::PiPiGnd || op1 == Operator::PiPiExc) &&
	      op2 == Operator::Sigma ){
      PiPiProjector proj_src = op1 == Operator::PiPiGnd ?  PiPiProjector::A1momSet111 :  PiPiProjector::A1momSet311;
      v = computePiPiToSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), raw_data.PiPiBubble(op1,op1), proj_src, resampler);
    }else if( op1 == Operator::Sigma && op2 == Operator::Sigma ){
      v = computeSigmaVacSub<resampledCorrelationFunctionType>(raw_data.SigmaBubble(), resampler);
    }else assert(0);
    
    if(timeslice_avg_vac_sub) v = timesliceAvg(v);
    return v;
  }
  template<typename binResampler>
  void generatedResampledData(const RawData &raw_data, const binResampler &resampler, const int Lt, const int tsep_pipi, 
			      const bool do_vacuum_subtraction = true, const bool timeslice_avg_vac_sub = false, bool do_fold=true){

    const static std::vector<std::pair<Operator,Operator> > rp = {  
      {Operator::PiPiGnd, Operator::PiPiGnd}, {Operator::PiPiGnd,Operator::PiPiExc}, {Operator::PiPiExc,Operator::PiPiExc}, 
      {Operator::PiPiGnd,Operator::Sigma}, {Operator::PiPiExc,Operator::Sigma}, {Operator::Sigma,Operator::Sigma} 
    };

    for(auto it=rp.begin();it!=rp.end();it++)
      if(raw_data.haveData(it->first, it->second)){
  	contains.insert({it->first, it->second});

  	auto &corr = correlator(it->first, it->second);

  	//Resample	
  	corr = binResample<resampledCorrelationFunctionType>(raw_data.correlator(it->first,it->second), resampler);

  	//Vacuum subtract
  	if(do_vacuum_subtraction) corr = corr - 
  				    computeVacSub(raw_data, it->first, it->second, resampler, tsep_pipi, timeslice_avg_vac_sub);

  	//Fold
  	if(do_fold) corr = fold( 
				corr, 
				foldOffsetMultiplier(it->first,it->second) * tsep_pipi
				 );

      } 
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
inline void write(HDF5writer &wr, const ResampledData<T> &v, const std::string &nm){
  v.write(wr, nm);
}
template<typename T>
inline void read(HDF5reader &rd, ResampledData<T> &v, const std::string &nm){
  v.read(rd, nm);
}


void saveCheckpoint(const ResampledData<jackknifeCorrelationFunctionD> &data_j,
		    const ResampledData<doubleJackknifeCorrelationFunctionD> &data_dj, 
		    const ResampledData<blockDoubleJackknifeCorrelationFunctionD> &data_bdj, 
		    const bool do_dj, const bool do_bdj,
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  data_j.write(wr, "j_data");
  if(do_dj) data_dj.write(wr, "dj_data");
  if(do_bdj) data_bdj.write(wr, "bdj_data");
}


void loadCheckpoint(ResampledData<jackknifeCorrelationFunctionD> &data_j,
		    ResampledData<doubleJackknifeCorrelationFunctionD> &data_dj, 
		    ResampledData<blockDoubleJackknifeCorrelationFunctionD> &data_bdj,
		    const bool do_dj, const bool do_bdj,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  data_j.read(rd, "j_data");
  if(do_dj) data_dj.read(rd, "dj_data");
  if(do_bdj) data_bdj.read(rd, "bdj_data");
}



void saveCheckpoint(const ResampledData<jackknifeCorrelationFunctionD> &data_j,
		    const ResampledData<doubleJackknifeCorrelationFunctionD> &data_dj, 
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  data_j.write(wr, "j_data");
  data_dj.write(wr, "dj_data");
}


void loadCheckpoint(ResampledData<jackknifeCorrelationFunctionD> &data_j,
		    ResampledData<doubleJackknifeCorrelationFunctionD> &data_dj, 		    
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  data_j.read(rd, "j_data");
  data_dj.read(rd, "dj_data");
}


void saveCheckpoint(const ResampledData<jackknifeCorrelationFunctionD> &data_j,
		    const std::string &file){
  std::cout << "Saving data checkpoint to " << file << std::endl;
  HDF5writer wr(file);
  data_j.write(wr, "j_data");
}


void loadCheckpoint(ResampledData<jackknifeCorrelationFunctionD> &data_j,
		    const std::string &file){
  std::cout << "Reading data checkpoint from " << file << std::endl;
  HDF5reader rd(file);
  data_j.read(rd, "j_data");
}





struct DataDescr{
  Operator op1;
  Operator op2;
  int t;
  DataDescr(Operator op1, Operator op2, int t): op1(op1),op2(op2),t(t){}
};

template<typename CorrelationFunctionType> 
std::vector<DataDescr> getFitDataElemIdx(const ResampledData<CorrelationFunctionType> &data_j,
					 const std::vector<Operator> &ops, const int tsep_pipi, const int t_min, const int t_max,
					 const Filters &filters, const bool use_filters){
  std::vector<DataDescr> keep;
  
  for(int i=0;i<ops.size();i++){
    for(int j=i;j<ops.size();j++){
      std::ostringstream nm; nm << ops[i] << " " << ops[j];

      for(int t=t_min;t<=t_max;t++){
	bool skip = false;
	std::string reason;
	if(use_filters)
	  for(int f=0;f<filters.filters.size();f++)
	    if(filters.filters[f].filterOut(ops[i],ops[j],t,data_j.correlator(ops[i],ops[j]).value(t),&reason)){
	      skip = true;
	      std::cout << "Skipping " << nm.str() << " t=" << t << " as: " << reason << std::endl;
	    }
	if(skip)
	  continue;

	keep.push_back(DataDescr(ops[i],ops[j],t));
      }
    }
  }
  return keep;
}


template<typename DistributionType>
void filterData(correlationFunction<SimFitCoordGen,  DistributionType> &corr_comb,
		const ResampledData<correlationFunction<double,DistributionType> > &data,
		const std::vector<DataDescr> &keep, 
		const std::map< std::pair<Operator,Operator>, SubFitFuncParameterMap > &subfit_pmaps, const int tsep_pipi){

  for(int i=0;i<keep.size();i++){
    std::unordered_map<std::string, std::string> const* pmap = &subfit_pmaps.find({keep[i].op1,keep[i].op2})->second;    
    SimFitCoordGen coord(keep[i].t, pmap, foldOffsetMultiplier(keep[i].op1,keep[i].op2)*tsep_pipi);
    corr_comb.push_back(coord, data.correlator(keep[i].op1,keep[i].op2).value(keep[i].t));
  }
}


#endif

#ifndef _FIT_SIMPLE_READ_DATA_H_
#define _FIT_SIMPLE_READ_DATA_H_

#include<fit_simple/data_info.h>

struct Parser{
  virtual void setup(rawDataCorrelationFunctionD &into, const int nsample, const int Lt) const = 0;
  virtual void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const int Lt) const = 0;  
};
struct ParseStandard: public Parser{ //expect Lt lines with format <t> <re> <im>.  Discards the imaginary part
  void setup(rawDataCorrelationFunctionD &into, const int nsample, const int Lt) const{
    into.resize(Lt);
    for(int i=0;i<Lt;i++) into.value(i).resize(nsample);    
  }    
  
  void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const int Lt) const{
    for(int t=0;t<Lt;t++){
      int tt; double im; double &re = into.value(t).sample(sample);
      is >> tt >> re >> im;
      assert(tt == t);
      assert(!is.fail());
      into.coord(t) = t;
    }
  }
};
struct ParseMultiSourceAverage: public Parser{ //expect Lt lines with format <tsrc> <tsep> <re> <im>.  Discards the imaginary part
  void setup(rawDataCorrelationFunctionD &into, const int nsample, const int Lt) const{
    into.resize(Lt);
    for(int t=0;t<Lt;t++){
      into.coord(t) = t;
      into.value(t).resize(nsample);
    }
  }    
  
  void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const int Lt) const{
    for(int t=0;t<Lt;t++)
      into.value(t).sample(sample) = 0.;
    
    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int t=0;t<Lt;t++){
	int tsrct, tt; double re, im;
	is >> tsrct >> tt >> re >> im;
	assert(tsrct == tsrc);
	assert(tt == t);
	assert(!is.fail());

	into.value(t).sample(sample) += re / double(Lt);	
      }
  } 
};

struct ParseMultiSourceAverageImag: public Parser{ //expect Lt lines with format <tsrc> <tsep> <re> <im>.  Discards the real part
  void setup(rawDataCorrelationFunctionD &into, const int nsample, const int Lt) const{
    into.resize(Lt);
    for(int t=0;t<Lt;t++){
      into.coord(t) = t;
      into.value(t).resize(nsample);
    }
  }    
  
  void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const int Lt) const{
    for(int t=0;t<Lt;t++)
      into.value(t).sample(sample) = 0.;
    
    for(int tsrc=0;tsrc<Lt;tsrc++)
      for(int t=0;t<Lt;t++){
	int tsrct, tt; double re, im;
	is >> tsrct >> tt >> re >> im;
	assert(tsrct == tsrc);
	assert(tt == t);
	assert(!is.fail());

	into.value(t).sample(sample) += im / double(Lt);	
      }
  } 
};

Parser* parserFactory(const ParserType p){
  switch(p){
  case ParserType::ParserStandard:
    return new ParseStandard;
  case ParserType::ParserMultiSourceAverage:
    return new ParseMultiSourceAverage;
  case ParserType::ParserMultiSourceAverageImag:
    return new ParseMultiSourceAverageImag;
  default:
    error_exit(std::cout << "readData: Unknown parser " << p << std::endl);
  };
}

void applyOperation(rawDataCorrelationFunctionD &to, const std::string &operation){
  if(operation == "") return;

  expressionAST ast = mathExpressionParse(operation);
  assert(ast.containsSymbol("x"));
  for(int i=0;i<to.size();i++)
    for(int s=0;s<to.value(i).size();s++){
      ast.assignSymbol("x",to.value(i).sample(s));
      to.value(i).sample(s) = ast.value();
    }  
}

template<typename CorrelationFunctionType>
void applyTimeDep(CorrelationFunctionType &to, const TimeDependence tdep, const int Lt){
  if(tdep == TimeDependence::TimeDepNormal) return;

  CorrelationFunctionType cp(to);
  if(tdep == TimeDependence::TimeDepReflect){
    for(int t=0;t<Lt;t++){
      int trefl = (Lt - t) % Lt;
      to.value(trefl) = cp.value(t);
    }
  }else if(tdep == TimeDependence::TimeDepFold){
    for(int t=1;t<Lt;t++)
      to.value(t) = ( cp.value(t) + cp.value(Lt-t) )/2.;
  }else if(tdep == TimeDependence::TimeDepAntiFold){
    for(int t=1;t<Lt;t++)
      to.value(t) = ( cp.value(t) - cp.value(Lt-t) )/2.;
  }else{
    error_exit(std::cout << "applyTimeDep unknown time dependence " << tdep << std::endl);
  }  
}



void readData(rawDataCorrelationFunctionD &into, const DataInfo &data_info, const int Lt, const int traj_start, const int traj_inc, const int traj_lessthan){
  std::size_t off = data_info.file_fmt.find("%d");
  if(off == std::string::npos) error_exit(std::cout << "readData expect file_fmt to contain a '%d', instead got " << data_info.file_fmt << std::endl);
  Parser* parser = parserFactory(data_info.parser);

  int nsample = (traj_lessthan - traj_start)/traj_inc;
  parser->setup(into, nsample, Lt);
#pragma omp parallel for
  for(int s=0;s<nsample;s++){
    std::ostringstream os; os << traj_start + s*traj_inc;
    std::string filename = data_info.file_fmt;
    filename.replace(off,2,os.str());
    std::cout << "Parsing " << filename << std::endl;
    std::ifstream is(filename.c_str());
    if(!is.good()) error_exit(std::cout << "readData failed to read file " << filename << std::endl);
    parser->parse(into, is, s, Lt);
  }  
  delete parser;
  
  applyOperation(into, data_info.operation);
  applyTimeDep(into, data_info.time_dep,Lt);
}

template<typename CorrelationFunctionType>
void applyCombination(CorrelationFunctionType &to, const std::vector<CorrelationFunctionType> &from, const Combination comb){
  if(comb == Combination::CombinationAverage){
    to = from[0];
    for(int i=1;i<from.size();i++) to = to + from[i];
    to = to / double(from.size());
  }else if(comb == Combination::CombinationAminusB){
    if(from.size()!=2) error_exit(std::cout << "applyCombination error: CombinationAminusB requires two channels\n");
    to = from[0] - from[1];
  }else if(comb == Combination::CombinationAdivB){
    if(from.size()!=2) error_exit(std::cout << "applyCombination error: CombinationAdivB requires two channels\n");
    to = from[0]/from[1];
  }else{
    error_exit(std::cout << "applyCombination unknown combination " << comb << std::endl);
  }
}

inline void bin(rawDataCorrelationFunctionD &data, const int bin_size){
  if(bin_size == 1) return;
  for(int d=0;d<data.size();d++) data.value(d) = data.value(d).bin(bin_size);
}

template<typename DistributionType>
inline DistributionType removeSamplesInRange(const DistributionType &raw, const int start, const int lessthan){
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

#endif

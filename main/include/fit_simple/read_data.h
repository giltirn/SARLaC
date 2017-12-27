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
    for(int i=0;i<Lt;i++) into.value(i).resize(nsample);    
  }    
  
  void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const int Lt) const{
    for(int t=0;t<Lt;t++) into.value(t).sample(sample) = 0.;

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

Parser* parserFactory(const ParserType p){
  switch(p){
  case ParserStandard:
    return new ParseStandard;
  case ParserMultiSourceAverage:
    return new ParseMultiSourceAverage;
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

void applyTimeDep(rawDataCorrelationFunctionD &to, const TimeDependence tdep, const int Lt){
  if(tdep == TimeDepNormal) return;

  rawDataCorrelationFunctionD cp(to);
  if(tdep == TimeDepReflect){
    for(int t=0;t<Lt;t++){
      int trefl = (Lt - t) % Lt;
      to.value(trefl) = cp.value(t);
    }
  }else if(tdep == TimeDepFold){
    for(int t=1;t<Lt;t++)
      to.value(t) = ( cp.value(t) + cp.value(Lt-t) )/2.;
  }else if(tdep == TimeDepAntiFold){
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

void applyCombination(rawDataCorrelationFunctionD &to, const std::vector<rawDataCorrelationFunctionD> &from, const Combination comb){
  if(comb == CombinationAverage){
    to = from[0];
    for(int i=1;i<from.size();i++) to = to + from[i];
    to = to / double(from.size());
  }else{
    error_exit(std::cout << "applyCombination unknown combination " << comb << std::endl);
  }
}


#endif

#ifndef _FIT_SIMPLE_READ_DATA_H_
#define _FIT_SIMPLE_READ_DATA_H_



struct Parser{
  virtual void setup(rawDataCorrelationFunctionD &into, const int nsample, const Args &args) const = 0;
  virtual void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const Args &args) const = 0;  
};
struct ParseStandard: public Parser{ //expect Lt lines with format <t> <re> <im>.  Discards the imaginary part
  void setup(rawDataCorrelationFunctionD &into, const int nsample, const Args &args) const{
    into.resize(args.Lt);
    for(int i=0;i<args.Lt;i++) into.value(i).resize(nsample);    
  }    
  
  void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const Args &args) const{
    for(int t=0;t<args.Lt;t++){
      int tt; double im; double &re = into.value(t).sample(sample);
      is >> tt >> re >> im;
      assert(tt == t);
      assert(!is.fail());
      into.coord(t) = t;
    }
  }
};
struct ParseMultiSourceAverage: public Parser{ //expect Lt lines with format <tsrc> <tsep> <re> <im>.  Discards the imaginary part
  void setup(rawDataCorrelationFunctionD &into, const int nsample, const Args &args) const{
    into.resize(args.Lt);
    for(int i=0;i<args.Lt;i++) into.value(i).resize(nsample);    
  }    
  
  void parse(rawDataCorrelationFunctionD &into, std::istream &is, const int sample, const Args &args) const{
    for(int t=0;t<args.Lt;t++) into.value(t).sample(sample) = 0.;

    for(int tsrc=0;tsrc<args.Lt;tsrc++)
      for(int t=0;t<args.Lt;t++){
	int tsrct, tt; double re, im;
	is >> tsrct >> tt >> re >> im;
	assert(tsrct == tsrc);
	assert(tt == t);
	assert(!is.fail());

	into.value(t).sample(sample) += re / double(args.Lt);	
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

void applyTimeDep(rawDataCorrelationFunctionD &to, const TimeDependence tdep, const Args &args){
  if(tdep == TimeDepNormal) return;

  rawDataCorrelationFunctionD cp(to);
  if(tdep == TimeDepReflect){
    const int Lt = args.Lt;
    for(int t=0;t<Lt;t++){
      int trefl = (Lt - t) % Lt;
      to.value(trefl) = cp.value(t);
    }
  }else if(tdep == TimeDepFold){
    const int Lt = args.Lt;
    for(int t=1;t<Lt;t++)
      to.value(t) = ( cp.value(t) + cp.value(Lt-t) )/2.;
  }else if(tdep == TimeDepAntiFold){
    const int Lt = args.Lt;
    for(int t=1;t<Lt;t++)
      to.value(t) = ( cp.value(t) - cp.value(Lt-t) )/2.;
  }else{
    error_exit(std::cout << "applyTimeDep unknown time dependence " << tdep << std::endl);
  }  
}



void readData(rawDataCorrelationFunctionD &into, const DataInfo &data_info, const Args &args){
  std::size_t off = data_info.file_fmt.find("%d");
  if(off == std::string::npos) error_exit(std::cout << "readData expect file_fmt to contain a '%d', instead got " << data_info.file_fmt << std::endl);
  Parser* parser = parserFactory(data_info.parser);

  int nsample = (args.traj_lessthan - args.traj_start)/args.traj_inc;
  parser->setup(into, nsample, args);
  for(int s=0;s<nsample;s++){
    std::ostringstream os; os << s;
    std::string filename = data_info.file_fmt;
    filename.replace(off,2,os.str());
    std::ifstream is(filename.c_str());
    assert(is.good());
    parser->parse(into, is, s, args);
  }  
  delete parser;
  
  applyOperation(into, data_info.operation);
  applyTimeDep(into, data_info.time_dep, args);
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

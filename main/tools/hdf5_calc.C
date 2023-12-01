#include<config.h>

#ifdef HAVE_HDF5
#include<iostream>
#include<cassert>
#include<sstream>
#include<array>

#include<distribution.h>
#include<parser.h>

using namespace CPSfit;

//A program to perform calculations using distributions loaded from canonical format hdf5 containers

struct symbolInfo{
  std::string symbol;
  std::string filename;
  std::vector<int> indices;

  //Info for reading the hdf5 containers
  DistributionTypeEnum type;
  int vector_depth;
  
  void parseIndices(const std::string &ind){
    indices.resize(0);    
    std::istringstream buffer(ind);
    std::string line;
    while(std::getline(buffer, line, ',')){
      int i;
      std::stringstream ss; ss << line; ss >> i;
      indices.push_back(i);
    }
  }
    
  symbolInfo(){}
  symbolInfo(const std::string &s, const std::string &f, const std::string &i): symbol(s), filename(f){
    parseIndices(i);
    getTypeInfo(type,vector_depth,filename);
    if(indices.size() != vector_depth) error_exit(std::cout << "For symbol " << s << " of filename " << f << " expected an index string of size " << vector_depth << std::endl);
  }
};


template<typename D>
void readSymbol(D &into, const symbolInfo &info){
  if(info.vector_depth == 1){
    std::vector<D> tmp;
    readParamsStandard(tmp, info.filename);
    into = tmp[info.indices[0]];
  }else if(info.vector_depth == 2){
    std::vector<std::vector<D> > tmp;
    readParamsStandard(tmp, info.filename);
    into = tmp[info.indices[0]][info.indices[1]];
  }else{
    error_exit(std::cout << "Do not currently support vector depths >2, got " << info.vector_depth << std::endl);
  }
}

struct runOptions{
  bool unite_superjack_layout;
  bool rename_nonmatching_superjack_ensembles; //if the symbols have ensemble named the same thing but with different sizes, perform a renaming such that they are treated as independent

  runOptions(){
    unite_superjack_layout = false;
    rename_nonmatching_superjack_ensembles = false;
  }
};


template<typename D>
void check(const std::vector<symbolInfo> &symbols, std::vector<D> &symbol_values, const runOptions &opt){
  for(int i=1;i<symbols.size();i++)
    if(symbol_values[i].size() != symbol_values[0].size())
      error_exit(std::cout << "All distributions must have the same size\n");
}

template<>
void check<superJackknifeDistribution<double> >(const std::vector<symbolInfo> &symbols, std::vector<superJackknifeDistribution<double> > &symbol_values, const runOptions &opt){
  if(opt.unite_superjack_layout){
    superJackknifeLayout* layout = new superJackknifeLayout;
    for(int i=0;i<symbols.size();i++)
      *layout = combine(*layout, symbol_values[i].getLayout());
    for(int i=0;i<symbols.size();i++)
      symbol_values[i].setLayout(*layout);
  }else{
    for(int i=1;i<symbols.size();i++)
      if(symbol_values[i].getLayout() != symbol_values[0].getLayout())
	error_exit(std::cout << "All superjackknife distributions must have the same layout\n");
  }
}

template<>
void check<superMultiDistribution<double> >(const std::vector<symbolInfo> &symbols, std::vector<superMultiDistribution<double> > &symbol_values, const runOptions &opt){
  if(opt.rename_nonmatching_superjack_ensembles){
    std::vector<superMultiLayout*> layouts(symbols.size());
    for(int i=0;i<symbols.size();i++){
      layouts[i] = new superMultiLayout(symbol_values[i].getLayout());
      symbol_values[i].setLayout(*layouts[i]);
    }

    std::map<std::string, std::set<int> > szmap;
    for(int i=0;i<symbols.size();i++){
      for(int e=0;e<layouts[i]->nEnsembles();e++){
	const std::string &nm = layouts[i]->ensTag(e);
	const int sz = layouts[i]->nSamplesEns(e);
	szmap[nm].insert(sz);
      }
    }
    for(auto mip=szmap.begin();mip!=szmap.end();mip++){
      if(mip->second.size() > 1){
	std::string tag = mip->first;
	for(int i=0;i<symbols.size();i++){
	  int ens = layouts[i]->ensIdx(tag);
	  if(ens != -1){
	    std::string newnm = stringize("%s_%d", tag.c_str(), layouts[i]->nSamplesEns(ens));
	    std::cout << "Renaming " << tag << " to " << newnm << " for symbol " << symbols[i].symbol << std::endl;
	    layouts[i]->setEnsTag(ens, newnm);
	  }
	}
      }
    }
  }

  if(opt.unite_superjack_layout){
    superMultiLayout* layout = new superMultiLayout;
    for(int i=0;i<symbols.size();i++)
      *layout = combine(*layout, symbol_values[i].getLayout());
    for(int i=0;i<symbols.size();i++)
      symbol_values[i].setLayout(*layout);
  }else{
    for(int i=1;i<symbols.size();i++)
      if(symbol_values[i].getLayout() != symbol_values[0].getLayout())
	error_exit(std::cout << "All superjackknife distributions must have the same layout\n");
  }
}


template<typename D>
void specDtype(const std::vector<symbolInfo> &symbols, const DistributionTypeEnum type, expressionAST &expr, const std::string &outfile, const runOptions &opt){
  std::vector<D> symbol_values(symbols.size());
  for(int i=0;i<symbols.size();i++)
    readSymbol(symbol_values[i],symbols[i]);

  check<D>(symbols,symbol_values, opt);

  for(int i=0;i<symbols.size();i++)
    std::cout << "Read symbol " << symbols[i].symbol << " with value " << symbol_values[i] << std::endl;
  
  D out(symbol_values[0]);
  zeroit(out);

  typedef iterate<D> it;
  for(int i=0;i<it::size(out);i++){
    for(int s=0;s<symbols.size();s++) expr.assignSymbol(symbols[s].symbol, it::at(i,symbol_values[s]));
    it::at(i, out) = expr.value();
  }

  std::cout << "Result: " << out << std::endl;

  writeParamsStandard(out,outfile);
}


void run(const std::vector<symbolInfo> &symbols, const DistributionTypeEnum type, expressionAST &expr, const std::string &outfile, const runOptions &opt){
  switch(type){
  case DistributionTypeEnum::Raw:
    specDtype<rawDataDistribution<double> >(symbols,type,expr,outfile,opt);  break;
  case DistributionTypeEnum::Jackknife:
    specDtype<jackknifeDistribution<double> >(symbols,type,expr,outfile,opt);  break;
  case DistributionTypeEnum::JackknifeC:
    specDtype<jackknifeCdistribution<double> >(symbols,type,expr,outfile,opt);  break;
  case DistributionTypeEnum::DoubleJackknife:
    specDtype<doubleJackknifeDistribution<double> >(symbols,type,expr,outfile,opt);  break;
  case DistributionTypeEnum::SuperJackknife:
    specDtype<superJackknifeDistribution<double> >(symbols,type,expr,outfile,opt);  break;    
  case DistributionTypeEnum::SuperMulti:
    specDtype<superMultiDistribution<double> >(symbols,type,expr,outfile,opt);  break;    
  case DistributionTypeEnum::Bootstrap:
    specDtype<bootstrapDistribution<double> >(symbols,type,expr,outfile,opt);  break;
  default:
    error_exit(std::cout << "run(...) unknown type " << type << std::endl);
  }
}
  

int main(const int argc, const char* argv[]){
  if(argc < 3) error_exit(std::cout << "Usage: <options> <out file> <expression> <symbol 1 name> <symbol 1 filename> <symbol 1 index expr> <symbol 2 name> ...." << std::endl);

  int argc_prune;
  std::vector<std::string> argv_prune;
  
  runOptions opt;
  for(int i=0;i<argc;i++){
    std::string arg(argv[i]);
    if(arg == "-unite_superjack_layout"){
      opt.unite_superjack_layout = true;
    }else if(arg == "-rename_nonmatching_superjack_ensembles"){
      opt.rename_nonmatching_superjack_ensembles = true;
    }else{
      argv_prune.push_back(arg);
    }
  }
  argc_prune = argv_prune.size();

  //First argument should be the output file name
  std::string outfile = argv_prune[1];
  
  //Second argument should be an expression
  expressionAST expr = mathExpressionParse(argv_prune[2]);

  //Remaining arguments should be sets of 3 strings: the symbol name followed by the filename and finally the index expression
  std::vector<symbolInfo> symbols;
  if( (argc_prune-3) % 3 != 0 ) error_exit(std::cout << "Symbol arguments must come in sets of 3\n");
  for(int i=0;i< (argc_prune-3)/3; i++){
    int base = i*3+3;
    symbols.push_back(symbolInfo(argv_prune[base],argv_prune[base+1],argv_prune[base+2]));
  }
  if(symbols.size() == 0) error_exit(std::cout << "Should have at least one symbol!\n");

  //All symbols must have the same distribution type
  for(int i=1;i<symbols.size();i++)
    if(symbols[i].type != symbols[0].type) error_exit(std::cout << "Symbol distribution types must be the same\n");
    
  run(symbols, symbols[0].type, expr, outfile, opt);
  return 0;
}


#else

int main(const int argc, const char* argv[]){
  std::cout << "Require HDF5\n";
  return 1;
}

#endif

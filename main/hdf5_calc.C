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

template<typename D>
void check(const std::vector<symbolInfo> &symbols, const std::vector<D> &symbol_values){
  for(int i=1;i<symbols.size();i++)
    if(symbol_values[i].size() != symbol_values[0].size())
      error_exit(std::cout << "All superjackknife distributions must have the same size\n");
}

template<>
void check<superJackknifeDistribution<double> >(const std::vector<symbolInfo> &symbols, const std::vector<superJackknifeDistribution<double> > &symbol_values){
  for(int i=1;i<symbols.size();i++)
    if(symbol_values[i].getLayout() != symbol_values[0].getLayout())
      error_exit(std::cout << "All superjackknife distributions must have the same layout\n");
}




template<typename D>
void specDtype(const std::vector<symbolInfo> &symbols, const DistributionTypeEnum type, expressionAST &expr, const std::string &outfile){
  std::vector<D> symbol_values(symbols.size());
  for(int i=0;i<symbols.size();i++)
    readSymbol(symbol_values[i],symbols[i]);

  check<D>(symbols,symbol_values);

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


void run(const std::vector<symbolInfo> &symbols, const DistributionTypeEnum type, expressionAST &expr, const std::string &outfile){
  switch(type){
  case DistributionTypeEnum::Raw:
    specDtype<rawDataDistribution<double> >(symbols,type,expr,outfile);  break;
  case DistributionTypeEnum::Jackknife:
    specDtype<jackknifeDistribution<double> >(symbols,type,expr,outfile);  break;
  case DistributionTypeEnum::JackknifeC:
    specDtype<jackknifeCdistribution<double> >(symbols,type,expr,outfile);  break;
  case DistributionTypeEnum::DoubleJackknife:
    specDtype<doubleJackknifeDistribution<double> >(symbols,type,expr,outfile);  break;
  case DistributionTypeEnum::SuperJackknife:
    specDtype<superJackknifeDistribution<double> >(symbols,type,expr,outfile);  break;    
  case DistributionTypeEnum::SuperMulti:
    specDtype<superMultiDistribution<double> >(symbols,type,expr,outfile);  break;    
  case DistributionTypeEnum::Bootstrap:
    specDtype<bootstrapDistribution<double> >(symbols,type,expr,outfile);  break;
  default:
    error_exit(std::cout << "run(...) unknown type " << type << std::endl);
  }
}
  

int main(const int argc, const char* argv[]){
  if(argc < 3) error_exit(std::cout << "Usage: <out file> <expression> <symbol 1 name> <symbol 1 filename> <symbol 1 index expr> <symbol 2 name> ...." << std::endl);

  //First argument should be the output file name
  std::string outfile = argv[1];
  
  //Second argument should be an expression
  expressionAST expr = mathExpressionParse(argv[2]);

  //Remaining arguments should be sets of 3 strings: the symbol name followed by the filename and finally the index expression
  std::vector<symbolInfo> symbols;
  if( (argc-3) % 3 != 0 ) error_exit(std::cout << "Symbol arguments must come in sets of 3\n");
  for(int i=0;i< (argc-3)/3; i++){
    int base = i*3+3;
    symbols.push_back(symbolInfo(argv[base],argv[base+1],argv[base+2]));
  }
  if(symbols.size() == 0) error_exit(std::cout << "Should have at least one symbol!\n");

  //All symbols must have the same distribution type
  for(int i=1;i<symbols.size();i++)
    if(symbols[i].type != symbols[0].type) error_exit(std::cout << "Symbol distribution types must be the same\n");
    
  run(symbols, symbols[0].type, expr, outfile);
  return 0;
}


#else

int main(const int argc, const char* argv[]){
  std::cout << "Require HDF5\n";
  return 1;
}

#endif

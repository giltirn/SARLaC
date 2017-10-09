#ifndef _FIT_WRAPPER_FREEZE_H_
#define _FIT_WRAPPER_FREEZE_H_

#include<parser.h>
#include<expression_parse.h>
#include<fit_wrapper.h>
#include<distribution_IO.h>

//Convenience functions and types for loading data for frozen fits

GENERATE_ENUM_AND_PARSER(FreezeDataReaderType, (UKfitXMLvectorReader)(HDF5fileReader) );

void readUKfitVectorEntry(jackknifeDistribution<double> &into, const std::string &filename, const int idx){
  XMLreader rd(filename);
  std::cout << "readUKfitVectorEntry reading file " << filename << std::endl;
  UKvalenceDistributionContainer<jackknifeDistribution<double> > con;
  read(rd,con,"data_in_file");
  std::cout << "Read " << con.Nentries << " distributions from file " << filename << std::endl;  
  into = con.list[idx];
}
#ifdef HAVE_HDF5
void readHDF5file(jackknifeDistribution<double> &into, const std::string &filename, const std::vector<int> &idx){
  DistributionTypeEnum type;
  int depth;
  getTypeInfo(type,depth,filename);
  if(type != Jackknife) error_exit(std::cout << "readHDF5file expected a jackknife distribution\n");

  if(depth == 1){
    assert(idx.size() == 1);
    std::vector<jackknifeDistribution<double> > tmp;
    readParamsStandard(tmp,filename);
    into = tmp[idx[0]];
  }else if(depth == 2){
    assert(idx.size() == 2);
    std::vector<std::vector<jackknifeDistribution<double> > > tmp;
    readParamsStandard(tmp,filename);
    into = tmp[idx[0]][idx[1]];
  }
}
#endif

//input_idx is the data file index: If the data file is a vector<vector> type you must provide 2 indices
#define FREEZE_PARAM_MEMBERS \
  (int, param_idx) \
  (FreezeDataReaderType, reader) \
  (std::string, filename) \
  (std::vector<int>, input_idx)	\
  (std::string, operation)

struct FreezeParam{
  GENERATE_MEMBERS(FREEZE_PARAM_MEMBERS);

  FreezeParam(): param_idx(0), reader(UKfitXMLvectorReader), filename("file.dat"), input_idx(1,0), operation(""){}
};
GENERATE_PARSER(FreezeParam, FREEZE_PARAM_MEMBERS);


#define FREEZE_PARAMS_MEMBERS			\
  (std::vector<FreezeParam>, sources)

struct FreezeParams{
  GENERATE_MEMBERS(FREEZE_PARAMS_MEMBERS);
  
  FreezeParams(): sources(1){}
};
GENERATE_PARSER(FreezeParams, FREEZE_PARAMS_MEMBERS);

#undef FREEZE_PARAMS_MEMBERS

template<typename FitFuncPolicies>
void readFrozenParams(fitter<FitFuncPolicies> &fitter, const std::string &freeze_file, const int nsample){
  typedef typename FitFuncPolicies::jackknifeFitParameters jackknifeFitParameters;
  jackknifeFitParameters values(nsample);
  
  std::vector<int> freeze;
  FreezeParams fparams;
  parse(fparams,freeze_file);
  
  for(int i=0;i<fparams.sources.size();i++){
    std::cout << "readFrozenParams loading freeze data for parameter " << fparams.sources[i].param_idx << std::endl;
    freeze.push_back(fparams.sources[i].param_idx);

    jackknifeDistribution<double> fval;
    if(fparams.sources[i].reader == UKfitXMLvectorReader) readUKfitVectorEntry(fval, fparams.sources[i].filename, fparams.sources[i].input_idx[0]);
    else if(fparams.sources[i].reader == HDF5fileReader){
#ifdef HAVE_HDF5
      readHDF5file(fval, fparams.sources[i].filename, fparams.sources[i].input_idx);
#else
      error_exit(std::cout << "readFrozenParams with reader " << HDF5fileReader << " : must have compiled with HDF5 enabled\n");
#endif    
    }else error_exit(std::cout << "readFrozenParams unknown reader " << fparams.sources[i].reader << std::endl);

    if(fval.size() != nsample) error_exit(std::cout << "readFrozenParams read jackknife of size " << fval.size() << ", expected " << nsample << std::endl);

    if(fparams.sources[i].operation != ""){
      expressionAST AST = mathExpressionParse(fparams.sources[i].operation);
      if(!AST.containsSymbol("x")) error_exit(std::cout << "readFrozenParams expect math expression to be a function of 'x', got \"" << fparams.sources[i].operation << "\"\n");
      for(int s=0;s<nsample;s++){
	AST.assignSymbol("x",fval.sample(s));
	fval.sample(s) = AST.value();
      }
    }

    std::cout << "readFrozenParams read " << fval << std::endl;
    
    for(int s=0;s<nsample;s++)
      values.sample(s)(fparams.sources[i].param_idx) = fval.sample(s);    
  }

  fitter.freeze(freeze, values);
}

#endif

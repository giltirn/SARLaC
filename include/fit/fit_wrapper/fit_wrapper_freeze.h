#ifndef _FIT_WRAPPER_FREEZE_H_
#define _FIT_WRAPPER_FREEZE_H_

//Convenience functions and types for loading data for frozen fits

#include<config.h>
#include<utils/macros.h>
#include<utils/utils.h>
#include<parser/parser.h>
#include<parser/expression_parse.h>
#include<fit/fit_wrapper/fitter.h>
#include<distribution/jackknife.h>
#include<distribution/bootstrap.h>
#include<distribution/distribution_hdf5io_conventional.h>

CPSFIT_START_NAMESPACE

GENERATE_ENUM_AND_PARSER(FreezeDataReaderType, (UKfitXMLvectorReader)(HDF5fileReader)(ConstantValue) );

//Read a jackknife from a UKfit-style bootxml
void readUKfitVectorEntry(jackknifeDistribution<double> &into, const std::string &filename, const int idx){
  XMLreader rd(filename);
  std::cout << "readUKfitVectorEntry reading file " << filename << std::endl;
  UKvalenceDistributionContainer<jackknifeDistribution<double> > con;
  read(rd,con,"data_in_file");
  std::cout << "Read " << con.Nentries << " distributions from file " << filename << std::endl;  
  into = con.list[idx];
}
void readUKfitVectorEntry(bootstrapDistribution<double> &into, const std::string &filename, const int idx){
  error_exit(std::cout << "readUKfitVectorEntry not implemented for bootstrap\n");
}

//Read a jackknife from a HDF5 file
template<typename DistributionType>
void readHDF5file(DistributionType &into, const std::string &filename, const std::vector<int> &idx){
#ifdef HAVE_HDF5
  DistributionTypeEnum type;
  int depth;
  getTypeInfo(type,depth,filename);

  DistributionTypeEnum expect = getTypeEnum(getDistributionTypeString<DistributionType>());

  if(type != expect) error_exit(std::cout << "readHDF5file expected a " << expect << " distribution\n");

  if(depth == 1){
    assert(idx.size() == 1);
    std::vector<DistributionType> tmp;
    readParamsStandard(tmp,filename);
    into = tmp[idx[0]];
  }else if(depth == 2){
    assert(idx.size() == 2);
    std::vector<std::vector<DistributionType> > tmp;
    readParamsStandard(tmp,filename);
    into = tmp[idx[0]][idx[1]];
  }
#else
  error_exit(std::cout << "readFrozenParams with reader " << HDF5fileReader << " : must have compiled with HDF5 enabled\n");
#endif
}


//input_idx is the data file index: If the data file is a vector<vector> type you must provide 2 indices
struct FreezeParam{
#define FREEZE_PARAM_MEMBERS \
  (int, param_idx) \
  (FreezeDataReaderType, reader) \
  (std::string, filename) \
  (std::vector<int>, input_idx)	\
  (std::string, operation)

  GENERATE_MEMBERS(FREEZE_PARAM_MEMBERS);

  FreezeParam(): param_idx(0), reader(FreezeDataReaderType::UKfitXMLvectorReader), filename("file.dat"), input_idx(1,0), operation(""){}
};
GENERATE_PARSER(FreezeParam, FREEZE_PARAM_MEMBERS);


struct FreezeParams{
#define FREEZE_PARAMS_MEMBERS			\
  (std::vector<FreezeParam>, sources)

  GENERATE_MEMBERS(FREEZE_PARAMS_MEMBERS);
  
  FreezeParams(): sources(1){}
};
GENERATE_PARSER(FreezeParams, FREEZE_PARAMS_MEMBERS);

#undef FREEZE_PARAMS_MEMBERS


//Apply a math expression to a jackknife of input data
template<typename DistributionType>
void applyOperation(DistributionType &fval, const std::string &operation, const FreezeDataReaderType reader){
  if(operation == ""){
    if(reader == FreezeDataReaderType::ConstantValue) error_exit(std::cout << "readFrozenParams with ConstantValue, require math expression, got \"" << operation << "\"\n");
    return;
  }
  typedef iterate<DistributionType> iter;

  const int nsample = iter::size(fval);
  expressionAST AST = mathExpressionParse(operation);

  if(reader == FreezeDataReaderType::ConstantValue){
    if(AST.nSymbols() != 0) error_exit(std::cout << "readFrozenParams with ConstantValue expects math expression with no symbols, got \"" << operation << "\"\n");
    for(int s=0;s<nsample;s++)
      iter::at(s, fval) = AST.value();
  }else{
    if(AST.nSymbols() != 1) error_exit(std::cout << "readFrozenParams with " << reader << " expects math expression with 1 symbol ('x'), got \"" << operation << "\"\n");
    else if(!AST.containsSymbol("x")) error_exit(std::cout << "readFrozenParams with " << reader << " expects math expression to be a function of 'x', got \"" << operation << "\"\n");

    for(int s=0;s<nsample;s++){
      AST.assignSymbol("x", iter::at(s, fval));
      iter::at(s, fval) = AST.value();
    }    
  }
}

//The main function - read and import the  frozen parameters. A struct "FreezeParams" is read in from "freeze_file" and used for perform the required actions
//For parameter types that don't have a default constructor the user should provide a pointer 'psetup' to a setup instance of the parameter type
template<typename FitFuncPolicies>
void readFrozenParams(fitter<FitFuncPolicies> &fitter, const std::string &freeze_file, const int nsample, typename FitFuncPolicies::baseFitFunc::ParameterType const* psetup = NULL){
  if(!fileExists(freeze_file)){
    FreezeParams templ;
    std::ofstream of("freeze_template.dat");
    of << templ;
    of.close();
    error_exit(std::cout << "Failed to read freeze file " << freeze_file << "; wrote template to freeze_template.dat\n");
  } 

  typedef typename FitFuncPolicies::FitParameterDistribution FitParameterDistribution;
  FitParameterDistribution values(nsample);
  if(psetup!=NULL) for(int s=0;s<nsample;s++) values.sample(s) = *psetup;
  
  std::vector<int> freeze;
  FreezeParams fparams;
  parse(fparams,freeze_file);
  
  for(int i=0;i<fparams.sources.size();i++){
    std::cout << "readFrozenParams loading freeze data for parameter " << fparams.sources[i].param_idx << std::endl;
    freeze.push_back(fparams.sources[i].param_idx);

    jackknifeDistribution<double> fval;

    //Different sources of data
    FreezeDataReaderType reader = fparams.sources[i].reader;
    
    if(reader == FreezeDataReaderType::UKfitXMLvectorReader){
      readUKfitVectorEntry(fval, fparams.sources[i].filename, fparams.sources[i].input_idx[0]);
    }else if(reader == FreezeDataReaderType::HDF5fileReader){
      readHDF5file(fval, fparams.sources[i].filename, fparams.sources[i].input_idx);
    }else if(reader == FreezeDataReaderType::ConstantValue){
      fval.resize(nsample);
    }else{
      error_exit(std::cout << "readFrozenParams unknown reader " << fparams.sources[i].reader << std::endl);
    }

    if(fval.size() != nsample) error_exit(std::cout << "readFrozenParams read jackknife of size " << fval.size() << ", expected " << nsample << std::endl);

    applyOperation(fval, fparams.sources[i].operation, reader);

    std::cout << "readFrozenParams read " << fval << std::endl;
    
    for(int s=0;s<nsample;s++)
      values.sample(s)(fparams.sources[i].param_idx) = fval.sample(s);    
  }

  fitter.freeze(freeze, values);
}

CPSFIT_END_NAMESPACE

#endif

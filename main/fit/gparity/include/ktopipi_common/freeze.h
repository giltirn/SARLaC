#ifndef _KTOPIPI_FREEZE_H_
#define _KTOPIPI_FREEZE_H_

#include<config.h>
#include<utils/macros.h>

#include<fit.h>

CPSFIT_START_NAMESPACE

#define KTOPIPI_FREEZE_PARAM_MEMBERS \
  (std::vector<int>, Qlist) \
  (int, param_idx) \
  (FreezeDataReaderType, reader) \
  (std::string, filename) \
  (int, input_idx) \
  (std::string, operation)

struct KtoPiPiFreezeParam{
  //If Qlist is left empty it is assumed the freeze is the same for all Q
  //operation can be any math expression. Use the variable 'x' to represent the data.  eg operation = sqrt(1e13 * x).  Use an empty string for no operation to be performed
  GENERATE_MEMBERS(KTOPIPI_FREEZE_PARAM_MEMBERS);

  KtoPiPiFreezeParam(): param_idx(0), reader(FreezeDataReaderType::UKfitXMLvectorReader), filename("file.dat"), input_idx(0), operation(""), Qlist(0){}
};
GENERATE_PARSER(KtoPiPiFreezeParam, KTOPIPI_FREEZE_PARAM_MEMBERS);


#define KTOPIPI_FREEZE_PARAMS_MEMBERS			\
  (std::vector<KtoPiPiFreezeParam>, sources)

struct KtoPiPiFreezeParams{
  GENERATE_MEMBERS(KTOPIPI_FREEZE_PARAMS_MEMBERS);
  
  KtoPiPiFreezeParams(): sources(1){}
};
GENERATE_PARSER(KtoPiPiFreezeParams, KTOPIPI_FREEZE_PARAMS_MEMBERS);

//Use Q=-1 for all Q
template<typename FitFuncPolicies>
void readFrozenParams(fitter<FitFuncPolicies> &fitter, const int Q, const std::string &freeze_file, const int nsample){
  typedef typename FitFuncPolicies::FitParameterDistribution FitParameterDistribution;
  FitParameterDistribution values(nsample);
  
  std::vector<int> freeze;
  KtoPiPiFreezeParams fparams;
  parse(fparams,freeze_file);
  
  for(int i=0;i<fparams.sources.size();i++){
    if(Q!=-1 && fparams.sources[i].Qlist.size() > 0 && std::find(fparams.sources[i].Qlist.begin(), fparams.sources[i].Qlist.end(), Q) == fparams.sources[i].Qlist.end()) continue;
    
    std::cout << "readFrozenParams loading freeze data for parameter " << fparams.sources[i].param_idx << std::endl;
    freeze.push_back(fparams.sources[i].param_idx);

    jackknifeDistribution<double> fval;
    if(fparams.sources[i].reader == FreezeDataReaderType::UKfitXMLvectorReader) readUKfitVectorEntry(fval, fparams.sources[i].filename, fparams.sources[i].input_idx);
    else if(fparams.sources[i].reader == FreezeDataReaderType::HDF5fileReader) readHDF5file(fval, fparams.sources[i].filename, std::vector<int>(1,fparams.sources[i].input_idx));
    else if(fparams.sources[i].reader == FreezeDataReaderType::ConstantValue) fval.resize(nsample);
    else error_exit(std::cout << "readFrozenParams unknown reader " << fparams.sources[i].reader << std::endl);

    if(fval.size() != nsample) error_exit(std::cout << "readFrozenParams read jackknife of size " << fval.size() << ", expected " << nsample << std::endl);

    if(fparams.sources[i].operation != ""){
      expressionAST AST = mathExpressionParse(fparams.sources[i].operation);
      if(fparams.sources[i].reader == FreezeDataReaderType::ConstantValue){
	if(AST.nSymbols() != 0) error_exit(std::cout << "readFrozenParams with ConstantValue expects math expression with no symbols, got \"" << fparams.sources[i].operation << "\"\n");
	for(int s=0;s<nsample;s++)
	  fval.sample(s) = AST.value();
      }else{
	if(!AST.containsSymbol("x")) error_exit(std::cout << "readFrozenParams expect math expression to be a function of 'x', got \"" << fparams.sources[i].operation << "\"\n");
	for(int s=0;s<nsample;s++){
	  AST.assignSymbol("x",fval.sample(s));
	  fval.sample(s) = AST.value();
	}
      }
    }

    std::cout << "readFrozenParams read " << fval << std::endl;
    
    for(int s=0;s<nsample;s++)
      values.sample(s)(fparams.sources[i].param_idx) = fval.sample(s);    
  }

  fitter.freeze(freeze, values);
}


struct FitKtoPiPiFreezeBase{
  //Q in (1..10)
  virtual std::vector< std::pair<int, jackknifeDistribution<double> > > getFreeze(const int Q, const int nsample) const = 0;
  virtual ~FitKtoPiPiFreezeBase(){}
};

//Same frozen parameters for all Q
struct FitKtoPiPiFreezeFixed: public FitKtoPiPiFreezeBase{
  std::vector< std::pair<int, jackknifeDistribution<double> > > pv;

  std::vector< std::pair<int, jackknifeDistribution<double> > > getFreeze(const int Q, const int nsample) const{ return pv; }
  void addFreeze(const int param, const jackknifeDistribution<double> &val){ pv.push_back( std::pair<int, jackknifeDistribution<double> >(param,val) ); }
};
	      
template<typename FitFuncPolicies>
void readFrozenParams(fitter<FitFuncPolicies> &fitter, const int Q, const FitKtoPiPiFreezeBase &freezer, const int nsample){
  typedef typename FitFuncPolicies::FitParameterDistribution FitParameterDistribution;
  FitParameterDistribution values(nsample);
  
  std::vector< std::pair<int, jackknifeDistribution<double> > > params_vals = freezer.getFreeze(Q, nsample);
  std::vector<int> freeze(params_vals.size());  

  for(int i=0; i<params_vals.size(); i++){
    std::cout << "Freezing param " << params_vals[i].first << " to " << params_vals[i].second << std::endl;
    freeze[i] = params_vals[i].first;
    
    for(int s=0;s<nsample;s++)
      values.sample(s)(freeze[i]) = params_vals[i].second.sample(s);    
  }
  fitter.freeze(freeze, values);
}











CPSFIT_END_NAMESPACE

#endif

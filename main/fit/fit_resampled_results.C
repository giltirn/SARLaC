//Perform a fit to a set of jackknife resampled data

#include<parser.h>
#include<common.h>
#include<fit.h>
#include<plot.h>
#include<containers.h>

using namespace SARLaC;

//Operation is a mathematical expression in the read value 'x' denoting an operation performed on said data post-read. E.g. "x" just reads the data and does nothing (also can use "" for this), "x*x" squares it and so on
#define DATA_INFO_MEMBERS				\
  ( std::string, filename )				\
  ( std::vector<int>, input_idx)			\
  ( std::string, operation )				\
  ( std::vector<double>, coord )

struct DataInfo{
  GENERATE_MEMBERS(DATA_INFO_MEMBERS);

  DataInfo(): input_idx(1,0), coord(1,0), filename("file.hdf5"), operation("x"){}
};
GENERATE_PARSER(DataInfo, DATA_INFO_MEMBERS);

GENERATE_ENUM_AND_PARSER(FitFuncType, (FParabola)(FDispnRelnPow4) );

//(FLinearSquareFourth)(FLatticeDispnReln)(FLatticeDispnRelnDiscE)

#define ARGS_MEMBERS				\
  ( std::vector<DataInfo>, data )		\
  ( FitFuncType, fitfunc)			\
  ( bool, correlated )

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);

  Args(): data(1), fitfunc(FitFuncType::FParabola), correlated(true){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


struct CMDline{
  bool load_frozen_fit_params;
  std::string load_frozen_fit_params_file;

  CMDline(){
    load_frozen_fit_params = false;
  }
  CMDline(const int argc, const char** argv, const int begin = 0): CMDline(){
    setup(argc,argv,begin);
  }
  
  void setup(const int argc, const char** argv, const int begin = 0){
    const int sz = argc-begin;
    std::vector<std::string> sargv(sz);
    for(int i=begin; i<argc; i++) sargv[i-begin] = std::string(argv[i]);

    int i = 0;
    while(i<sz){
      if(sargv[i] == "-load_frozen_fit_params"){
	load_frozen_fit_params = true;
	load_frozen_fit_params_file = sargv[i+1];

	if(load_frozen_fit_params_file == "TEMPLATE"){
	  std::cout << "Saving frozen fit params template file to freeze_template.args" << std::endl;
	  FreezeParams fp;
	  std::ofstream of("freeze_template.args");
	  of << fp;
	  of.close();
	  exit(0);
	}else if(!fileExists(load_frozen_fit_params_file)) error_exit(std::cout << "CMDline freeze data file " << load_frozen_fit_params_file << " does not exist!\n");
	i+=2;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};


typedef correlationFunction<std::vector<double>, jackknifeDistributionD> correlationFunctionType;

//Plot accessor for fit functions that just use entry 0 of the input coordinate vector
struct BasicAccessor{
  const correlationFunctionType &data;
  
  double x(const int i) const{ return data.coord(i)[0]; }
  double dxm(const int i) const{ return 0; }
  double dxp(const int i) const{ return 0; }

  double y(const int i) const{ return data.value(i).best(); }
  double dym(const int i) const{ return data.value(i).standardError(); }
  double dyp(const int i) const{ return data.value(i).standardError(); }
  
  double upper(const int i) const{ return y(i) + dyp(i); }
  double lower(const int i) const{ return y(i) - dym(i); } 

  int size() const{ return data.size(); }
  
  BasicAccessor(const correlationFunctionType &data): data(data){}
};

//Fit  A + B(x-C)^2
class FitParabola{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef std::vector<double> GeneralizedCoordinate;

  typedef BasicAccessor AccessorType;

  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    double dt = t[0]-p(2);
    return p(0) + p(1)*dt*dt;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    double dt = t[0]-p(2);
    ValueDerivativeType yderivs(3);
    yderivs(0) = 1.;
    yderivs(1) = dt*dt;
    yderivs(2) = -2*p(1)*dt;
    return yderivs;
  }

  inline int Nparams() const{ return 3; }

  inline ParameterType basicGuess() const{ ParameterType out(3); out(0) = out(1) = out(2) = 1; return out; }
};
//E^2 = A^2 + \sum_i (B p_i + C p_i^3 + O(p^5) )^2  = A^2 + \sum_i [ B^2*p^2_i + 2*B*C*p_i^4 + O(p^6) ]
//Relabel  B^2 -> B    2*B*C -> C     continuum values B=1 C=0    

//2 sin(x/2) ~ x - 1/24* x^3
//Lattice dispn reln \sum_i [ 2sin(p_i/2) ]^2 ~ \sum_i ( p_i -1/24 *p_i^3 )^2 ~ \sum_i [ p_i^2 -1/12 * p_i^4 ]     B=1 C=-1/12
class FitDispnRelnPow4{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
  typedef std::vector<double> GeneralizedCoordinate;

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    T A = b(p(0),0);
    T B = b(p(1),1);
    T C = b(p(2),2);

    T sum(0.);
    for(int i=0;i<3;i++) sum = sum + B*pow(x[i],2) + C*pow(x[i],4);
    return A*A + sum;
  }
  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d(3);
    for(int i=0;i<3;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  inline int Nparams() const{ return 3; }

  inline ParameterType basicGuess() const{ 
    ParameterType out(3); 
    out(0)=0.1;
    out(1)=1.;
    out(2)=0.;
    return out;
  }

  struct AccessorType{
    const correlationFunctionType &data;
  
    double x(const int i) const{ 
      double sum = 0;
      for(int d=0;d<3;d++) sum = sum + pow(data.coord(i)[d],2); 
      return sum;
    }
    double dxm(const int i) const{ return 0; }
    double dxp(const int i) const{ return 0; }
    
    double y(const int i) const{ return data.value(i).best(); }
    double dym(const int i) const{ return data.value(i).standardError(); }
    double dyp(const int i) const{ return data.value(i).standardError(); }
    
    double upper(const int i) const{ return y(i) + dyp(i); }
    double lower(const int i) const{ return y(i) - dym(i); } 
    
    int size() const{ return data.size(); }

    AccessorType(const correlationFunctionType &data): data(data){}
  };


};

// //Fit A + Bsin(x/2)^2
// class FitLatticeDispnReln{
// public:
//   typedef double ValueType;
//   typedef parameterVector<double> ParameterType;
//   typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
//   typedef double GeneralizedCoordinate;

//   ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
//     return p(0) + p(1)*pow(sin(t/2.),2);
//   }
//   ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
//     ValueDerivativeType yderivs(2);
//     yderivs(0) = 1.;
//     yderivs(1) = pow(sin(t/2.),2);
//     return yderivs;
//   }

//   inline int Nparams() const{ return 2; }

//   inline ParameterType basicGuess() const{ ParameterType out(2); out(0) = out(1) = 1; return out; }
// };

// //[ 2sinh(E/2) ]^2  = [ 2sinh(A/2) ]^2 + Bsin(x/2)^2
// //Data should be for E
// class FitLatticeDispnRelnDiscE{
// public:
//   typedef double ValueType;
//   typedef parameterVector<double> ParameterType;
//   typedef parameterVector<double> ValueDerivativeType; //derivative wrt parameters
//   typedef double GeneralizedCoordinate;

//   template<typename T, typename Boost>
//   T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
//     T A = b(p(0),0);
//     T B = b(p(1),1);
    
//     T v= sqrt( pow(2*sinh(A/2.),2) + B*pow(sin(x/2.),2) );
//     return 2. * asinh(v/2.);
//   }
//   inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
//     return eval<double>(x,p,[&](const double a, const int i){ return a; });
//   }
//   inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
//     ValueDerivativeType d(2);
//     for(int i=0;i<2;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
//     return d;
//   }

//   inline int Nparams() const{ return 2; }

//   inline ParameterType basicGuess() const{ ParameterType out(2); out(0) = out(1) = 1; return out; }
// };



template<typename FitFunc>
void plot(const FitFunc &fitfunc, const correlationFunctionType &data, 
	  const jackknifeDistribution<typename FitFunc::ParameterType> &params){
  const int ndata = data.size();
  const int nsample = data.value(0).size();
  assert(params.size() == nsample);

  correlationFunctionType fit_data(ndata);
  for(int i=0;i<ndata;i++){
    //std::cout << "Computing fit value for coord " << data.coord(i) << std::endl;

    jackknifeDistributionD val(nsample);
    for(int s=0;s<nsample;s++){
      //std::cout << "Sample " << s << " with params " << params.sample(s) << std::endl;
      val.sample(s) = fitfunc.value(data.coord(i),params.sample(s));
    } 
    fit_data.value(i) = val;
    fit_data.coord(i) = data.coord(i);
  }


  typedef MatPlotLibScriptGenerate::handleType Handle;
  typedef typename FitFunc::AccessorType Accessor;

  {
    MatPlotLibScriptGenerate plotter;
    typename MatPlotLibScriptGenerate::kwargsType plot_args;    
    Accessor a(data);
    Handle ah = plotter.plotData(a);
  
    Accessor band(fit_data);
    plot_args["alpha"] = 0.2;
    ah = plotter.errorBand(band, plot_args);
  
    std::cout << "Writing plot to 'plot.py'\n";  
    plotter.write("plot.py", "plot.pdf");
  }

  correlationFunctionType reldiff(ndata);
  for(int i=0;i<ndata;i++){
    reldiff.coord(i) = data.coord(i);
    reldiff.value(i) = (fit_data.value(i) - data.value(i))/data.value(i);
  }
  
  { 
    MatPlotLibScriptGenerate plotter;
    Accessor a(reldiff);
    Handle ah = plotter.plotData(a);
  
    plotter.write("plot_reldiff.py", "plot_reldiff.pdf");
  }

}



template<typename FitFunc>
void fitSpecFF(const correlationFunctionType &data, const FitFunc &fitfunc, const Args &args, const CMDline &cmdline){
  typedef typename composeFitPolicy<FitFunc,frozenFitFuncPolicy,correlatedFitPolicy>::type FitPolicies;
  
  const int nsample = data.value(0).size();
  const int ndata = data.size();
  
  NumericSquareMatrix<double> cov(ndata, [&](const int i, const int j){ return jackknifeDistributionD::covariance(data.value(i),data.value(j)); });
  NumericVector<jackknifeDistributionD> sigma_j(ndata, [&](const int i){ return jackknifeDistributionD(nsample,sqrt(cov(i,i))); });
  
  NumericSquareMatrix<double> corr(ndata, [&](const int i, const int j){ return cov(i,j)/sqrt(cov(i,i)*cov(j,j)); });
  NumericSquareMatrix<double> inv_corr(ndata);
  svd_inverse(inv_corr,corr);

  NumericSquareMatrix<jackknifeDistributionD> inv_corr_j;
  if(args.correlated)
    inv_corr_j = NumericSquareMatrix<jackknifeDistributionD>(ndata, [&](const int i, const int j){ return jackknifeDistributionD(nsample, inv_corr(i,j)); });
  else
    inv_corr_j = NumericSquareMatrix<jackknifeDistributionD>(ndata, [&](const int i, const int j){ return jackknifeDistributionD(nsample, i==j ? 1. : 0.); });

  std::cout << "Correlation matrix:\n" << corr << std::endl;
  std::cout << "Sigma: " << sigma_j << std::endl;

  MarquardtLevenbergParameters<double> mlparams;
  mlparams.delta_cost_min = 1e-8;

  typename FitFunc::ParameterType guess = fitfunc.basicGuess();

  fitter<FitPolicies> fitter(mlparams);
  fitter.importCostFunctionParameters(inv_corr_j,sigma_j);
  fitter.importFitFunc(fitfunc);

  if(cmdline.load_frozen_fit_params)
    readFrozenParams(fitter, cmdline.load_frozen_fit_params_file, nsample, &guess);
  else
    fitter.freeze({}, jackknifeDistribution<typename FitFunc::ParameterType>(nsample, guess));

  jackknifeDistribution<typename FitFunc::ParameterType> params(nsample,guess);
  
  jackknifeDistributionD chisq, chisq_per_dof;
  fitter.fit(params,chisq,chisq_per_dof,data);

  double dof = chisq.sample(0)/chisq_per_dof.sample(0);
  jackknifeDistributionD pvalue(nsample, [&](const int s){ return chiSquareDistribution::pvalue(dof, chisq.sample(s)); });
  
  std::cout << "Fit results\n";
  std::cout << "Params: " << params << std::endl;
  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "P-value: " << pvalue << std::endl;

  writeParamsStandard(chisq, "chisq.hdf5");
  writeParamsStandard(chisq_per_dof, "chisq_per_dof.hdf5");
  writeParamsStandard(params, "params.hdf5"); 
  writeParamsStandard(pvalue, "pvalue.hdf5");

  plot(fitfunc, data, params);
}

 void fit(const correlationFunctionType &data, const Args &args, const CMDline &cmdline){
  if(args.fitfunc == FitFuncType::FParabola){
    FitParabola fitfunc;
    fitSpecFF<FitParabola>(data,fitfunc,args,cmdline);
  }else if(args.fitfunc == FitFuncType::FDispnRelnPow4){
    FitDispnRelnPow4 fitfunc;
    fitSpecFF<FitDispnRelnPow4>(data,fitfunc,args,cmdline);

  // }else if(args.fitfunc == FitFuncType::FLinearSquareFourth){
  //   FitLinearSquareFourth fitfunc;
  //   fitSpecFF<FitLinearSquareFourth>(data,fitfunc,args,cmdline);
  // }else if(args.fitfunc == FitFuncType::FLatticeDispnReln){
  //   FitLatticeDispnReln fitfunc;
  //   fitSpecFF<FitLatticeDispnReln>(data,fitfunc,args,cmdline);
  // }else if(args.fitfunc == FitFuncType::FLatticeDispnRelnDiscE){
  //   FitLatticeDispnRelnDiscE fitfunc;
  //   fitSpecFF<FitLatticeDispnRelnDiscE>(data,fitfunc,args,cmdline);
  }else{
    assert(0);
  }
}
	 
void applyOperation(jackknifeDistribution<double> &fval, const std::string &operation){
  if(operation == "")
    return;
    
  const int nsample = fval.size();
  expressionAST AST = mathExpressionParse(operation);

  if(AST.nSymbols() != 1) error_exit(std::cout << "applyOperation expects math expression with 1 symbol ('x'), got \"" << operation << "\"\n");
  else if(!AST.containsSymbol("x")) error_exit(std::cout << "applyOperation expects math expression to be a function of 'x', got \"" << operation << "\"\n");
  
  for(int s=0;s<nsample;s++){
    AST.assignSymbol("x",fval.sample(s));
    fval.sample(s) = AST.value();
  }    
}


int main(const int argc, const char *argv[]){
  Args args;
  if(argc < 2){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);
  
  CMDline cmdline(argc,argv,2);

  const int ndata = args.data.size();
  correlationFunctionType data(ndata);
  int nsample;
  for(int i=0;i<ndata;i++){
    data.coord(i) = args.data[i].coord;
    readHDF5file(data.value(i), args.data[i].filename, args.data[i].input_idx);
    applyOperation(data.value(i),args.data[i].operation);

    if(i==0) nsample = data.value(i).size();
    else if(data.value(i).size() != nsample) error_exit(std::cout << "jackknife distribution from file " << args.data[i].filename << " has " << data.value(i).size() << " samples, expected "<< nsample << std::endl);

    std::cout << data.coord(i) << " " << data.value(i) << std::endl;
  }
  
  fit(data,args,cmdline);
  
  std::cout << "Done" << std::endl;

  return 0;
}



//Code to generate toy examples of the bootstrap p-value approach

#include<fit.h>
#include<random.h>
#include<common.h>
#include<containers.h>
#include<plot.h>
using namespace CPSfit;

class randomDataBase{
public:
  //Assumed to be thread safe
  virtual rawDataDistributionD generate(const int t, const int nsample) const = 0;
  virtual ~randomDataBase(){}
};
class randomDataGaussian: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
public:
  randomDataGaussian(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma): mu(mu), sigma(sigma){
    if(mu.size() != Lt || sigma.size() != Lt) error_exit(std::cout << "mu, sigma size must equal Lt" << std::endl);
  }
  randomDataGaussian(int Lt, double mu_all, double sigma_all): mu(Lt, mu_all), sigma(Lt, sigma_all){}
  
  rawDataDistributionD generate(const int t, const int nsample) const override{
    rawDataDistributionD out(nsample);
    gaussianRandom(out, mu[t], sigma[t], threadRNG());
    return out;
  }
};
class randomDataLogNormal: public randomDataBase{
  std::vector<double> mu;
  std::vector<double> sigma;
public:
  randomDataLogNormal(int Lt, const std::vector<double> &mu, const std::vector<double> &sigma): mu(mu), sigma(sigma){
    if(mu.size() != Lt || sigma.size() != Lt) error_exit(std::cout << "mu, sigma size must equal Lt" << std::endl);
  }
  randomDataLogNormal(int Lt, double mu_all, double sigma_all): mu(Lt, mu_all), sigma(Lt, sigma_all){}
  
  rawDataDistributionD generate(const int t, const int nsample) const override{
    rawDataDistributionD mult(nsample);
    gaussianRandom(mult, 0, 1, threadRNG());
    mult = exp(mu[t] + sigma[t]*mult);
    return mult;
  }
};


//Parameters for random data generators with just mu, sigma shared for all timeslices
#define RDATA_UNIFORM_NRMLIKE (double, mu)(double, sigma)
struct RdataUniformNrmLikeArgs{
  GENERATE_MEMBERS(RDATA_UNIFORM_NRMLIKE); 
  RdataUniformNrmLikeArgs(): mu(0.), sigma(1.){  }
};
GENERATE_PARSER( RdataUniformNrmLikeArgs, RDATA_UNIFORM_NRMLIKE);

#define RDATA_TIMEDEP_NRMLIKE (std::vector<double>, mu)(std::vector<double>, sigma)
struct RdataTimeDepNrmLikeArgs{
  GENERATE_MEMBERS(RDATA_TIMEDEP_NRMLIKE); 
  RdataTimeDepNrmLikeArgs(): mu(10,0.), sigma(10,1.){  }
};
GENERATE_PARSER( RdataTimeDepNrmLikeArgs, RDATA_TIMEDEP_NRMLIKE);




GENERATE_ENUM_AND_PARSER(DataGenStrategy, (NormalUniform)(NormalTimeDep)(LogNormalUniform) );

template<typename T>
void parseOrTemplate(T &args, const std::string &params_file, const std::string &template_file){
  if(params_file == "TEMPLATE"){
    std::ofstream of(template_file);
    (std::cout << "Outputting data generation template to '" << template_file << "' and exiting\n").flush();
    of << args;
    of.close();
    exit(0);
  } 
  parse(args, params_file);
}

std::unique_ptr<randomDataBase> dataGenStrategyFactory(DataGenStrategy strat, const std::string &params_file, const int Lt){
  if(strat == DataGenStrategy::NormalUniform){
    RdataUniformNrmLikeArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussian(Lt, args.mu, args.sigma));
  }else if(strat == DataGenStrategy::NormalTimeDep){
    RdataTimeDepNrmLikeArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataGaussian(Lt, args.mu, args.sigma));
  }else if(strat == DataGenStrategy::LogNormalUniform){
    RdataUniformNrmLikeArgs args; parseOrTemplate(args, params_file, "datagen_template.args");
    return std::unique_ptr<randomDataBase>(new randomDataLogNormal(Lt, args.mu, args.sigma));
  }else{
    error_exit(std::cout << "Invalid data generation strategy" << std::endl);
  }
}


class covMatStrategyBase{
public:
  virtual void compute(simpleSingleFitWrapper &fitter, const correlationFunction<double, rawDataDistributionD> &data) const = 0;
  virtual ~covMatStrategyBase(){};
};
class covMatStrategyCorrelated: public covMatStrategyBase{
public:
  void compute(simpleSingleFitWrapper &fitter, const correlationFunction<double, rawDataDistributionD> &data) const override{
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt);
    for(int t1=0;t1<Lt;t1++)
      for(int t2=0;t2<Lt;t2++)
	cov(t1,t2) = rawDataDistributionD::covariance( data.value(t1), data.value(t2) );
    fitter.importCovarianceMatrix(cov);
  }  
};
class covMatStrategyUncorrelated: public covMatStrategyBase{
public:
  void compute(simpleSingleFitWrapper &fitter, const correlationFunction<double, rawDataDistributionD> &data) const override{
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt, 0.);
    for(int t=0;t<Lt;t++)
      cov(t,t) = rawDataDistributionD::covariance( data.value(t), data.value(t) );
    fitter.importCovarianceMatrix(cov);
  }  
};

#define CMAT_PARAMS_CUTOFF (double, cutoff)
struct covMatStrategyCutoffArgs{
  GENERATE_MEMBERS(CMAT_PARAMS_CUTOFF); 
  covMatStrategyCutoffArgs(): cutoff(0.005){  }
};
GENERATE_PARSER( covMatStrategyCutoffArgs, CMAT_PARAMS_CUTOFF);

//cf https://arxiv.org/pdf/1101.2248.pdf
//Remove contribution of eigenvectors
class covMatStrategyCutoff: public covMatStrategyBase{
public:
  double cutoff; 

  covMatStrategyCutoff(const std::string &args_file){
    covMatStrategyCutoffArgs args;
    parseOrTemplate(args, args_file, "cov_cutoff_template.args");
    cutoff = args.cutoff;
  }

  void compute(simpleSingleFitWrapper &fitter, const correlationFunction<double, rawDataDistributionD> &data) const override{
    int Lt = data.size();
    NumericSquareMatrix<double> cov(Lt);
    for(int t1=0;t1<Lt;t1++)
      for(int t2=0;t2<Lt;t2++)
	cov(t1,t2) = rawDataDistributionD::covariance( data.value(t1), data.value(t2) );

    //Remove the contribution of eigenvalues L < cutoff from covariance matrix
    std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
    std::vector<double> evals(Lt);
    
    GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);
    
    std::cout << "Cutoff at " << cutoff << std::endl;

    NumericSquareMatrix<double> inv_cov(Lt, 0.);
    std::cout << "Evals : ";
    for(int l=0;l<Lt;l++){
      std::cout << evals[l];
      if(evals[l]>cutoff){
	std::cout << "*";

	for(int t=0;t<Lt;t++)
	  for(int u=0;u<Lt;u++)
	    inv_cov(t,u) += evecs[l](t)*evecs[l](u)*(1./evals[l]);
      }
      std::cout << " ";
    }
    std::cout << std::endl;

    std::vector<double> sigma_dummy(Lt,1.0);
    fitter.importInverseCorrelationMatrix(inv_cov, sigma_dummy);
  }  
};



GENERATE_ENUM_AND_PARSER(CovMatStrategy, (Correlated)(Uncorrelated)(Cutoff) );

//args_file need only be provided for those strategies with tunable arguments
std::unique_ptr<covMatStrategyBase> covMatStrategyFactory(CovMatStrategy strat, const std::string &args_file){
  if(strat == CovMatStrategy::Correlated){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyCorrelated);
  }else if(strat == CovMatStrategy::Uncorrelated){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyUncorrelated);
  }else if(strat == CovMatStrategy::Cutoff){
    return std::unique_ptr<covMatStrategyBase>(new covMatStrategyCutoff(args_file));
  }else{
    error_exit(std::cout << "Invalid covariance matrix strategy" << std::endl);
  }
}


GENERATE_ENUM_AND_PARSER(FitFuncType, (FConstant));

std::unique_ptr<genericFitFuncBase> fitFuncFactory(FitFuncType type){
  if(type == FitFuncType::FConstant){
    FitConstant fitfunc;
    return std::unique_ptr<genericFitFuncBase>(new simpleFitFuncWrapper<FitConstant>(fitfunc));
  }else{
    error_exit(std::cout << "Invalid fit function" << std::endl);
  }
}

//Perform some kind of pre-analysis
struct Args;

struct preAnalysisBase{
  virtual void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const = 0;
  virtual ~preAnalysisBase(){}
};
struct preAnalysisNone: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{};
};
struct preAnalysisCovMatEvals: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{}
};

GENERATE_ENUM_AND_PARSER(preAnalysisType, (None)(CovMatEvals));

std::unique_ptr<preAnalysisBase> preAnalysisFactory(preAnalysisType type){
  if(type == preAnalysisType::None){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisNone);
  }else if(type == preAnalysisType::CovMatEvals){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisCovMatEvals);
  }else{
    error_exit(std::cout << "Invalid pre-analysis type" << std::endl);
  }
}

//Note: cov_strat_params_file only required for those strategies that have tunable parameters
struct Args{
#define ARGS_MEM (int, nsample)(int, Lt)(int, ntest)(DataGenStrategy, data_strat)(FitFuncType, fitfunc)(CovMatStrategy, cov_strat)(std::string, cov_strat_params_file)(preAnalysisType, preanalysis)(MarquardtLevenbergParameters<double>, MLparams)
  GENERATE_MEMBERS(ARGS_MEM);
  
  Args(): nsample(200), Lt(30), ntest(5000), fitfunc(FitFuncType::FConstant), cov_strat(CovMatStrategy::Correlated), data_strat(DataGenStrategy::NormalUniform), cov_strat_params_file(""), preanalysis(preAnalysisType::None){
    MLparams.verbose = true;
    MLparams.lambda_factor = 1.2;
    MLparams.dampening_matrix = MLdampeningMatrix::Unit;
    MLparams.max_iter = 50000;
  }
};
GENERATE_PARSER( Args, ARGS_MEM );



struct CMDline{
  bool write_data;

  CMDline(){
    write_data = false;
  }
  CMDline(const int argc, const char** argv, const int begin = 0): CMDline(){
    setup(argc,argv,begin);
  }
  
  void setup(const int argc, const char** argv, const int begin = 0){
    const int sz = argc-begin;
    if(sz <= 0) return;

    std::vector<std::string> sargv(sz);
    for(int i=begin; i<argc; i++) sargv[i-begin] = std::string(argv[i]);

    int i = 0;
    while(i<sz){
      if(sargv[i] == "-write_data"){
	write_data = true;
	i++;
      }else{
	error_exit(std::cout << "Error: unknown argument \"" << sargv[i] << "\"\n");
      }
    }
  }
};

double estimatePvalue(const double q2, const std::vector<double> &q2s){
  double c=0, n=q2s.size();
  for(double v : q2s)
    if(v > q2) c+=1.;
  return c/n;
}

int main(const int argc, const char** argv){
  RNG.initialize(1234);
  threadRNG.initialize(5678);

  CMDline cmdline(argc,argv,3);

  Args args;
  if(argc < 3){
    std::ofstream of("template.args");
    (std::cout << "No parameter file provided: writing template to 'template.args' and exiting\n").flush();
    of << args;
    return 1;
  }    
  
  parse(args, argv[1]);
  int nsample = args.nsample;
  int Lt = args.Lt;
  int ntest = args.ntest;
  std::unique_ptr<genericFitFuncBase> ffunc = fitFuncFactory(args.fitfunc);
  std::unique_ptr<covMatStrategyBase> covgen = covMatStrategyFactory(args.cov_strat, args.cov_strat_params_file);
  std::unique_ptr<randomDataBase> dgen = dataGenStrategyFactory(args.data_strat, argv[2], Lt);

  //Run any pre-analysis
  std::unique_ptr<preAnalysisBase> preanalysis = preAnalysisFactory(args.preanalysis);
  preanalysis->run(args, *covgen, *dgen);

  //----------------------------------------------
  //Generate distribution for true data
  //----------------------------------------------
  int dof = Lt - ffunc->Nparams();

  std::vector<double> q2_dist_true(ntest);
  std::vector<double> mean_dist_true(Lt*ntest);

  std::vector<correlationFunction<double, rawDataDistributionD> > wr_data(cmdline.write_data ? ntest+1 : 0);
  
#pragma omp parallel for
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data(Lt);
    correlationFunction<double, double> data_means(Lt);
    for(int t=0;t<Lt;t++){
      data.coord(t) = t;
      data.value(t) = dgen->generate(t,nsample);

      data_means.coord(t) = t;
      data_means.value(t) = data.value(t).mean();

      mean_dist_true[t+Lt*test] = data_means.value(t);
    }

    simpleSingleFitWrapper fitter(*ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen->compute(fitter, data);

    parameterVector<double> params(1,0.);
    double q2, q2_per_dof; int dof;
    assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
    q2_dist_true[test] = q2;

    if(cmdline.write_data) wr_data[test] = data;
  }      

  //----------------------------------------------
  //Repeat with bootstrap
  //----------------------------------------------
  std::vector<double> q2_dist_boot(ntest);

  {
    correlationFunction<double, rawDataDistributionD> orig_data(Lt);
    correlationFunction<double, double> orig_data_means(Lt);

    for(int t=0;t<Lt;t++){
      orig_data.coord(t) = t;
      orig_data.value(t) = dgen->generate(t,nsample);

      orig_data_means.coord(t) = t;
      orig_data_means.value(t) = orig_data.value(t).mean();      
    }
    if(cmdline.write_data) wr_data[ntest] = orig_data;

    //Get the fit value for the parameter from the original ensemble (for recentering)
    double fit_value;
    {
      simpleSingleFitWrapper fitter(*ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
      covgen->compute(fitter, orig_data);

      parameterVector<double> params(1,0.);
      double q2, q2_per_dof; int dof;
      assert(fitter.fit(params,q2,q2_per_dof,dof, orig_data_means));
      fit_value = params[0];
    }
      
    std::vector<std::vector<int> > rtable = resampleTable(RNG, nsample, ntest);
   
#pragma omp parallel for
    for(int test=0;test<ntest;test++){
      correlationFunction<double, rawDataDistributionD> data(Lt);
      correlationFunction<double, double> data_means(Lt);
      for(int t=0;t<Lt;t++){
	rawDataDistributionD &dd = data.value(t);
	dd.resize(nsample);
	for(int s=0;s<nsample;s++){
	  dd.sample(s) = orig_data.value(t).sample(rtable[test][s]) //resample
	    + fit_value - orig_data_means.value(t); //recenter
	} 

	data.coord(t) = t;
	data_means.coord(t) = t;
	data_means.value(t) = dd.mean();
      }
      simpleSingleFitWrapper fitter(*ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
      covgen->compute(fitter, data);

      parameterVector<double> params(1,0.);
      double q2, q2_per_dof; int dof;
      assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
      q2_dist_boot[test] = q2;
    }      

    if(cmdline.write_data){ //write resample table
      std::ofstream f("rtable.dat");
      for(int e=0;e<ntest;e++){
	for(int s=0;s<nsample;s++)
	  f << rtable[e][s] << " ";
	f << std::endl;
      }
    }
  }//bootstrap analysis

  if(cmdline.write_data){ //write data
    std::ofstream f("data.dat");
    f << std::setprecision(16);
    for(int e=0;e<ntest+1;e++){ //original ensemble for bootstrap analysis is the last (e==ntest) ensemble
      for(int t=0;t<Lt;t++){
	f << e << " " << t;
	for(int s=0;s<nsample;s++)
	  f << " " << wr_data[e].value(t).sample(s);
	f << std::endl;
      }
    }
  }

  //----------------------------------------------
  //Plot results
  //----------------------------------------------

  //Generate histograms
  {
    MatPlotLibScriptGenerate plot;
    struct acc{
      const std::vector<double> &d;
      acc(const std::vector<double> &d): d(d){}
      double y(const int i) const{ return d[i]; }
      int size() const{ return d.size(); }
    };
    typename MatPlotLibScriptGenerate::kwargsType kwargs;
    kwargs["density"] = true;
    kwargs["alpha"] = 0.4;
    kwargs["bins"] = 60;
    auto htrue = plot.histogram(acc(q2_dist_true),kwargs,"true");
    plot.setLegend(htrue, R"(${\\rm true}$)");

    kwargs["color"] = 'c';
    auto hboot = plot.histogram(acc(q2_dist_boot),kwargs,"boot");
    plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

    double q2_max = *std::max_element(q2_dist_true.begin(),q2_dist_true.end());
    int npt = 200;
  
    struct T2data : public CurveDataAccessorBase<double>{
      int npt;
      double delta;
      int p;
      int n;

      T2data(int npt, double delta, int p, int n): npt(npt), delta(delta), p(p), n(n){}

      double x(const int i) const override{ return i*delta; }
      double y(const int i) const override{ return TsquareDistribution::PDF(x(i),p,n); }
      int size() const override{ return npt; }
    };
    kwargs.clear();
    kwargs["color"] = "b";
    auto hT2 = plot.errorBand(T2data(npt, q2_max/(npt-1), dof, nsample-1),kwargs,"T2");
    plot.setLegend(hT2, R"($T^2$)");

    struct ChisqData : public CurveDataAccessorBase<double>{
      int npt;
      double delta;
      int p;

      ChisqData(int npt, double delta, int p): npt(npt), delta(delta), p(p){}

      double x(const int i) const override{ return 1e-3 + i*delta; }
      double y(const int i) const override{ return chiSquareDistribution::PDF(x(i),p); }
      int size() const override{ return npt; }
    };
    kwargs.clear();
    kwargs["color"] = "g";
    auto hchi2 = plot.errorBand(ChisqData(npt, q2_max/(npt-1), dof),kwargs,"chi2");
    plot.setLegend(hchi2, R"($\\chi^2$)");

    plot.setXlabel(R"($q^2$)");
    plot.setYlabel(R"(${\cal F}(q^2)$)");
    plot.createLegend();

    plot.write("q2_dist.py","q2_dist.pdf");
  }

  //Plot true and bootstrap p-value as a function of q2
  {
    int npt = 200;
    double q2_max = *std::max_element(q2_dist_true.begin(),q2_dist_true.end()) * 1.10;

    std::vector<double> q2vals(npt);
    std::vector<double> ptrue(npt);
    std::vector<double> pboot(npt);
    std::vector<double> pT2(npt);
    std::vector<double> pchi2(npt);
    
    double dq2 = q2_max/(npt-1);

    for(int i=0;i<npt;i++){
      double q2 = dq2*i;
      q2vals[i] = q2;
      ptrue[i] = estimatePvalue(q2, q2_dist_true);
      pboot[i] = estimatePvalue(q2, q2_dist_boot);
      pchi2[i] = chiSquareDistribution::pvalue(dof, q2);
      pT2[i] = TsquareDistribution::pvalue(q2, dof, nsample-1);
    }

    struct acc : public CurveDataAccessorBase<double>{
      const std::vector<double> &xx;
      const std::vector<double> &yy;
      acc(const std::vector<double> &xx, const std::vector<double> &yy): xx(xx), yy(yy){}

      double x(const int i) const override{ return xx[i]; }
      double y(const int i) const override{ return yy[i]; }
      int size() const override{ return xx.size(); }
    };

    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;

      kwargs["color"] = "b";
      auto htrue = plot.errorBand(acc(q2vals,ptrue), kwargs, "true");
      plot.setLegend(htrue, R"(${\\rm true}$)");
      kwargs["color"] = "r";
      auto hboot = plot.errorBand(acc(q2vals,pboot), kwargs, "boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");
      kwargs["color"] = "g";
      auto hT2 = plot.errorBand(acc(q2vals,pT2), kwargs, "T2");
      plot.setLegend(hT2, R"($T^2$)");
      kwargs["color"] = "m";
      auto hchi2 = plot.errorBand(acc(q2vals,pchi2), kwargs, "chi2");
      plot.setLegend(hchi2, R"($\\chi^2$)");
      plot.setXlabel(R"($q^2$)");
      plot.setYlabel(R"(${\rm p-value}$)");
      plot.createLegend();
      plot.write("pvalue.py","pvalue.pdf");
    }
    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;

      kwargs["color"] = "b";
      auto htrue = plot.errorBand(acc(ptrue,ptrue), kwargs, "true");
      plot.setLegend(htrue, R"(${\\rm true}$)");
      kwargs["color"] = "r";
      auto hboot = plot.errorBand(acc(ptrue,pboot), kwargs, "boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");
      kwargs["color"] = "g";
      auto hT2 = plot.errorBand(acc(ptrue,pT2), kwargs, "T2");
      plot.setLegend(hT2, R"($T^2$)");
      kwargs["color"] = "m";
      auto hchi2 = plot.errorBand(acc(ptrue,pchi2), kwargs, "chi2");
      plot.setLegend(hchi2, R"($\\chi^2}$)");
      plot.setXlabel(R"(${\rm p-value\ (true)}$)");
      plot.setYlabel(R"(${\rm p-value\ (est.)}$)");
      plot.createLegend();
      plot.write("pest_v_ptrue.py","pest_v_ptrue.pdf");
    }



  }


  //Plot distribution of means
  {
    MatPlotLibScriptGenerate plot_mean;
    typename MatPlotLibScriptGenerate::kwargsType kwargs;
    kwargs["density"] = true;
    kwargs["bins"] = 60;
    struct acc{
      const std::vector<double> &d;
      acc(const std::vector<double> &d): d(d){}
      double y(const int i) const{ return d[i]; }
      int size() const{ return d.size(); }
    };
    plot_mean.histogram(acc(mean_dist_true),kwargs);
  
    plot_mean.write("means_true.py","means_true.pdf");
  }


  std::cout << "Done\n";
  return 0;
}


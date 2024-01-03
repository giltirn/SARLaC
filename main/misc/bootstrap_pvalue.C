//Code to generate toy examples of the bootstrap p-value approach
#include<tuple>
#include<fit.h>
#include<random.h>
#include<common.h>
#include<containers.h>
#include<plot.h>
using namespace SARLaC;

#include <bootstrap_pvalue/cmdline.h>
#include <bootstrap_pvalue/enum.h>
#include <bootstrap_pvalue/args.h>
#include <bootstrap_pvalue/utils.h>
#include <bootstrap_pvalue/data_gen.h>
#include <bootstrap_pvalue/cov_mat.h>
#include <bootstrap_pvalue/fitfunc.h>
#include <bootstrap_pvalue/preanalysis.h>


//data_means_out :  output unrecentered data means
void bootstrapAnalyze(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Args &args, const std::string &write_rtable = "", std::vector<correlationFunction<double, double> > *data_means_out = nullptr){
  int nsample = args.nsample;
  int Lt = args.Lt;
  int nblock = nsample / args.block_size;
  int nsample_reduced = nblock * args.block_size;
  int ntest = args.ntest;
  assert(q2_into.size() == ntest);

  correlationFunction<double, double> orig_data_means(Lt);

  for(int t=0;t<Lt;t++){
    orig_data_means.coord(t) = t;
    orig_data_means.value(t) = orig_data.value(t).mean();      
  }

  //Get the fit value for the parameter from the original ensemble (for recentering
  parameterVector<double> orig_fit_params(ffunc.Nparams(),0.);
  {
    simpleSingleFitWrapper fitter(ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen.compute(fitter, orig_data);
    
    double q2, q2_per_dof; int dof;
    assert(fitter.fit(orig_fit_params,q2,q2_per_dof,dof, orig_data_means));
  }    

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
  assert(rtable.size() == ntest && rtable[0].size() == nsample_reduced);
   
  if(data_means_out) data_means_out->resize(ntest);

#pragma omp parallel for
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data(Lt);
    correlationFunction<double, double> data_means(Lt), data_means_unrecentered(Lt);
    for(int t=0;t<Lt;t++){
      double fit_value = ffunc.value(generalContainer(double(t)), orig_fit_params);
      data.coord(t) = data_means.coord(t) = t;
      rawDataDistributionD shift(nsample_reduced, fit_value - orig_data_means.value(t));

      rawDataDistributionD &dd = data.value(t);
      dd = resampledEnsemble(orig_data.value(t), test, rtable);

      data_means_unrecentered.value(t) = dd.mean();

      //Recenter
      dd = dd + shift;
      data_means.value(t) = dd.mean();
    }
    simpleSingleFitWrapper fitter(ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen.compute(fitter, data);

    parameterVector<double> params(ffunc.Nparams(),0.);
    double q2, q2_per_dof; int dof;
    assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
    q2_into[test] = q2;

    if(data_means_out) (*data_means_out)[test] = std::move(data_means_unrecentered);
  }

  if(write_rtable.size()){
    std::ofstream f(write_rtable);
    for(int e=0;e<ntest;e++){
      for(int s=0;s<nblock;s++)
	f << rtable[e][s] << " ";
      f << std::endl;
    }
  }
} 


void bootstrapAnalyzeResiduals(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Args &args, const std::string &write_rtable = ""){
  int nsample = args.nsample;
  int Lt = args.Lt;
  int nblock = nsample / args.block_size;
  int nsample_reduced = nblock * args.block_size;
  int ntest = args.ntest;
  assert(q2_into.size() == ntest);

  correlationFunction<double, double> orig_data_means(Lt);

  for(int t=0;t<Lt;t++){
    orig_data_means.coord(t) = t;
    orig_data_means.value(t) = orig_data.value(t).mean();      
  }

  //Get the fit value for the parameter from the original ensemble (for recentering)
  parameterVector<double> orig_fit_params(ffunc.Nparams(),0.);
  {
    simpleSingleFitWrapper fitter(ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen.compute(fitter, orig_data);
    double q2, q2_per_dof; int dof;
    assert(fitter.fit(orig_fit_params,q2,q2_per_dof,dof, orig_data_means));
  }    

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
  assert(rtable.size() == ntest && rtable[0].size() == nsample_reduced);
   
  correlationFunction<double, rawDataDistributionD> orig_data_resids(Lt);
  for(int t=0;t<Lt;t++){
    orig_data_resids.coord(t) = t;
    double fitval = ffunc.value(generalContainer(double(t)), orig_fit_params);
    orig_data_resids.value(t) = orig_data.value(t) - rawDataDistributionD(nsample, fitval);
  }

#pragma omp parallel for
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data(Lt);
    correlationFunction<double, double> data_means(Lt);
    for(int t=0;t<Lt;t++){
      rawDataDistributionD &dd = data.value(t);
      dd.resize(nsample_reduced);
      double shift = orig_data_resids.value(t).mean();
      for(int s=0;s<nsample_reduced;s++){
	dd.sample(s) = orig_data_resids.value(t).sample(rtable[test][s])
	  - shift; //recenter
      } 

      data.coord(t) = data_means.coord(t) = t;
      data_means.value(t) = dd.mean();
    }
    NumericSquareMatrix<double> cov = covgen.compute(data);
    NumericSquareMatrix<double> inv_cov(cov); svd_inverse(inv_cov, cov);

    double &q2 = q2_into[test]; 
    q2=0.;
    for(int t=0;t<Lt;t++)
      for(int u=0;u<Lt;u++)
	q2 += data_means.value(t) * inv_cov(t,u) * data_means.value(u);
  }

  if(write_rtable.size()){
    std::ofstream f(write_rtable);
    for(int e=0;e<ntest;e++){
      for(int s=0;s<nblock;s++)
	f << rtable[e][s] << " ";
      f << std::endl;
    }
  }
} 



int main(const int argc, const char** argv){
  CMDline cmdline(argc,argv,3);

  RNG.initialize(1234);
  threadRNG.initialize(5678);

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
  if(args.preanalysis.size() != args.preanalysis_params_file.size()) error_exit(std::cout << "Require as many preanalysis params files as there are preanalyses specified" << std::endl);

  for(int i=0;i<args.preanalysis.size();i++){
    std::unique_ptr<preAnalysisBase> preanalysis = preAnalysisFactory(args.preanalysis[i]);
    preanalysis->run(args, args.preanalysis_params_file[i], *covgen, *dgen, *ffunc);
  }
  if(cmdline.exit_after_preanalysis) return 0;

  //----------------------------------------------
  //Generate distribution for true data
  //----------------------------------------------
  int dof = Lt - ffunc->Nparams();

  std::vector<double> q2_dist_true(ntest);
  std::vector<double> mean_dist_true(Lt*ntest);

  std::vector<correlationFunction<double, rawDataDistributionD> > wr_data(cmdline.write_data ? ntest+args.norig_ens : 0);
  
#pragma omp parallel for schedule(static)
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data = dgen->generate(Lt,nsample);
    correlationFunction<double, double> data_means(Lt);
    for(int t=0;t<Lt;t++){
      data_means.coord(t) = t;
      data_means.value(t) = data.value(t).mean();

      mean_dist_true[t+Lt*test] = data_means.value(t);
    }

    simpleSingleFitWrapper fitter(*ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen->compute(fitter, data);

    parameterVector<double> params(ffunc->Nparams(),0.);
    double q2, q2_per_dof; int dof;
    assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
    q2_dist_true[test] = q2;

    if(cmdline.write_data) wr_data[test] = data;
  }      

  //---------------------------------------------------------------
  //Repeat with bootstrap for norig_ens separate original ensembles
  //---------------------------------------------------------------
  std::vector< std::vector<double> > q2_dist_boot(args.norig_ens, std::vector<double>(ntest));
  for(int o=0;o<args.norig_ens;o++){
    correlationFunction<double, rawDataDistributionD> orig_data = dgen->generate(Lt,nsample);
    if(cmdline.write_data) wr_data[ntest+o] = orig_data;
    bootstrapAnalyze(q2_dist_boot[o], orig_data, *covgen, *ffunc, args, cmdline.write_data ? "rtable."+std::to_string(o)+".dat" : "");
  }

  if(cmdline.write_data){ //write data
    std::ofstream f("data.dat");
    f << std::setprecision(16);
    for(int e=0;e<ntest+args.norig_ens;e++){ //original ensemble for bootstrap analysis is the last norig_ens ensemble (e==ntest+o  for o=0..norig_ens-1) 
      for(int t=0;t<Lt;t++){
	f << e << " " << t;
	for(int s=0;s<nsample;s++)
	  f << " " << wr_data[e].value(t).sample(s);
	f << std::endl;
      }
    }
  }



  //---------------------------------------------------------------
  //Repeat with bootstrap for norig_ens bootstrap resampled ensembles in place of real original ensembles for error estimation
  //---------------------------------------------------------------
  std::vector< std::vector<double> > q2_dist_dbl_boot(args.norig_ens, std::vector<double>(ntest));  
  correlationFunction<double, rawDataDistributionD> dbl_boot_base_data = dgen->generate(Lt,nsample);

  std::vector<std::vector<int> > dbl_boot_rtable = generateResampleTable(nsample, args.norig_ens, args.bootstrap_strat, args.block_size, threadRNG);
  int nblock = nsample / args.block_size;
  int nsample_reduced = nblock * args.block_size;

  assert(dbl_boot_rtable.size() == args.norig_ens && dbl_boot_rtable[0].size() == nsample_reduced);

  for(int o=0;o<args.norig_ens;o++){
    correlationFunction<double, rawDataDistributionD> orig_data(Lt);
    for(int t=0;t<Lt;t++){
      orig_data.coord(t) = t;
      rawDataDistributionD &dd = orig_data.value(t);
      dd.resize(nsample_reduced);
      for(int s=0;s<nsample_reduced;s++)
	dd.sample(s) = dbl_boot_base_data.value(t).sample(dbl_boot_rtable[o][s]);
    }

    bootstrapAnalyze(q2_dist_dbl_boot[o], orig_data, *covgen, *ffunc, args);
  }



  //---------------------------------------------------------------
  //Repeat with bootstrap residuals version for norig_ens separate original ensembles
  //---------------------------------------------------------------
  std::vector< std::vector<double> > q2_resid_dist_boot(args.norig_ens, std::vector<double>(ntest));

  for(int o=0;o<args.norig_ens;o++){
    correlationFunction<double, rawDataDistributionD> orig_data = dgen->generate(Lt,nsample);
    bootstrapAnalyzeResiduals(q2_resid_dist_boot[o], orig_data, *covgen, *ffunc, args);
  }

  //----------------------------------------------
  //Plot results
  //----------------------------------------------

  //Generate histograms (only for first bootstrap original ensemble)
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
    auto hboot = plot.histogram(acc(q2_dist_boot[0]),kwargs,"boot");
    plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

    kwargs["color"] = "tab:gray";
    auto hboot_resid = plot.histogram(acc(q2_resid_dist_boot[0]),kwargs,"boot_resid");
    plot.setLegend(hboot_resid, R"(${\\rm bootstrap resid}$)");


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
      double y(const int i) const override{ return chiSquareDistribution::PDF(p,x(i)); }
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
    std::vector<rawDataDistributionD> pboot(npt, rawDataDistributionD(args.norig_ens) ); //variation over actual original ensembles
    std::vector<rawDataDistributionD> pdbl_boot(npt, rawDataDistributionD(args.norig_ens) ); //variation over bootstrap resampled ensembles in place of original ensembles
    std::vector<double> pboot_ens0(npt); //bootstrap p-value from just first original ensemble
    std::vector<rawDataDistributionD> pboot_resid(npt, rawDataDistributionD(args.norig_ens) ); //variation over actual original ensembles for residuals version
    std::vector<double> pT2(npt);
    std::vector<double> pchi2(npt);
    
    double dq2 = q2_max/(npt-1);

    for(int i=0;i<npt;i++){
      double q2 = dq2*i;
      q2vals[i] = q2;
      ptrue[i] = estimatePvalue(q2, q2_dist_true);
      pchi2[i] = chiSquareDistribution::pvalue(dof, q2);
      pT2[i] = TsquareDistribution::pvalue(q2, dof, nsample-1);
      for(int o=0;o<args.norig_ens;o++){
	pboot[i].sample(o) = estimatePvalue(q2, q2_dist_boot[o]);
	pdbl_boot[i].sample(o) = estimatePvalue(q2, q2_dist_dbl_boot[o]);
	pboot_resid[i].sample(o) = estimatePvalue(q2, q2_resid_dist_boot[o]);
      }
      pboot_ens0[i] = pboot[i].sample(0);
    }

    struct acc : public CurveDataAccessorBase<double>{
      const std::vector<double> &xx;
      const std::vector<double> &yy;
      acc(const std::vector<double> &xx, const std::vector<double> &yy): xx(xx), yy(yy){}

      double x(const int i) const override{ return xx[i]; }
      double y(const int i) const override{ return yy[i]; }
      int size() const override{ return xx.size(); }
    };

    class acc_werr{
      const std::vector<double> &xx;
      const std::vector<rawDataDistributionD> &yy;

    public:
      acc_werr(const std::vector<double> &x, const std::vector<rawDataDistributionD> &y): xx(x), yy(y){}

      double x(const int i) const{ return xx[i]; }
      double y(const int i) const{ return yy[i].mean(); }
      double dxm(const int i) const{ return 0; }
      double dxp(const int i) const{ return 0; }
      double dym(const int i) const{ return yy[i].standardDeviation(); }
      double dyp(const int i) const{ return yy[i].standardDeviation(); }

      double upper(const int i) const{ return yy[i].mean() + yy[i].standardDeviation(); }
      double lower(const int i) const{ return yy[i].mean() - yy[i].standardDeviation(); }
      
      int size() const{ return xx.size(); }
    };

    class acc_wsep_err{
      const std::vector<double> &xx;
      const std::vector<double> &yy;
      const std::vector<rawDataDistributionD> &yy_err;

    public:
      acc_wsep_err(const std::vector<double> &x, const std::vector<double> &y,  const std::vector<rawDataDistributionD> &y_err): xx(x), yy(y), yy_err(y_err){}

      double x(const int i) const{ return xx[i]; }
      double y(const int i) const{ return yy[i]; }
      double dxm(const int i) const{ return 0; }
      double dxp(const int i) const{ return 0; }
      double dym(const int i) const{ return yy_err[i].standardDeviation(); }
      double dyp(const int i) const{ return yy_err[i].standardDeviation(); }

      double upper(const int i) const{ return yy[i] + yy_err[i].standardDeviation(); }
      double lower(const int i) const{ return yy[i] - yy_err[i].standardDeviation(); }
      
      int size() const{ return xx.size(); }
    };


    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;

      kwargs["color"] = "b";
      auto htrue = plot.errorBand(acc(q2vals,ptrue), kwargs, "true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = "r";
      auto kwb = kwargs; kwb["alpha"]=0.3;
      auto hboot = plot.errorBand(acc_werr(q2vals,pboot), kwb, "boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      kwb["color"] = "tab:pink";
      auto hdbl_boot = plot.errorBand(acc_wsep_err(q2vals,pboot_ens0,pdbl_boot), kwb, "dbl_boot");
      plot.setLegend(hdbl_boot, R"(${\\rm dbl. bootstrap}$)");

      kwb["color"] = "tab:gray";
      auto hboot_resid = plot.errorBand(acc_werr(q2vals,pboot_resid), kwb, "boot_resid");
      plot.setLegend(hboot_resid, R"(${\\rm bootstrap resid.}$)");

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
      auto kwb = kwargs; kwb["alpha"]=0.3;
      auto hboot = plot.errorBand(acc_werr(ptrue,pboot), kwb, "boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      kwb["color"] = "tab:pink";
      auto hdbl_boot = plot.errorBand(acc_wsep_err(ptrue,pboot_ens0,pdbl_boot), kwb, "dbl_boot");
      plot.setLegend(hdbl_boot, R"(${\\rm dbl. bootstrap}$)");

      kwb["color"] = "tab:gray";
      auto hboot_resid = plot.errorBand(acc_werr(ptrue,pboot_resid), kwb, "boot_resid");
      plot.setLegend(hboot_resid, R"(${\\rm bootstrap resid.}$)");

      kwargs["color"] = "g";
      auto hT2 = plot.errorBand(acc(ptrue,pT2), kwargs, "T2");
      plot.setLegend(hT2, R"($T^2$)");

      kwargs["color"] = "m";
      auto hchi2 = plot.errorBand(acc(ptrue,pchi2), kwargs, "chi2");

      plot.setLegend(hchi2, R"($\\chi^2}$)");
      plot.setXlabel(R"(${\rm p-value\ (true)}$)");
      plot.setYlabel(R"(${\rm p-value\ (est.)}$)");

      kwargs.clear();
      kwargs["loc"] = "upper left";
      plot.createLegend(kwargs);
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


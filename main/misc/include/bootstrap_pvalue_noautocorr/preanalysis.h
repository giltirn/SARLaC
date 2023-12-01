#pragma once

struct preAnalysisBase{
  virtual void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const = 0;
  virtual ~preAnalysisBase(){}
};
struct preAnalysisNone: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{};
};
struct preAnalysisCovMatEvals: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{
    int nsample = args.nsample;
    int Lt = args.Lt;
    int ntest = args.ntest;
    
    //Generate data for the evals of the true covariance matrix
    std::vector<double> evals_true(ntest*Lt);
    std::vector<double> evals_lo_true(ntest);
    std::vector<double> evals_hi_true(ntest);
#pragma omp parallel for
    for(int test=0;test<ntest;test++){
      correlationFunction<double, rawDataDistributionD> data = datagen.generate(Lt,nsample);

      NumericSquareMatrix<double> cov = covgen.compute(data);

      std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
      std::vector<double> evals(Lt);
    
      GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);

      for(int t=0;t<Lt;t++) evals_true[t+Lt*test] = evals[t];
      evals_hi_true[test] = evals[0];
      evals_lo_true[test] = evals[Lt-1];
    }

    correlationFunction<double, rawDataDistributionD> orig_data = datagen.generate(Lt,nsample);

    std::vector<std::vector<int> > rtable = resampleTable(RNG, nsample, ntest);

    std::vector<double> evals_boot(ntest*Lt);
    std::vector<double> evals_lo_boot(ntest);
    std::vector<double> evals_hi_boot(ntest);

#pragma omp parallel for
    for(int test=0;test<ntest;test++){
      correlationFunction<double, rawDataDistributionD> data(Lt);
      for(int t=0;t<Lt;t++){
	rawDataDistributionD &dd = data.value(t);
	dd.resize(nsample);
	for(int s=0;s<nsample;s++){
	  dd.sample(s) = orig_data.value(t).sample(rtable[test][s]); //recentering not needed for covariance matrix
	} 
	data.coord(t) = t;
      }
      NumericSquareMatrix<double> cov = covgen.compute(data);

      std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
      std::vector<double> evals(Lt);
    
      GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);

      for(int t=0;t<Lt;t++) evals_boot[t+Lt*test] = evals[t];
      evals_hi_boot[test] = evals[0];
      evals_lo_boot[test] = evals[Lt-1];
    }      
    
    //Generate histograms
    struct acc{
      const std::vector<double> &d;
      acc(const std::vector<double> &d): d(d){}
      double y(const int i) const{ return d[i]; }
      int size() const{ return d.size(); }
    };

    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["density"] = true;
      kwargs["alpha"] = 0.4;
      kwargs["bins"] = 60;
      auto htrue = plot.histogram(acc(evals_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(acc(evals_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($\lambda$)");
      plot.setYlabel(R"($\rho(\lambda)$)");
      plot.createLegend();

      plot.write("eval_density.py","eval_density.pdf");
    }

    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["alpha"] = 0.4;
      kwargs["bins"] = 60;
      auto htrue = plot.histogram(acc(evals_hi_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(acc(evals_hi_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($\lambda_{\rm hi}$)");
      plot.setYlabel(R"($f(\lambda_{\rm hi})$)");
      plot.createLegend();

      plot.write("eval_hi_hist.py","eval_hi_hist.pdf");
    }
    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["alpha"] = 0.4;
      kwargs["bins"] = 60;
      auto htrue = plot.histogram(acc(evals_lo_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(acc(evals_lo_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($\lambda_{\rm lo}$)");
      plot.setYlabel(R"($f(\lambda_{\rm lo})$)");
      plot.createLegend();

      plot.write("eval_lo_hist.py","eval_lo_hist.pdf");
    }
    struct accNrm{
      double N;
      const std::vector<double> &d;
      accNrm(double N, const std::vector<double> &d): N(N), d(d){}
      double y(const int i) const{ return N*d[i]; }
      int size() const{ return d.size(); }
    };
    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["alpha"] = 0.4;
      kwargs["bins"] = 60;
      auto htrue = plot.histogram(accNrm(nsample,evals_hi_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(accNrm(nsample,evals_hi_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($N\lambda_{\rm hi}$)");
      plot.setYlabel(R"($f(N\lambda_{\rm hi})$)");
      plot.createLegend();

      plot.write("eval_hi_nrm_hist.py","eval_hi_nrm_hist.pdf");
    }
    {
      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["alpha"] = 0.4;
      kwargs["bins"] = 60;
      auto htrue = plot.histogram(accNrm(nsample,evals_lo_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(accNrm(nsample,evals_lo_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($N\lambda_{\rm lo}$)");
      plot.setYlabel(R"($f(N\lambda_{\rm lo})$)");
      plot.createLegend();

      plot.write("eval_lo_nrm_hist.py","eval_lo_nrm_hist.pdf");
    }

    {
      struct accCondNum{
	double N;
	const std::vector<double> &hi;
	const std::vector<double> &lo;
	accCondNum(double N, const std::vector<double> &hi, const std::vector<double> &lo): N(N), hi(hi), lo(lo){}
	double y(const int i) const{ return hi[i]/lo[i]; }
	int size() const{ return hi.size(); }
      };

      MatPlotLibScriptGenerate plot;
      typename MatPlotLibScriptGenerate::kwargsType kwargs;
      kwargs["alpha"] = 0.4;
      kwargs["bins"] = 60;
      auto htrue = plot.histogram(accCondNum(nsample,evals_hi_true,evals_lo_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(accCondNum(nsample,evals_hi_boot,evals_lo_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($\lambda_{\rm hi}/\lambda_{\rm lo}$)");
      plot.setYlabel(R"($f(\lambda_{\rm hi}/\lambda_{\rm lo})$)");
      plot.createLegend();

      plot.write("cond_num_hist.py","cond_num_hist.pdf");
    }

  }
};

struct preAnalysisCorrMatEvals: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{
    int nsample = args.nsample;
    int Lt = args.Lt;
    int ntest = args.ntest;
    
    //Generate data for the evals of the true covariance matrix
    std::vector<double> evals_true(ntest*Lt);
#pragma omp parallel for
    for(int test=0;test<ntest;test++){
      correlationFunction<double, rawDataDistributionD> data = datagen.generate(Lt,nsample);

      NumericSquareMatrix<double> cov = covgen.compute(data);

      std::vector<double> sigma(Lt);
      for(int t=0;t<Lt;t++) sigma[t] = sqrt(cov(t,t));

      for(int t=0;t<Lt;t++)
	for(int u=0;u<Lt;u++) cov(t,u) = cov(t,u)/sigma[t]/sigma[u];

      std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
      std::vector<double> evals(Lt);
    
      GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);

      for(int t=0;t<Lt;t++) evals_true[t+Lt*test] = evals[t];
    }

    correlationFunction<double, rawDataDistributionD> orig_data = datagen.generate(Lt,nsample);

    std::vector<std::vector<int> > rtable = resampleTable(RNG, nsample, ntest);

    std::vector<double> evals_boot(ntest*Lt);
#pragma omp parallel for
    for(int test=0;test<ntest;test++){
      correlationFunction<double, rawDataDistributionD> data(Lt);
      for(int t=0;t<Lt;t++){
	rawDataDistributionD &dd = data.value(t);
	dd.resize(nsample);
	for(int s=0;s<nsample;s++){
	  dd.sample(s) = orig_data.value(t).sample(rtable[test][s]); //recentering not needed for covariance matrix
	} 
	data.coord(t) = t;
      }
      NumericSquareMatrix<double> cov = covgen.compute(data);

      std::vector<double> sigma(Lt);
      for(int t=0;t<Lt;t++) sigma[t] = sqrt(cov(t,t));

      for(int t=0;t<Lt;t++)
	for(int u=0;u<Lt;u++) cov(t,u) = cov(t,u)/sigma[t]/sigma[u];

      std::vector< NumericVector<double> > evecs(Lt, NumericVector<double>(Lt) );
      std::vector<double> evals(Lt);
    
      GSLsymmEigenSolver< NumericVector<double>, NumericSquareMatrix<double> >::symmetricMatrixSolve(evecs, evals, cov);

      for(int t=0;t<Lt;t++) evals_boot[t+Lt*test] = evals[t];
    }      
    
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
      auto htrue = plot.histogram(acc(evals_true),kwargs,"true");
      plot.setLegend(htrue, R"(${\\rm true}$)");

      kwargs["color"] = 'c';
      auto hboot = plot.histogram(acc(evals_boot),kwargs,"boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      plot.setXlabel(R"($\lambda$)");
      plot.setYlabel(R"($\rho(\lambda)$)");
      plot.createLegend();

      plot.write("eval_density.py","eval_density.pdf");
    }
  }
};

struct preAnalysisStandardError: public preAnalysisBase{
  static void distProp(double &mean, double &std_dev, const std::vector<parameterVector<double> > &vals, int param){
    double s=0,s2=0;
    for(auto const &v : vals){
      s += v[param];
      s2 += v[param]*v[param];
    }
    mean = s/vals.size();
    std_dev = sqrt( s2/vals.size() - mean*mean );
  }

  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{   
    //----------------------------------------------
    //Generate distribution for true data
    //----------------------------------------------
    int nsample = args.nsample;
    int Lt = args.Lt;
    int ntest = args.ntest;
    int dof = Lt - fitfunc.Nparams();

    std::vector<parameterVector<double> > fitval_dist_true(ntest);
  
#pragma omp parallel for
    for(int test=0;test<ntest;test++){
      correlationFunction<double, rawDataDistributionD> data = datagen.generate(Lt,nsample);
      correlationFunction<double, double> data_means(Lt);
      for(int t=0;t<Lt;t++){
	data_means.coord(t) = t;
	data_means.value(t) = data.value(t).mean();
      }

      simpleSingleFitWrapper fitter(fitfunc, MinimizerType::MarquardtLevenberg, args.MLparams);
      covgen.compute(fitter, data);

      parameterVector<double> params(fitfunc.Nparams(),0.);
      double q2, q2_per_dof; int dof;
      assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
      fitval_dist_true[test] = params;
    }      

    //---------------------------------------------------------------
    //Repeat with bootstrap for norig_ens separate original ensembles
    //---------------------------------------------------------------
    std::vector< std::vector<parameterVector<double> > > fitval_dist_boot(args.norig_ens, std::vector<parameterVector<double> >(ntest));

    for(int o=0;o<args.norig_ens;o++){
      correlationFunction<double, rawDataDistributionD> orig_data = datagen.generate(Lt,nsample);     
      int nblock = nsample / args.block_size;
      int nsample_reduced = nblock * args.block_size;
      std::vector<std::vector<int> > rtable = resampleTable(threadRNG, nblock, ntest);
   
#pragma omp parallel for
      for(int test=0;test<ntest;test++){
	correlationFunction<double, rawDataDistributionD> data(Lt);
	correlationFunction<double, double> data_means(Lt);
	for(int t=0;t<Lt;t++){
	  rawDataDistributionD &dd = data.value(t);
	  dd.resize(nsample_reduced);
	  for(int b=0;b<nblock;b++){
	    for(int bs=0;bs<args.block_size;bs++)
	      dd.sample(bs + args.block_size * b) = orig_data.value(t).sample(bs + args.block_size * rtable[test][b]); //block resample
	  } 

	  data.coord(t) = t;
	  data_means.coord(t) = t;
	  data_means.value(t) = dd.mean();
	}
	simpleSingleFitWrapper fitter(fitfunc, MinimizerType::MarquardtLevenberg, args.MLparams);
	covgen.compute(fitter, data);

	parameterVector<double> params(fitfunc.Nparams(),0.);
	double q2, q2_per_dof; int dof;
	assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
	fitval_dist_boot[o][test] = params;
      }      
    }//bootstrap analysis

    for(int p=0;p<fitfunc.Nparams();p++){
      double mu_true, std_dev_true;
      distProp(mu_true,std_dev_true,fitval_dist_true,p);
      std::vector<double> mu_boot(args.norig_ens), std_dev_boot(args.norig_ens);
      for(int i=0;i<args.norig_ens;i++)  distProp(mu_boot[i],std_dev_boot[i],fitval_dist_boot[i],p);

      std::cout << "Param " << p << " true fit results: " << mu_true << " +- " << std_dev_true << std::endl;
      double s_mu=0, s2_mu=0, s_se=0, s2_se=0;
      for(int i=0;i<args.norig_ens;i++){ 
	std::cout << "Param " << p << " bootstrap orig ens " << i << " fit results: " << mu_boot[i] << " +- " << std_dev_boot[i] << std::endl;
	
	s_mu += mu_boot[i];
	s2_mu += mu_boot[i]*mu_boot[i];
	
	s_se += std_dev_boot[i];
	s2_se += std_dev_boot[i]*std_dev_boot[i];
      }
      double mu_boot_origens_avg = s_mu/args.norig_ens;
      double mu_boot_origens_stddev = sqrt( s2_mu/args.norig_ens - mu_boot_origens_avg*mu_boot_origens_avg );
    
      double se_boot_origens_avg = s_se/args.norig_ens;
      double se_boot_origens_stddev = sqrt( s2_se/args.norig_ens - se_boot_origens_avg*se_boot_origens_avg );
    
      std::cout << "Param " << p << " bootstrap variation over orig ens:  cen=" << mu_boot_origens_avg << " +- " << mu_boot_origens_stddev
		<< " err=" << se_boot_origens_avg << " +- " << se_boot_origens_stddev << std::endl;

      //Generate histograms (only for first bootstrap original ensemble)
      {
	MatPlotLibScriptGenerate plot;
	struct acc{
	  int p;
	  const std::vector<parameterVector<double> > &d;
	  acc(const std::vector<parameterVector<double> > &d, int p): d(d), p(p){}
	  double y(const int i) const{ return d[i][p]; }
	  int size() const{ return d.size(); }
	};
	typename MatPlotLibScriptGenerate::kwargsType kwargs;
	kwargs["density"] = true;
	kwargs["alpha"] = 0.4;
	kwargs["bins"] = 60;
	auto htrue = plot.histogram(acc(fitval_dist_true,p),kwargs,"true");
	plot.setLegend(htrue, R"(${\\rm true}$)");
      
	kwargs["color"] = 'c';
	auto hboot = plot.histogram(acc(fitval_dist_boot[0],p),kwargs,"boot");
	plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

	plot.setXlabel(R"($a$)");
	plot.setYlabel(R"(${\cal F}(a)$)");
	plot.createLegend();

	std::string stub = "fitval_p" + std::to_string(p) + "_dist";
	plot.write(stub+".py",stub+".pdf");
      }

    }//p
  }//run
};


std::unique_ptr<preAnalysisBase> preAnalysisFactory(preAnalysisType type){
  if(type == preAnalysisType::None){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisNone);
  }else if(type == preAnalysisType::CovMatEvals){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisCovMatEvals);
  }else if(type == preAnalysisType::CorrMatEvals){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisCorrMatEvals);
  }else if(type == preAnalysisType::StandardError){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisStandardError);
  }else{
    error_exit(std::cout << "Invalid pre-analysis type" << std::endl);
  }
}

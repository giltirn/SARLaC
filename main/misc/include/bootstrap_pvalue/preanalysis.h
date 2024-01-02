#pragma once

struct preAnalysisBase{
  virtual void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const = 0;
  virtual ~preAnalysisBase(){}
};
struct preAnalysisNone: public preAnalysisBase{
  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{};
};
struct preAnalysisCovMatEvals: public preAnalysisBase{
  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{
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
  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{
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

  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{   
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

//Fit with the autocorrelation-avoiding covariance matrix for both (block) jackknife and (non-overlapping block) bootstrap. Ignore the covgen input
struct preAnalysisFitAutoCorrAvoid: public preAnalysisBase{

  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{   
    //Generate raw data
    int nsample = args.nsample;
    int Lt = args.Lt;
    int ntest = args.ntest;
    int dof = Lt - fitfunc.Nparams();

    int norig_ens = args.norig_ens; //to get an idea of the variations in the error

    std::vector<std::vector<bootstrapDistributionD> > params_b_all;
    std::vector<std::vector<jackknifeDistributionD> > params_j_all;
    std::vector<bootstrapDistributionD> chisq_b_all;
    std::vector<jackknifeDistributionD> chisq_j_all;
    std::vector<bool> conv_j_all(norig_ens,true), conv_b_all(norig_ens,true);

    for(int orig_ens=0; orig_ens < norig_ens; orig_ens++){
      correlationFunction<double, rawDataDistributionD> data = datagen.generate(Lt,nsample);
    
      //Block bootstrap (NBB)  with autocorrelation-ignoring covariance matrix
      std::vector<std::vector<int> > rtable = nonoverlappingBlockResampleTable(RNG,nsample,args.block_size,ntest); //ntest == nboots
      int nsample_reduced = rtable[0].size(); //allow for truncation to match multiple of block size

      bootstrapInitType binit(ntest);
      bootJackknifeInitType bjinit(nsample_reduced,nsample,ntest);
      correlationFunction<double, bootstrapDistributionD> data_rb(Lt);
      correlationFunction<double, bootJackknifeDistributionD> data_rbj(Lt);
    
      for(int t=0;t<Lt;t++){
	data_rb.coord(t) = data_rbj.coord(t) = t;
	data_rb.value(t).resize(binit);
	data_rb.value(t).resample(data.value(t), rtable);
	data_rbj.value(t).resize(bjinit);
	data_rbj.value(t).resample(data.value(t), rtable);
      }
    
      //Block double-jackknife with autocorrelation-ignoring covariance matrix
      std::pair<int,int> bdjinit(nsample, args.block_size);
      int jinit = nsample / args.block_size; //block jackknife  == binned jackknife
      correlationFunction<double, jackknifeDistributionD> data_rj(Lt);
      correlationFunction<double, blockDoubleJackknifeDistributionD> data_rdj(Lt);
      for(int t=0;t<Lt;t++){
	data_rj.coord(t) = data_rdj.coord(t) = t;
	data_rj.value(t).resize(jinit);
	data_rj.value(t).resample(data.value(t).bin(args.block_size, true));
	data_rdj.value(t).resize(bdjinit);
	data_rdj.value(t).resample(data.value(t));
      }

      MarquardtLevenbergParameters<double> mlparams(args.MLparams);
      mlparams.exit_on_convergence_fail = false;

      simpleFitWrapper<bootstrapDistributionD> fit_b(fitfunc, MinimizerType::MarquardtLevenberg, mlparams);
      simpleFitWrapper<jackknifeDistributionD> fit_j(fitfunc, MinimizerType::MarquardtLevenberg, mlparams);

      fit_b.generateCovarianceMatrix(data_rbj);
      fit_j.generateCovarianceMatrix(data_rdj);

      std::vector<bootstrapDistributionD> params_b(fitfunc.Nparams(), bootstrapDistributionD(binit));
      std::vector<jackknifeDistributionD> params_j(fitfunc.Nparams(), jackknifeDistributionD(jinit));
    
      int d; 
      bootstrapDistributionD chisq_b(binit), chisq_per_dof_b(binit);
      jackknifeDistributionD chisq_j(jinit), chisq_per_dof_j(jinit);
      std::cout << "Performing bootstrap fit" << std::endl << std::flush;
      bool conv_b = fit_b.fit(params_b, chisq_b, chisq_per_dof_b, d, data_rb);
      std::cout << "Performing jackknife fit" << std::endl << std::flush;
      bool conv_j = fit_j.fit(params_j, chisq_j, chisq_per_dof_j, d, data_rj);
      if(!conv_b){
	std::cout << "Bootstrap fit did not converge" << std::endl;
	for(int i=0;i<fitfunc.Nparams();i++) params_b[i].zero();
	conv_b_all[orig_ens] = false;
      }
      if(!conv_j){
	std::cout << "Jackknife fit did not converge" << std::endl;
	for(int i=0;i<fitfunc.Nparams();i++) params_j[i].zero();
	conv_j_all[orig_ens] = false;
      }
      params_b_all.push_back(std::move(params_b));
      params_j_all.push_back(std::move(params_j));
      chisq_b_all.push_back(std::move(chisq_b));
      chisq_j_all.push_back(std::move(chisq_j));
    }

    //Compute the true error
    //(use bootstrap containers so that the error is the standard deviation)
    bootstrapInitType binit_true(ntest);
    std::vector<bootstrapDistributionD> params_true(fitfunc.Nparams(), bootstrapDistributionD(binit_true));
    bootstrapDistributionD chisq_true(binit_true);

#pragma omp parallel for
    for(int orig_ens=0; orig_ens < ntest; orig_ens++){
      correlationFunction<double, rawDataDistributionD> data = datagen.generate(Lt,nsample);
      correlationFunction<double, double> data_means(Lt);
      for(int t=0;t<Lt;t++){
	data_means.coord(t) = t;
	data_means.value(t) = data.value(t).mean();
      }
      NumericSquareMatrix<double> cov(Lt);
      for(int t=0;t<Lt;t++)
	for(int u=t;u<Lt;u++)
	  cov(t,u) = cov(u,t) = rawDataDistributionD::covariance(data.value(t),data.value(u));
	  
      simpleSingleFitWrapper fitter(fitfunc, MinimizerType::MarquardtLevenberg, args.MLparams);
      fitter.importCovarianceMatrix(cov);
      
      parameterVector<double> params(fitfunc.Nparams(),0.);
      double q2, q2_per_dof; int dof;
      assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));

      for(int p=0;p<fitfunc.Nparams();p++) params_true[p].sample(orig_ens) = params[p];
      chisq_true.sample(orig_ens) = q2;
    }

    std::cout << "Fit results,  jackknife : bootstrap" << std::endl;
    for(int o=0;o<norig_ens;o++){
      std::cout << "Orig ens " << o << std::endl;
      for(int i=0;i<fitfunc.Nparams();i++){
	std::cout << i << " " << params_j_all[o][i] << " : " << params_b_all[o][i] << std::endl;
      }
    }

    writeParamsStandard(params_b_all,"fit_params_b.hdf5");
    writeParamsStandard(params_j_all,"fit_params_j.hdf5");
    writeParamsStandard(chisq_b_all,"chisq_b.hdf5");
    writeParamsStandard(chisq_j_all,"chisq_j.hdf5");
    
    //Compute the average error and how much it varies over original ensembles
    std::ofstream err_b("err_b.dat"), err_j("err_j.dat");    
    for(int i=0;i<fitfunc.Nparams();i++){
      rawDataDistributionD b, j;
      for(int o=0;o<norig_ens;o++){
	if(conv_b_all[o]) b.push_back(params_b_all[o][i].standardError());
	if(conv_j_all[o]) j.push_back(params_j_all[o][i].standardError());
      }
      //Only count those that converged in average
      if(b.size()) err_b << i << " " << b.mean() << " " << b.standardDeviation() << " " << b.size() << std::endl;
      else err_b << i << " NA NA 0" << std::endl;

      if(j.size()) err_j << i << " " << j.mean() << " " << j.standardDeviation() << " " << j.size() << std::endl;
      else err_j << i << " NA NA 0" << std::endl;
    }

    //True error
    std::cout << "Fit results,  true" << std::endl;
    for(int i=0;i<fitfunc.Nparams();i++){
      params_true[i].best() = params_true[i].mean();
      std::cout << i << " " << params_true[i] << std::endl;
    }
    writeParamsStandard(params_true,"fit_params_true.hdf5");
    chisq_true.best() = chisq_true.mean();
    writeParamsStandard(chisq_true,"chisq_true.hdf5");

    std::ofstream err_true("err_true.dat");
    for(int i=0;i<fitfunc.Nparams();i++){
      err_true << i << " " << params_true[i].standardError() << std::endl;
    }
  }//run
};

#define ARGS (int, delta_max)(int, bin_size)(int, nsample)
struct preAnalysisTauIntArgs{
  GENERATE_MEMBERS(ARGS);
  preAnalysisTauIntArgs(): delta_max(10), bin_size(1), nsample(50000){}
};
GENERATE_PARSER( preAnalysisTauIntArgs, ARGS );
#undef ARGS

struct preAnalysisTauInt: public preAnalysisBase{

  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{   
    std::cout << "Analysing integrated autocorrelation time for data generator" << std::endl;
    preAnalysisTauIntArgs pargs;
    if(params_file == "TEMPLATE"){
      std::ofstream of("preanalysis_tau_int_template.args");
      of << pargs;
      return;
    }
    parse(pargs, params_file);

    auto data = datagen.generate(1, pargs.nsample);
    
    std::cout << "Computing error bars using bin/resample of products ( v[s] - <v> )( v[s + delta] - <v> ) with bin size " << pargs.bin_size << std::endl;
    int nbin = (pargs.nsample - pargs.delta_max)/pargs.bin_size;
    std::vector<std::vector<int> > rtable = resampleTable(RNG, nbin);
    AutoCorrelationOptions opt;
    auto tau_int = integratedAutocorrelationMulti(pargs.delta_max, pargs.bin_size, pargs.delta_max, data.value(0), rtable, opt);

    for(int d=0;d<=pargs.delta_max;d++) std::cout << d << " " << tau_int[d] << std::endl;

    std::cout << "Generating plot" << std::endl;
    MatPlotLibScriptGenerate plot;
    
    struct Accessor{
      const std::vector<bootstrapDistributionD> &vec;
      Accessor(const std::vector<bootstrapDistributionD> &vec): vec(vec){}
      
      inline double x(const int i) const{ return i; }
      inline double upper(const int i) const{ return vec[i].confidenceRegion().second; }
      inline double lower(const int i) const{ return vec[i].confidenceRegion().first; }
      inline int size() const{ return vec.size(); }
    };
    
    Accessor acc(tau_int);
    plot.errorBand(acc);
    plot.setXlabel(R"($\Delta_{\rm cut}$)");
    plot.setYlabel(R"($\tau_{\rm int}(\Delta_{\rm cut}$))");

    plot.write("plot_tau_int.py","plot_tau_int.pdf");
  }
};


//Analyze for bias in the means of block resampled ensembles
struct preAnalysisBlockBootstrapMeanBias : public preAnalysisBase{

  void run(const Args &args, const std::string &params_file, const covMatStrategyBase &covgen, const randomDataBase &datagen, genericFitFuncBase &fitfunc) const override{   
    //Generate raw data
    int nsample = args.nsample;
    int ntest = args.ntest;

    //Do many times and plot a histogram of the bias
    struct acc{
      const std::vector<double> &d;
      acc(const std::vector<double> &d): d(d){}
      double y(const int i) const{ return d[i]; }
      int size() const{ return d.size(); }
    };

    int norig_ens = args.norig_ens;
    std::vector<double> bias(norig_ens);
    double bias_avg = 0.;
    for(int e=0;e<norig_ens;e++){
      std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
      correlationFunction<double, rawDataDistributionD> orig_data = datagen.generate(1,nsample);
      double orig_mean = orig_data.value(0).mean();
      double rmeans_avg = 0.;
      std::vector<double> rmeans(ntest);
      for(int test=0;test<ntest;test++){
	rawDataDistributionD rdata = resampledEnsemble(orig_data.value(0),test,rtable);
	rmeans[test] = rdata.mean();
	rmeans_avg += rmeans[test];
      }
      rmeans_avg /= ntest;
      bias[e] = rmeans_avg - orig_mean;
      bias_avg += bias[e];
    }
    bias_avg /= norig_ens;

    MatPlotLibScriptGenerate plot_bias;
    typename MatPlotLibScriptGenerate::kwargsType kwargs;
    kwargs["density"] = true;
    kwargs["bins"] = 60;

    plot_bias.histogram(acc(bias),kwargs);
    kwargs.clear();
    kwargs["color"] = 'k';
    plot_bias.verticalLine(0.,kwargs);
    kwargs["color"] = 'b';
    plot_bias.verticalLine(bias_avg,kwargs);
    
    std::string fname = "bootstrap_means_bias";
    plot_bias.write(fname + ".py",fname + ".pdf");
  }
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
  }else if(type == preAnalysisType::FitAutoCorrAvoid){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisFitAutoCorrAvoid);
  }else if(type == preAnalysisType::TauInt){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisTauInt);
  }else if(type == preAnalysisType::BlockBootstrapMeanBias){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisBlockBootstrapMeanBias);
  }else{
    error_exit(std::cout << "Invalid pre-analysis type" << std::endl);
  }
}

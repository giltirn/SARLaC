#pragma once

struct preAnalysisBase{
  virtual void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const = 0;
  virtual ~preAnalysisBase(){}
};
struct preAnalysisNone: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{};
};
struct preAnalysisCovMatEvals: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{
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



  }
};

struct preAnalysisCorrMatEvals: public preAnalysisBase{
  void run(const Args &args, const covMatStrategyBase &covgen, const randomDataBase &datagen) const override{
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



std::unique_ptr<preAnalysisBase> preAnalysisFactory(preAnalysisType type){
  if(type == preAnalysisType::None){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisNone);
  }else if(type == preAnalysisType::CovMatEvals){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisCovMatEvals);
  }else if(type == preAnalysisType::CorrMatEvals){
    return std::unique_ptr<preAnalysisBase>(new preAnalysisCorrMatEvals);
  }else{
    error_exit(std::cout << "Invalid pre-analysis type" << std::endl);
  }
}

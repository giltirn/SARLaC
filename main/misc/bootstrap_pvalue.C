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

struct Model{
  const genericFitFuncBase &ffunc;
  parameterVector<double> params;
  
  double value(const int t) const{ return ffunc.value(generalContainer(double(t)), params); }

  Model(const genericFitFuncBase &ffunc, const parameterVector<double> &params): ffunc(ffunc), params(params){}
};

//For data drawn from the underlying distribution, center the data on the model
void centerEnsembleOnModel(correlationFunction<double, rawDataDistributionD> &data, const Model &model, const randomDataBase &dgen){
  int nsample = data.value(0).size();
  std::vector<double> true_means = dgen.populationTimesliceMeans();
  for(int t=0;t<data.size();t++) data.value(t) = data.value(t) + rawDataDistributionD(nsample, model.value(t)-true_means[t]);
}

struct bootstrapAnalyzeOpts{
  std::string write_rtable; //write the rtable if non-empty string
  bool write_rtable_include_q2val; //include the corresponding q^2 value as the first entry on each line of the rtable
  bool write_rtable_sort_q2val; //sort the rows of the rtable over the rows by q^2 in ascending order
  bool write_rtable_sort_samples; //sort the columns of the rtable in ascending order, making it easier to see repeats
  std::vector<correlationFunction<double, double> > *data_means_out; //output the data means if non-null

  bootstrapAnalyzeOpts(): write_rtable(""), data_means_out(nullptr), write_rtable_include_q2val(false), write_rtable_sort_q2val(false), write_rtable_sort_samples(false){}
};

//data_means_out :  output unrecentered data means
void bootstrapAnalyze(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Model &model, const Args &args, const bootstrapAnalyzeOpts &opts = bootstrapAnalyzeOpts()){
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

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
  assert(rtable.size() == ntest && rtable[0].size() == nsample_reduced);
   
  if(opts.data_means_out) opts.data_means_out->resize(ntest);

#pragma omp parallel for
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data(Lt);
    correlationFunction<double, double> data_means(Lt), data_means_unrecentered(Lt);
    for(int t=0;t<Lt;t++){
      double fit_value = model.value(t);
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

    if(opts.data_means_out) (*opts.data_means_out)[test] = std::move(data_means_unrecentered);
  }

  if(opts.write_rtable.size()){
    std::vector<int> tests_sorted(ntest); 
    if(opts.write_rtable_sort_q2val){
      for(int e=0;e<ntest;e++) tests_sorted[e] = e;
      std::sort(tests_sorted.begin(), tests_sorted.end(), [&](const int a, const int b){ return q2_into[a] < q2_into[b]; });
    }

    std::ofstream f(opts.write_rtable);
    for(int e=0;e<ntest;e++){
      int ee = opts.write_rtable_sort_q2val ? tests_sorted[e] : e;

      if(opts.write_rtable_include_q2val)
	f << q2_into[ee] << " ";

      if(args.bootstrap_strat == BootResampleTableType::NonOverlappingBlock){
	//output block indices for easier reading!
	std::vector<int> bidx(nblock);
	for(int s=0;s<nblock;s++)
	  bidx[s] = rtable[ee][s*args.block_size]/args.block_size;

	if(opts.write_rtable_sort_samples) std::sort(bidx.begin(),bidx.end());

	for(int s=0;s<nblock;s++)
	  f << bidx[s] << " ";

      }else{
	if(opts.write_rtable_sort_samples)  std::sort(rtable[ee].begin(),rtable[ee].end());

	for(int s=0;s<rtable[ee].size();s++)
	  f << rtable[ee][s] << " ";
      }
 
      f << std::endl;
    }
  }
} 

void bootstrapAnalyzeResiduals(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Model &model, const Args &args, const std::string &write_rtable = ""){
  std::cout << "Bootstrap analysis using residuals" << std::endl;
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

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
  assert(rtable.size() == ntest && rtable[0].size() == nsample_reduced);
   
  correlationFunction<double, rawDataDistributionD> orig_data_resids(Lt);
  for(int t=0;t<Lt;t++){
    orig_data_resids.coord(t) = t;
    orig_data_resids.value(t) = orig_data.value(t) - rawDataDistributionD(nsample, model.value(t));
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
      for(int s=0;s<rtable[e].size();s++)
	f << rtable[e][s] << " ";
      f << std::endl;
    }
  }
} 

void bootstrapAnalyzeResidualsDiag(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Model &model, const Args &args, const std::string &write_rtable = ""){
  std::cout << "Bootstrap analysis using residuals *with diagonalization*" << std::endl;
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
  NumericSquareMatrix<double> orig_data_cov;
  {
    simpleSingleFitWrapper fitter(ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen.compute(fitter, orig_data);
    orig_data_cov = fitter.getCovarianceMatrix();
  }    
  std::vector<NumericVector<double> > orig_cov_evecs(Lt, NumericVector<double>(Lt));
  std::vector<double> orig_cov_evals(Lt);
  symmetricMatrixEigensolve(orig_cov_evecs, orig_cov_evals, orig_data_cov);

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
  assert(rtable.size() == ntest && rtable[0].size() == nsample_reduced);
   
  correlationFunction<double, rawDataDistributionD> orig_data_resids_unrot(Lt);
  for(int t=0;t<Lt;t++){
    orig_data_resids_unrot.coord(t) = t;
    double fitval = model.value(t);
    orig_data_resids_unrot.value(t) = orig_data.value(t) - rawDataDistributionD(nsample, fitval);
  }  
  correlationFunction<double, rawDataDistributionD> orig_data_resids(Lt);
  for(int i=0;i<Lt;i++){
    orig_data_resids.coord(i) = i;
    orig_data_resids.value(i) = rawDataDistributionD(nsample, 0.);
    for(int t=0;t<Lt;t++)
      orig_data_resids.value(i) = orig_data_resids.value(i) + orig_cov_evecs[i](t)*orig_data_resids_unrot.value(t);
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
      for(int s=0;s<rtable[e].size();s++)
	f << rtable[e][s] << " ";
      f << std::endl;
    }
  }
} 


void bootstrapAnalyzeResidualsDiagEvals(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Model &model, const Args &args, const std::string &write_rtable = ""){
  std::cout << "Bootstrap analysis using residuals *with diagonalization*" << std::endl;
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
  NumericSquareMatrix<double> orig_data_cov;
  {
    simpleSingleFitWrapper fitter(ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen.compute(fitter, orig_data);
    orig_data_cov = fitter.getCovarianceMatrix();
  }    
  std::vector<NumericVector<double> > orig_cov_evecs(Lt, NumericVector<double>(Lt));
  std::vector<double> orig_cov_evals(Lt);
  symmetricMatrixEigensolve(orig_cov_evecs, orig_cov_evals, orig_data_cov);

  std::vector<std::vector<int> > rtable = generateResampleTable(nsample, ntest, args.bootstrap_strat, args.block_size, threadRNG);
  assert(rtable.size() == ntest && rtable[0].size() == nsample_reduced);
   
  correlationFunction<double, rawDataDistributionD> orig_data_resids_unrot(Lt);
  for(int t=0;t<Lt;t++){
    orig_data_resids_unrot.coord(t) = t;
    double fitval = model.value(t);
    orig_data_resids_unrot.value(t) = orig_data.value(t) - rawDataDistributionD(nsample, fitval);
  }  
  correlationFunction<double, rawDataDistributionD> orig_data_resids(Lt);
  for(int i=0;i<Lt;i++){
    orig_data_resids.coord(i) = i;
    orig_data_resids.value(i) = rawDataDistributionD(nsample, 0.);
    for(int t=0;t<Lt;t++){
      orig_data_resids.value(i) = orig_data_resids.value(i) + orig_cov_evecs[i](t)*orig_data_resids_unrot.value(t);
    }
    orig_data_resids.value(i) = orig_data_resids.value(i) * (1./sqrt(orig_cov_evals[i])); //normalize such that q^2 = \sum_i |r|^2
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
    double &q2 = q2_into[test]; 
    q2=0.;
    for(int t=0;t<Lt;t++)
      q2 += pow(data_means.value(t),2);
  }

  if(write_rtable.size()){
    std::ofstream f(write_rtable);
    for(int e=0;e<ntest;e++){
      for(int s=0;s<rtable[e].size();s++)
	f << rtable[e][s] << " ";
      f << std::endl;
    }
  }
} 




void independentEnsAnalyze(std::vector<double> &q2_into, const correlationFunction<double, rawDataDistributionD> &orig_data, randomDataBase &dgen, const covMatStrategyBase &covgen, const genericFitFuncBase &ffunc, const Model &model, const Args &args){
  int nsample = args.nsample;
  int Lt = args.Lt;
  int ntest = args.ntest;
  assert(q2_into.size() == ntest);

#pragma omp parallel for
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data = dgen.generate(Lt, nsample);(Lt);
    centerEnsembleOnModel(data,model,dgen);
    correlationFunction<double, double> data_means(Lt);
    for(int t=0;t<Lt;t++){
      data_means.coord(t) = t;
      data_means.value(t) = data.value(t).mean();
    }
    simpleSingleFitWrapper fitter(ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen.compute(fitter, data);

    parameterVector<double> params(ffunc.Nparams(),0.);
    double q2, q2_per_dof; int dof;
    assert(fitter.fit(params,q2,q2_per_dof,dof,data_means));
    q2_into[test] = q2;
  }
} 



int main(const int argc, const char** argv){
  CMDline cmdline(argc,argv,3);

  RNG.initialize(cmdline.seed);
  threadRNG.initialize(cmdline.seed_thr);

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
  int dof = Lt - ffunc->Nparams();

  //Run any pre-analysis
  if(args.preanalysis.size() != args.preanalysis_params_file.size()) error_exit(std::cout << "Require as many preanalysis params files as there are preanalyses specified" << std::endl);

  for(int i=0;i<args.preanalysis.size();i++){
    std::unique_ptr<preAnalysisBase> preanalysis = preAnalysisFactory(args.preanalysis[i]);
    preanalysis->run(args, args.preanalysis_params_file[i], *covgen, *dgen, *ffunc);
  }
  if(cmdline.exit_after_preanalysis) return 0;

  //------------------------------------------------------------------------------------------------------------------------------------------
  //Generate the base original ensemble from which we obtain the model that we are interested in computing the null distribution for
  //------------------------------------------------------------------------------------------------------------------------------------------
  correlationFunction<double, rawDataDistributionD> base_orig_ens = dgen->generate(Lt, nsample);
  parameterVector<double> base_fit_params(ffunc->Nparams(),0.);
  double base_q2; //q^2 from the base original ensemble
  {
    correlationFunction<double,double> data_means(Lt);
    for(int t=0;t<Lt;t++){
      data_means.coord(t) = t;
      data_means.value(t) = base_orig_ens.value(t).mean();
    }
    simpleSingleFitWrapper fitter(*ffunc, MinimizerType::MarquardtLevenberg, args.MLparams);
    covgen->compute(fitter, base_orig_ens);
    double q2_per_dof; int dof;
    assert(fitter.fit(base_fit_params,base_q2,q2_per_dof,dof, data_means));
    std::cout << "Base fit q^2=" << base_q2 << std::endl;
  }      
  Model the_model(*ffunc,base_fit_params);
  
  //------------------------------------------------------------------------------------------------------------------------------------------
  //Generate the true distribution. Data must be centered on the specific model not the true population center
  //------------------------------------------------------------------------------------------------------------------------------------------
  std::vector<double> q2_dist_true(ntest);
  std::vector<double> mean_dist_true(Lt*ntest);

#pragma omp parallel for schedule(static)
  for(int test=0;test<ntest;test++){
    correlationFunction<double, rawDataDistributionD> data = dgen->generate(Lt, nsample);
    if(cmdline.recenter_orig_ens) centerEnsembleOnModel(data,the_model,*dgen);
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
  }      
  
  //Compute the q^2 distribution using the bootstrap for the base original ensemble with recentering around the model
  std::vector<double> q2_dist_boot(ntest);
  bootstrapAnalyzeOpts bopts;

  if(cmdline.output_bootstrap_q2sorted_rtable){
    bopts.write_rtable = "boot_rtable_q2sorted.dat";
    bopts.write_rtable_include_q2val = true;
    bopts.write_rtable_sort_q2val = true;
    bopts.write_rtable_sort_samples = true;
  }
  
  bootstrapAnalyze(q2_dist_boot, base_orig_ens, *covgen, *ffunc, the_model, args, bopts);
  //independentEnsAnalyze(q2_dist_boot, base_orig_ens, *dgen, *covgen, *ffunc, the_model, args);

  //------------------------------------------------------------------------------------------------------------------------------
  //Repeat with bootstrap for norig_ens separate original ensembles centered on the model
  //------------------------------------------------------------------------------------------------------------------------------
  std::vector< std::vector<double> > q2_dist_boot_var(args.norig_ens, std::vector<double>(ntest));
  for(int o=0;o<args.norig_ens;o++){
    correlationFunction<double, rawDataDistributionD> data = dgen->generate(Lt,nsample);
    if(cmdline.recenter_orig_ens) centerEnsembleOnModel(data,the_model,*dgen);
    bootstrapAnalyze(q2_dist_boot_var[o], data, *covgen, *ffunc, the_model, args);
    //independentEnsAnalyze(q2_dist_boot_var[o], data, *dgen, *covgen, *ffunc, the_model, args);
  }

  //---------------------------------------------------------------
  //Repeat with bootstrap for norig_ens bootstrap resampled ensembles in place of real original ensembles for error estimation
  //Still should center them on the model
  //---------------------------------------------------------------
  std::vector< std::vector<double> > q2_dist_dbl_boot(args.norig_ens, std::vector<double>(ntest));  

  std::vector<std::vector<int> > dbl_boot_rtable = generateResampleTable(nsample, args.norig_ens, args.bootstrap_strat, args.block_size, threadRNG);
  int nblock = nsample / args.block_size;
  int nsample_reduced = nblock * args.block_size;

  assert(dbl_boot_rtable.size() == args.norig_ens && dbl_boot_rtable[0].size() == nsample_reduced);

  for(int o=0;o<args.norig_ens;o++){
    correlationFunction<double, rawDataDistributionD> data(Lt);
    for(int t=0;t<Lt;t++){
      data.coord(t) = t;
      rawDataDistributionD &dd = data.value(t);
      dd.resize(nsample_reduced);
      double shift = the_model.value(t) - base_orig_ens.value(t).mean();
      for(int s=0;s<nsample_reduced;s++)
	dd.sample(s) = base_orig_ens.value(t).sample(dbl_boot_rtable[o][s]) + shift;
    }
    bootstrapAnalyze(q2_dist_dbl_boot[o], data, *covgen, *ffunc, the_model, args);
  }


  //Compute the q^2 distribution using the bootstrap for the base original ensemble with recentering around the model
  std::vector<double> q2_resid_dist_boot(ntest);
  if(cmdline.bootstrap_resid_diagonalize) bootstrapAnalyzeResidualsDiag(q2_resid_dist_boot, base_orig_ens, *covgen, *ffunc, the_model, args);
  else if(cmdline.bootstrap_resid_diagonalize_evals) bootstrapAnalyzeResidualsDiagEvals(q2_resid_dist_boot, base_orig_ens, *covgen, *ffunc, the_model, args);
  else bootstrapAnalyzeResiduals(q2_resid_dist_boot, base_orig_ens, *covgen, *ffunc, the_model, args);

  //------------------------------------------------------------------------------------------------------------------------------
  //Repeat with bootstrap for norig_ens separate original ensembles centered on the model
  //------------------------------------------------------------------------------------------------------------------------------
  std::vector< std::vector<double> > q2_resid_dist_boot_var(args.norig_ens, std::vector<double>(ntest));
  for(int o=0;o<args.norig_ens;o++){
    correlationFunction<double, rawDataDistributionD> data = dgen->generate(Lt,nsample);
    if(cmdline.recenter_orig_ens) centerEnsembleOnModel(data,the_model,*dgen);

    if(cmdline.bootstrap_resid_diagonalize) bootstrapAnalyzeResidualsDiag(q2_resid_dist_boot_var[o], data, *covgen, *ffunc, the_model, args);
    else if(cmdline.bootstrap_resid_diagonalize_evals) bootstrapAnalyzeResidualsDiagEvals(q2_resid_dist_boot_var[o], data, *covgen, *ffunc, the_model, args);
    else bootstrapAnalyzeResiduals(q2_resid_dist_boot_var[o], data, *covgen, *ffunc, the_model, args);
  }

  //----------------------------------------------
  //Plot results
  //----------------------------------------------
  setGSLerrorHandlerThrow();

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

    kwargs["color"] = "tab:gray";
    auto hboot_resid = plot.histogram(acc(q2_resid_dist_boot),kwargs,"boot_resid");
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
      double y(const int i) const override{
	double v = 0.;
	try{
	  v= TsquareDistribution::PDF(x(i),p,n);
	}catch(const std::exception &e){
	  std::cout << "T2data caught error: " << e.what() << std::endl;
	}
	return v;
      }
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
      double y(const int i) const override{
	double v=0;
	try{
	  v= chiSquareDistribution::PDF(p,x(i)); 
	}catch(const std::exception &e){
	  std::cout << "ChisqData caught error: " << e.what() << std::endl;
	}
	return v;
      }
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
    std::vector<double> pboot(npt); //bootstrap p-value from the base original ensemble
    std::vector<double> pboot_resid(npt); //bootstrap p-value by residuals from the base original ensemble

    std::vector<rawDataDistributionD> pboot_var(npt, rawDataDistributionD(args.norig_ens) ); //variation over actual original ensembles
    std::vector<rawDataDistributionD> pdbl_boot_var(npt, rawDataDistributionD(args.norig_ens) ); //variation over bootstrap resampled ensembles in place of original ensembles
    std::vector<rawDataDistributionD> pboot_resid_var(npt, rawDataDistributionD(args.norig_ens) ); //variation over actual original ensembles for residuals version
    std::vector<double> pT2(npt);
    std::vector<double> pchi2(npt);
    
    double dq2 = q2_max/(npt-1);

    for(int i=0;i<npt;i++){
      double q2 = dq2*i;
      q2vals[i] = q2;

      try{
	pchi2[i] = chiSquareDistribution::pvalue(dof, q2);
      }catch(const std::exception &e){
	std::cout << "WARNING chi2 caught error: " << e.what() << std::endl;
	pchi2[i] = 0;
      }
      try{	
	pT2[i] = TsquareDistribution::pvalue(q2, dof, nsample-1);
      }catch(const std::exception &e){
	std::cout << "WARNING T2 caught error: " << e.what() << std::endl;
	pT2[i] = 0;
      }

      ptrue[i] = estimatePvalue(q2, q2_dist_true);
      pboot[i] = estimatePvalue(q2, q2_dist_boot);
      pboot_resid[i] = estimatePvalue(q2, q2_resid_dist_boot);
      
      for(int o=0;o<args.norig_ens;o++){
	pboot_var[i].sample(o) = estimatePvalue(q2, q2_dist_boot_var[o]);
	pdbl_boot_var[i].sample(o) = estimatePvalue(q2, q2_dist_dbl_boot[o]);
	pboot_resid_var[i].sample(o) = estimatePvalue(q2, q2_resid_dist_boot_var[o]);
      }
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
      auto hboot = plot.errorBand(acc_wsep_err(q2vals,pboot,pboot_var), kwb, "boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      kwb["color"] = "tab:pink";
      auto hdbl_boot = plot.errorBand(acc_wsep_err(q2vals,pboot,pdbl_boot_var), kwb, "dbl_boot");
      plot.setLegend(hdbl_boot, R"(${\\rm dbl. bootstrap}$)");

      kwb["color"] = "tab:gray";
      auto hboot_resid = plot.errorBand(acc_wsep_err(q2vals,pboot_resid,pboot_resid_var), kwb, "boot_resid");
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
      auto hboot = plot.errorBand(acc_wsep_err(ptrue,pboot,pboot_var), kwb, "boot");
      plot.setLegend(hboot, R"(${\\rm bootstrap}$)");

      kwb["color"] = "tab:pink";
      auto hdbl_boot = plot.errorBand(acc_wsep_err(ptrue,pboot,pdbl_boot_var), kwb, "dbl_boot");
      plot.setLegend(hdbl_boot, R"(${\\rm dbl. bootstrap}$)");

      kwb["color"] = "tab:gray";
      auto hboot_resid = plot.errorBand(acc_wsep_err(ptrue,pboot_resid,pboot_resid_var), kwb, "boot_resid");
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

  //Compute p-value of original fit
  double pT2 = TsquareDistribution::PDF(base_q2,dof,nsample-1);
  double pchi2 = chiSquareDistribution::PDF(dof,base_q2);
  double ptrue = estimatePvalue(base_q2, q2_dist_true);
  double pboot = estimatePvalue(base_q2, q2_dist_boot);
  double pboot_resid = estimatePvalue(base_q2, q2_resid_dist_boot);

  std::cout << "Computing p-values for initial fit, q^2=" << base_q2 << ":" << std::endl;
  std::cout << "T^2: " << pT2 << std::endl;
  std::cout << "chi^2: " << pchi2 << std::endl;
  std::cout << "true: " << ptrue << std::endl;
  std::cout << "boot: " << pboot << std::endl;
  std::cout << "boot-resid: " << pboot_resid << std::endl;

  std::cout << "Done\n";
  return 0;
}


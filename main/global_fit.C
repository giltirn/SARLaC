#include<fit.h>
#include<distribution.h>
#include<data_series.h>
#include<utils.h>
#include<containers.h>
#include<tensors.h>
#include<random.h>

using namespace CPSfit;

#include<global_fit/utils.h>
#include<global_fit/data_params.h>
#include<global_fit/read_reweighted.h>
#include<global_fit/fitfunc_general.h>
#include<global_fit/def_binding.h>
#include<global_fit/globalfit_params.h>
#include<global_fit/fitfunc_analytic.h>

template<typename ParameterType>
inline void freezeToConstant(std::vector<int> &freeze_params, ParameterType &freeze_const, const std::vector<int> &to_freeze, const double val){
  for(int i=0;i<to_freeze.size();i++){
    freeze_params.push_back(to_freeze[i]);
    freeze_const(to_freeze[i]) = val;
  }
}

void set32Ilayout(superJackknifeLayout &layout){
  layout.addEnsemble("32I_ml0.004", 300/4);
  layout.addEnsemble("32I_ml0.006", 312/4);
  layout.addEnsemble("32I_ml0.008", 252/4);
}



void read32I(DataSeriesType &data, const superJackknifeLayout &layout){
  std::vector<int> rwidx_keep = {0,2,4,6};  //this ensemble has 9 reweightings, we wish to keep 4. Original unreweighted data on index 0
  ReweightInfo rwinfo(0.03,0.026,-0.0005,0.03);

  const double NA = 0.;

  int lat = 0;

  //Load mres
  std::vector<superJackknifeDistribution<double> > mres_rw = readReweightedMresQDP("/home/ckelly/projects/32nt64_fullanalysis/chiral/mres/jackknife/reweighted/boot/mres_extrapcoords.boot",  
										   rwidx_keep, 1);
  //  These ensembles have different tags than we want, so remap them
  std::map<std::string, std::string> ens_map = {{"0.004", "32I_ml0.004"}, {"0.006", "32I_ml0.006"}, {"0.008","32I_ml0.008"}};
  for(int i=0;i<mres_rw.size();i++) mres_rw[i] = remapLayout(mres_rw[i], layout, ens_map);

  //Load ensemble data
  std::vector<double> ml_values = {0.004, 0.006, 0.008};
  std::vector<std::string> enstags = {"32I_ml0.004","32I_ml0.006","32I_ml0.008"};
  for(int mli=0;mli<ml_values.size();mli++){
    double ml = ml_values[mli];
    std::string ens = enstags[mli];

    //mpi
    std::vector<double> mx_values = {0.008, 0.006, 0.004, 0.002};
    std::vector<double> my_values = {0.008, 0.006, 0.004, 0.002};
    for(int mxi=0;mxi<mx_values.size();mxi++){
      double mx = mx_values[mxi];
      for(int myi=0;myi<my_values.size();myi++){
	double my = my_values[myi];

	if(my < mx) continue;
	  
	std::string file = stringize("/home/ckelly/projects/32nt64_fullanalysis/m%s/fpi/jackknife/reweighted/boot/fit_%s_%s_12_52.boot",dstr(ml).c_str(),dstr(my).c_str(),dstr(mx).c_str());
	DataParams base(lat,mpi2,mx,my,ml,NA);
	readAddReweighedDataQDP(data,file,0,base,ens,rwinfo,rwidx_keep,mres_rw, [&](jackknifeDistribution<double> &v){ v = v*v; });
      }
    }

    //mK (nb my is the heavier quark)
    mx_values = {0.008, 0.006, 0.004, 0.002};
    my_values = {0.03, 0.025};
    for(int mxi=0;mxi<mx_values.size();mxi++){
      double mx = mx_values[mxi];
      for(int myi=0;myi<my_values.size();myi++){
	double my = my_values[myi];
	  
	std::string file = stringize("/home/ckelly/projects/32nt64_fullanalysis/m%s/fpi/jackknife/reweighted/boot/fit_%s_%s_12_52.boot",dstr(ml).c_str(),dstr(my).c_str(),dstr(mx).c_str());
	DataParams base(lat,mK2,mx,my,ml,NA);
	readAddReweighedDataQDP(data,file,0,base,ens,rwinfo,rwidx_keep,mres_rw, [&](jackknifeDistribution<double> &v){ v = v*v; });
      }
    }

    //mOmega
    {
      std::string file = stringize("/home/ckelly/projects/32nt64_fullanalysis/m%s/momega/jackknife/reweighted/boot/fit_0.03_7_13.boot",dstr(ml).c_str());
      DataParams base(lat,mOmega,0.0,0.03,ml,NA);
      readAddReweighedDataQDP(data,file,1,base,ens,rwinfo,rwidx_keep,mres_rw);
    }

  }
}



void set24Ilayout(superJackknifeLayout &layout){
  layout.addEnsemble("24I_ml0.005", 202/2);
  layout.addEnsemble("24I_ml0.01", 178/2);
}



void read24I(DataSeriesType &data, const superJackknifeLayout &layout){
  std::vector<int> rwidx_keep = {0,9,18,27}; 
  ReweightInfo rwinfo(0.04,0.033,-0.00025,0.04);

  const double NA = 0.;

  int lat = 1;

  //Load mres
  std::vector<superJackknifeDistribution<double> > mres_rw = readReweightedMresQDP("/home/ckelly/projects/24nt64_fullanalysis/chiral/mres/jackknife/reweighted/boot/mres_extrapcoords.boot",  
										   rwidx_keep, 1);
  //  These ensembles have different tags than we want, so remap them
  std::map<std::string, std::string> ens_map = {{"0.005", "24I_ml0.005"}, {"0.01", "24I_ml0.01"}};
  for(int i=0;i<mres_rw.size();i++) mres_rw[i] = remapLayout(mres_rw[i], layout, ens_map);


  //Load ensemble data
  std::vector<double> ml_values = {0.005, 0.01};
  std::vector<std::string> enstags = {"24I_ml0.005","24I_ml0.01"};
  for(int mli=0;mli<ml_values.size();mli++){
    double ml = ml_values[mli];
    std::string ens = enstags[mli];

    //mpi
    std::vector<double> mx_values = {0.01,0.005,0.001};
    std::vector<double> my_values = {0.01,0.005,0.001};
    for(int mxi=0;mxi<mx_values.size();mxi++){
      double mx = mx_values[mxi];
      for(int myi=0;myi<my_values.size();myi++){
	double my = my_values[myi];

	if(my < mx) continue;

	std::string file = stringize("/home/ckelly/projects/24nt64_fullanalysis/%s/fpi/jackknife/reweighted/boot/fit_%s_%s_10_50.boot",dstr(ml).c_str(),dstr(my).c_str(),dstr(mx).c_str());
	DataParams base(lat,mpi2,mx,my,ml,NA);
	readAddReweighedDataQDP(data,file,0,base,ens,rwinfo,rwidx_keep,mres_rw, [&](jackknifeDistribution<double> &v){ v = v*v; });
      }
    }

    //mK (nb my is the heavier quark)
    mx_values = {0.01,0.005,0.001};
    my_values = {0.04, 0.03};
    for(int mxi=0;mxi<mx_values.size();mxi++){
      double mx = mx_values[mxi];
      for(int myi=0;myi<my_values.size();myi++){
	double my = my_values[myi];
	  
	std::string file = stringize("/home/ckelly/projects/24nt64_fullanalysis/%s/fpi/jackknife/reweighted/boot/fit_%s_%s_10_50.boot",dstr(ml).c_str(),dstr(my).c_str(),dstr(mx).c_str());
	DataParams base(lat,mK2,mx,my,ml,NA);
	readAddReweighedDataQDP(data,file,0,base,ens,rwinfo,rwidx_keep,mres_rw, [&](jackknifeDistribution<double> &v){ v = v*v; });
      }
    }
    
    //mOmega
    {
      std::string file = stringize("/home/ckelly/projects/24nt64_fullanalysis/%s/momega/jackknife/reweighted/boot/fit_0.04_5_11.boot",dstr(ml).c_str());
      DataParams base(lat,mOmega,0.0,0.04,ml,NA);
      readAddReweighedDataQDP(data,file,1,base,ens,rwinfo,rwidx_keep,mres_rw);
    }

  }
}

  

int main(void){
  RNG.initialize(1234);
  
  superJackknifeLayout layout;
  set32Ilayout(layout);
  set24Ilayout(layout);
  
  DataSeriesType data;
  read32I(data,layout);
  read24I(data,layout);


  AnalyticGlobalFitPolicies::fitFuncSetupType fparams;
  fparams.nlat = 2;
  fparams.a2_term = true;
  fparams.mh0 = 0.03;

  typedef AnalyticGlobalFit FitFunc;
  FitFunc gfit(fparams);
  
  typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, uncorrelatedFitPolicy, globalFitTypes>::type FitPolicy;

  MarquardtLevenbergParameters<double> ml_params;
  ml_params.verbose = true;
  ml_params.lambda_factor = 2;
  ml_params.lambda_start = 100;

  fitter<FitPolicy> fit(ml_params);
  
  NumericVector<superJackknifeDistribution<double> > sigma(data.size(), [&](const int i){ return superJackknifeDistribution<double>(layout, data.value(i).standardError()); } );

  fit.importCostFunctionParameters(sigma);
  fit.importFitFunc(gfit);

  typename FitFunc::ParameterType guess = gfit.guess();

  //Set up freezes
  std::vector<int> freeze_params;
  typename FitFunc::ParameterType freeze_const(guess);
  
  //Scaling parameters on 0th latt should be 1.
  freezeToConstant(freeze_params, freeze_const, {freeze_const.paramIdx(&freeze_const.Zl[0]), freeze_const.paramIdx(&freeze_const.Zh[0]), freeze_const.paramIdx(&freeze_const.Ra[0])}, 1.);

  // //Temporary: freeze a^2 scaling terms
  // freezeToConstant(freeze_params, freeze_const, {freeze_const.paramIdx(&freeze_const.params.mpi2_ca), freeze_const.paramIdx(&freeze_const.params.mK2_ca), freeze_const.paramIdx(&freeze_const.params.mOmega_ca)}, 0.);

  // //Temporary: freeze mOmega valence heavy (only one value)
  // freezeToConstant(freeze_params, freeze_const, {freeze_const.paramIdx(&freeze_const.params.mOmega_cvh)}, 0.);


  typename FitPolicy::FitParameterDistribution freeze(layout, freeze_const);

  fit.freeze(freeze_params, freeze);

  std::cout << "Freeze parameters " << freeze_params << " to " << freeze_const << std::endl;


  //Perform the fit
  typename FitPolicy::FitParameterDistribution params(layout, guess);

  superJackknifeDistribution<double> chisq(layout, 0.);
  superJackknifeDistribution<double> chisq_per_dof(layout, 0.);

  fit.fit(params, chisq, chisq_per_dof, data);

  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;
  std::cout << "Results: " << params << std::endl;
}

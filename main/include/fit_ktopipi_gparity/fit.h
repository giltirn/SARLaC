#ifndef _KTOPIPI_FIT_H_
#define _KTOPIPI_FIT_H_

//Fit to each Q independently
template<typename FitFunc, template<typename> class corrUncorrFitPolicy>
struct fit_corr_uncorr{
  static std::vector<jackknifeDistribution<typename FitFunc::Params> > fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
									   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,
									   const Args &args, const CMDline &cmdline){
    const int nsample = A0_fit_j[0].value(0).size();
    typedef typename FitFunc::Params Params;
    std::vector<Params> guess(10);
    std::vector<jackknifeDistribution<Params> > fit_params(10);

    //Perform the fit
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
    FitFunc fitfunc;
  
    for(int q=0;q<10;q++){
      std::cout << "Starting fit for Q=" << q+1 << std::endl;
      fitter<FitPolicies> fitter;
      fitter.importFitFunc(fitfunc);

      //Compute the correlation matrix / weights
      importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, A0_fit_dj[q]);

      readFrozenParams(fitter, q+1, cmdline, nsample);
  
      jackknifeDistribution<Params> &params = fit_params[q];
      params = jackknifeDistribution<Params>(nsample, guess[q]);
      jackknifeDistributionD chisq;
      jackknifeDistributionD chisq_per_dof;
      fitter.fit(params, chisq, chisq_per_dof, A0_fit_j[q]);

      int nparams = params.sample(0).size();
      for(int p=0;p<nparams;p++){
	jackknifeDistributionD tmp(nsample, [&](const int s){ return params.sample(s)(p); }); 
	std::cout << params.sample(0).memberName(p) << ": " << tmp << std::endl;
      }

      std::cout << "Q" << q+1 << " Chisq: " << chisq << std::endl;
      std::cout << "Q" << q+1 << " Chisq/dof: " << chisq_per_dof << std::endl;
    }

    return fit_params;
  }
};

//Simultaneous fit in 10-basis
template<template<typename> class corrUncorrFitPolicy>
struct fit_corr_uncorr<FitKtoPiPiSim<10>, corrUncorrFitPolicy>{
  static jackknifeDistribution<FitKtoPiPiSim<10>::Params> fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
							      const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,
							      const Args &args, const CMDline &cmdline){
    const int nsample = A0_fit_j[0].value(0).size();
    typedef FitKtoPiPiSim<10> FitFunc;
    typedef FitFunc::Params Params;
    Params guess;
    
    //Gather the data into single containers
    int sz = 0;
    for(int q=0;q<10;q++) sz += A0_fit_j[q].size();
    std::vector<std::pair<int,int> > map(sz);
    int a=0;
    for(int q=0;q<10;q++)
      for(int i=0;i<A0_fit_j[q].size();i++)
	map[a++] = std::pair<int,int>(q,i);
        
    correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD> A0_fit_j_all(sz,
										    [&](const int s){
										      const amplitudeDataCoord &c = A0_fit_j[map[s].first].coord(map[s].second);
										      const jackknifeDistributionD &v = A0_fit_j[map[s].first].value(map[s].second);
										      amplitudeDataCoordSim cc(c, map[s].first);
										      return typename correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD>::ElementType(cc,v);
										    });
    correlationFunction<amplitudeDataCoordSim, doubleJackknifeA0StorageType> A0_fit_dj_all(sz,
											  [&](const int s){
											    const amplitudeDataCoord &c = A0_fit_dj[map[s].first].coord(map[s].second);
											    const doubleJackknifeA0StorageType &v = A0_fit_dj[map[s].first].value(map[s].second);
											    amplitudeDataCoordSim cc(c, map[s].first);
											    return typename correlationFunction<amplitudeDataCoordSim, doubleJackknifeA0StorageType>::ElementType(cc,v);
											  });
    //Perform the fit
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
    FitFunc fitfunc;
  
    std::cout << "Starting fit for all Q simultaneously" << std::endl;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);

    //Compute the correlation matrix / weights
    importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, A0_fit_dj_all);

    readFrozenParams(fitter, -1, cmdline, nsample);

    distributionPrint<jackknifeDistribution<Params> >::printer(new ktopipiParamsPrinter<FitFunc>);
    
    jackknifeDistribution<Params> params(nsample, guess);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, A0_fit_j_all);

    std::cout << "Params: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

    return params;
  }
};


//Simultaneous fit in chiral
template<template<typename> class corrUncorrFitPolicy>
struct fit_corr_uncorr<FitKtoPiPiSim<7>, corrUncorrFitPolicy>{
  static jackknifeDistribution<FitKtoPiPiSim<10>::Params> fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
							      const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,
							      const Args &args, const CMDline &cmdline){
    const int nsample = A0_fit_j[0].value(0).size();
    typedef FitKtoPiPiSim<7> FitFunc;
    typedef FitFunc::Params Params;
    Params guess;

    const int nfit_q = A0_fit_j[0].size();
    for(int i=1;i<10;i++) assert(A0_fit_j[i].size() == nfit_q);

    //Convert the data to the chiral basis and gather into single containers
    
    //Convert Q123 -> Q'123
    static const double Q123rot[3][3] = {  { 3    ,  2,    -1     },
					   { 2./5 , -2./5,  1./5  },
					   {-3./5,   3./5,  1./5  } };

    jackknifeDistributionD zero_j(nsample,0.);
    doubleJackknifeDistributionD zero_dj(nsample,0.);
    correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD> A0_fit_j_all;
    correlationFunction<amplitudeDataCoordSim, doubleJackknifeA0StorageType> A0_fit_dj_all;
      
    for(int d=0;d<nfit_q;d++){
      for(int i=1;i<10;i++) assert(A0_fit_j[i].coord(d) == A0_fit_j[0].coord(d)); //should all be at the same coordinate, just different q
      
      std::vector<jackknifeDistributionD> Qprime_j(7, zero_j);
      std::vector<doubleJackknifeA0StorageType> Qprime_dj(7, zero_dj);
#define Q(i,j) Q123rot[i-1][j-1]
      
#define MOj(i) Qprime_j[i-1]
#define MIj(i) A0_fit_j[i-1].value(d)
#define MOdj(i) Qprime_dj[i-1]
#define MIdj(i) A0_fit_dj[i-1].value(d)      
      
      for(int i=1;i<=3;i++){
	MOj(i) = Q(i,1)*MIj(1) + Q(i,2)*MIj(2) + Q(i,3)*MIj(3);
	MOdj(i) = Q(i,1)*MIdj(1) + Q(i,2)*MIdj(2) + Q(i,3)*MIdj(3);
      }
      for(int i=4;i<=7;i++){
	MOj(i) = MIj(i+1); //5->4  6->5 etc
	MOdj(i) = MIdj(i+1);
      }
#undef MOj
#undef MOdj
#undef MIj
#undef MIdj
#undef Q

      for(int q=0;q<7;q++){      
	amplitudeDataCoordSim cc(A0_fit_j[0].coord(d), q);	
	A0_fit_j_all.push_back(correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD>::ElementType(cc,Qprime_j[q]));
	A0_fit_dj_all.push_back(correlationFunction<amplitudeDataCoordSim, doubleJackknifeA0StorageType>::ElementType(cc,Qprime_dj[q]));
      }
    }

    //Perform the fit
    typedef typename composeFitPolicy<FitFunc, frozenFitFuncPolicy, corrUncorrFitPolicy>::type FitPolicies;
    FitFunc fitfunc;
  
    std::cout << "Starting fit for all Q simultaneously" << std::endl;
    fitter<FitPolicies> fitter;
    fitter.importFitFunc(fitfunc);

    //Compute the correlation matrix / weights
    importCostFunctionParameters<corrUncorrFitPolicy,FitPolicies> prepare(fitter, A0_fit_dj_all);

    readFrozenParams(fitter, -1, cmdline, nsample);

    distributionPrint<jackknifeDistribution<Params> >::printer(new ktopipiParamsPrinter<FitFunc>);
    
    jackknifeDistribution<Params> params(nsample, guess);
    jackknifeDistributionD chisq;
    jackknifeDistributionD chisq_per_dof;
    fitter.fit(params, chisq, chisq_per_dof, A0_fit_j_all);

    std::cout << "Params in chiral basis: " << params << std::endl;
    std::cout << "Chisq: " << chisq << std::endl;
    std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

    //Convert back to 10 basis
    jackknifeDistribution<FitKtoPiPiSim<10>::Params> out(nsample);
    
    //Convert Q'123 -> Q123
    static const double Q123invrot[3][3] = {  {1./5,   1,   0},
					      {1./5,   0,   1},
					      {  0 ,   3,   2} };
    for(int s=0;s<nsample;s++){
      out.sample(s).zero();
      for(int i=0;i<4;i++) out.sample(s)(i) = params.sample(s)(i);
      
#define MO(i) out.sample(s)(i+4-1)
#define MI(i) params.sample(s)(i+4-1)
#define Qinv(i,j) Q123invrot[i-1][j-1]

      for(int i=1;i<=3;i++)
	MO(i) = Qinv(i,1)*MI(1) + Qinv(i,2)*MI(2) + Qinv(i,3)*MI(3);

      MO(4) = MO(2) + MO(3) - MO(1); //Q4 = Q2 + Q3 - Q1    [Lehner, Sturm, arXiv:1104.4948 eq 9]
      for(int i=5;i<=8;i++) MO(i) = MI(i-1); //4->5 5->6 etc
	
      MO(9) = 3./2*MO(1)  -1./2*MO(3); //Q9 = 3/2 Q1 - 1/2 Q3
      MO(10) = 1./2*MO(1)  -1./2*MO(3) + MO(2); //Q10 = 1/2(Q1 - Q3) + Q2

#undef MO
#undef MI
#undef Qinv
    }

    std::cout << "Params in 10 basis: " << out << std::endl;
    
    return out;
  }
};





template<typename FitFunc>
inline typename fitReturnType<FitFunc>::type fit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
						 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,
						 const Args &args, const CMDline &cmdline){
  return args.correlated ?
    fit_corr_uncorr<FitFunc,correlatedFitPolicy>::fit(A0_fit_j,A0_fit_dj,args,cmdline) :
    fit_corr_uncorr<FitFunc,uncorrelatedFitPolicy>::fit(A0_fit_j,A0_fit_dj,args,cmdline);
}

void getFitData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
		std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,		
		const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_j,
		const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_dj,
		const int tmin_k_op, const int tmin_op_pi){
  //Separate out the data in the desired fit range
  for(int q=0;q<10;q++){
    for(int d=0;d<A0_j[q].size();d++){
      const int t = int(A0_j[q].coord(d).t);
      const int tsep_k_pi = A0_j[q].coord(d).tsep_k_pi;
      const int tsep_op_pi = tsep_k_pi - t;
      if(t <= tsep_k_pi && t >= tmin_k_op && tsep_op_pi >= tmin_op_pi){
	A0_fit_j[q].push_back(A0_j[q].coord(d), A0_j[q].value(d));
	A0_fit_dj[q].push_back(A0_dj[q].coord(d), A0_dj[q].value(d));
      }
    }
  }
}
inline void getFitData(std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_fit_j,
		std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_fit_dj,		
		const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_j,
		const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_dj,
		const Args &args){
  return getFitData(A0_fit_j, A0_fit_dj, A0_j, A0_dj, args.tmin_k_op, args.tmin_op_pi);
}


template<typename FitFunc>
inline void fitAndPlotFF(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
			 const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
			 const Args &args, const CMDline &cmdline){
  //Extract the data we are going to fit
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_fit_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_fit_dj(10);
  getFitData(A0_fit_j,A0_fit_dj,A0_all_j,A0_all_dj,args);
  
  std::cout << "Including " << A0_fit_j[0].size() << " data points in fit\n";
  for(int q=0;q<10;q++){
    std::cout << "For Q" << q+1 << std::endl;
    for(int i=0;i<A0_fit_j[q].size();i++)
      std::cout << A0_fit_j[q].coord(i) << " : " << A0_fit_j[q].value(i) << std::endl;
  }
  
  typename fitReturnType<FitFunc>::type fit_params = fit<FitFunc>(A0_fit_j,A0_fit_dj,args,cmdline);

  plotFF<FitFunc>::plot(A0_all_j, fit_params, args, cmdline);

#ifdef HAVE_HDF5
  writeParamsStandard(fit_params, "params.hdf5");
#endif
}

inline void fitAndPlot(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &A0_all_j,
		       const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &A0_all_dj,
		       const Args &args, const CMDline &cmdline){
  switch(args.fitfunc){
  case FitSeparate:
    return fitAndPlotFF<FitKtoPiPi>(A0_all_j,A0_all_dj,args,cmdline);
  case FitSimultaneous:
    return fitAndPlotFF<FitKtoPiPiSim<10> >(A0_all_j,A0_all_dj,args,cmdline);
  case FitSimultaneousChiralBasis:
    return fitAndPlotFF<FitKtoPiPiSim<7> >(A0_all_j,A0_all_dj,args,cmdline);
  case FitSeparateWithConstant:
    return fitAndPlotFF<FitKtoPiPiWithConstant>(A0_all_j,A0_all_dj,args,cmdline);
  case FitSeparateTwoExp:
    return fitAndPlotFF<FitKtoPiPiTwoExp>(A0_all_j,A0_all_dj,args,cmdline);
  default:
    error_exit(std::cout << "fitAndPlot(..) Unknown fit function " << args.fitfunc << std::endl);
  }
}
  

#endif

#ifndef _ANALYZE_KTOPIPI_GPARITY_READ_DATA_H_
#define _ANALYZE_KTOPIPI_GPARITY_READ_DATA_H_

//The data we need to compute the A0 and epsilon', in superJackknife format
struct superJackknifeData{
  typedef superJackknifeDistribution<double> superJackD;
  
  superJackknifeLayout layout;
  
  superJackD ainv_sj;
  superJackD omega_expt_sj;
  superJackD mod_eps_sj;
  superJackD ImA2_lat_sj;
  superJackD ReA2_lat_sj;
  superJackD ReA2_expt_sj;
  superJackD ReA0_expt_sj;
  superJackD delta_2_sj;
  
  NumericTensor<superJackD,1> M_lat_sj;
  NumericTensor<superJackD,2> NPR_sj;
  superJackD Epi_sj;
  superJackD Epipi_sj;
  superJackD mpi_sj;
  superJackD mK_sj;  

  superJackknifeData(const Args &args){
    //Load the matrix elements
    NumericTensor<jackknifeDistributionD,1> M_lat_j;
    {
      std::vector<std::vector<jackknifeDistributionD> > fit_params;
      readParamsStandard(fit_params, args.fit_results.M_lat.file);
      M_lat_j = NumericTensor<jackknifeDistributionD,1>({10}, [&](const int *i){ return fit_params[*i][args.fit_results.M_lat.idx]; } );
    }
    std::cout << "Unrenormalized lattice matrix elements:\n";
    for(int q=0;q<10;q++)
      std::cout << "Q" << q+1 << " = " << M_lat_j(&q) << std::endl;
    
    //Load the inputs needed for the Lellouch Luscher factor
    jackknifeDistributionD mK_j; readFromHDF5(mK_j,args.fit_results.mK,"mK");
    jackknifeDistributionD Epi_j; readFromHDF5(Epi_j,args.fit_results.Epi,"Epi");
    jackknifeDistributionD Epipi_j; readFromHDF5(Epipi_j,args.fit_results.Epipi,"Epipi");
    
    //Compute the pion mass  
    jackknifeDistributionD ppi2_j(Epi_j.size(),0.);
    for(int i=0;i<3;i++)
      if(args.twists[i])
	ppi2_j = ppi2_j + jackknifeDistributionD(Epi_j.size(), pow( args.lattice_dispersion_reln ? sin(M_PI/args.L) : M_PI/args.L, 2 ));  
    std::cout << "Pion p^2 = " << ppi2_j << std::endl;
    
    jackknifeDistributionD mpi_j = sqrt( Epi_j*Epi_j - ppi2_j );
    std::cout << "m_pi = " << mpi_j << std::endl;
    
    //Load the NPR matrix
    NumericTensor<jackknifeDistributionD,2> NPR_j = loadNPR(args.renormalization.file);
    std::cout << "NPR matrix:\n" << NPR_j << std::endl;
  
    //Read the various other inputs as superjackknife
    readFromXML(ainv_sj,args.other_inputs.ainv,"a^{-1}");
    readFromXML(omega_expt_sj,args.other_inputs.omega_expt,"omega_expt");
    readFromXML(mod_eps_sj,args.other_inputs.mod_eps,"mod_eps");
    readFromXML(ImA2_lat_sj,args.other_inputs.ImA2_lat,"ImA2_lat");
    readFromXML(ReA2_lat_sj,args.other_inputs.ReA2_lat,"ReA2_lat");
    readFromXML(ReA2_expt_sj,args.other_inputs.ReA2_expt,"ReA2_expt");
    readFromXML(ReA0_expt_sj,args.other_inputs.ReA0_expt,"ReA0_expt");
    readFromXML(delta_2_sj,args.other_inputs.delta_2_lat,"delta_2");

    //Create a superjackknife layout for treating all these data consistently
    std::vector<superJackD*> sjack_inputs = {&ainv_sj, &omega_expt_sj, &mod_eps_sj, &ImA2_lat_sj, &ReA2_lat_sj, &ReA2_expt_sj, &ReA0_expt_sj, &delta_2_sj};
    layout.addEnsemble("Main", M_lat_j({0}).size());
    layout.addEnsemble("NPR", NPR_j({0,0}).size() );
    for(int i=0;i<sjack_inputs.size();i++) layout = combine(layout, sjack_inputs[i]->getLayout());
    
    std::cout << "Created a superJackknife layout with size " << layout.nSamplesTotal() << " comprising " << layout.nEnsembles() << " ensembles:\n";
    for(int i=0;i<layout.nEnsembles();i++) std::cout << '\t' << layout.ensTag(i) << " of size " << layout.nSamplesEns(i) << std::endl;
  
    //Upcast everything to a superjackknife
    M_lat_sj = NumericTensor<superJackD,1>({10}, [&](const int *i){ return superJackD(layout,"Main",M_lat_j(i)); } );
    Epi_sj = superJackD(layout,"Main",Epi_j);
    Epipi_sj = superJackD(layout,"Main",Epipi_j);
    mpi_sj = superJackD(layout,"Main",mpi_j);
    mK_sj = superJackD(layout,"Main",mK_j);
    
    NPR_sj = NumericTensor<superJackD,2>({7,7}, [&](const int *ij){ return superJackD(layout, "NPR", NPR_j(ij)); } );
    for(int i=0;i<sjack_inputs.size();i++) sjack_inputs[i]->setLayout(layout);
  }
  

};


#endif

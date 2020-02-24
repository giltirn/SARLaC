#ifndef _ANALYZE_KTOPIPI_GPARITY_READ_DATA_H_
#define _ANALYZE_KTOPIPI_GPARITY_READ_DATA_H_

//The data we need to compute the A0 and epsilon', in superJackknife format
struct superMultiData{
  typedef superMultiDistribution<double> superMultiD;
  
  superMultiLayout layout;
  
  superMultiD ainv_sj;
  superMultiD omega_expt_sj;
  superMultiD mod_eps_sj;
  superMultiD ImA2_lat_sj;
  superMultiD ReA2_lat_sj;
  superMultiD ReA2_expt_sj;
  superMultiD ReA0_expt_sj;
  superMultiD delta_2_sj;
  
  NumericTensor<superMultiD,1> M_lat_sj;
  NumericTensor<superMultiD,2> NPR_sj;
  superMultiD Epi_sj;
  superMultiD Epipi_sj;
  superMultiD mpi_sj;
  superMultiD mK_sj;  

  superMultiData(const Args &args){
    //Main ensemble results can be either bootstrap or jackknife
    //Load the matrix elements
    {
      assert(args.fit_results.M_lat.idx.size() == 10); //must be in 10 basis
      DistributionTypeEnum type;
      int depth;
      getTypeInfo(type,depth,args.fit_results.M_lat.file);
      if(depth == 2){
	for(int i=0;i<10;i++) assert(args.fit_results.M_lat.idx[i].size() == 2); //pair of indices for the vectors
	std::vector<std::vector<superMultiD> > fit_params;
	readFromHDF5(fit_params, args.fit_results.M_lat.file, "Main");
	M_lat_sj = NumericTensor<superMultiD,1>({10}, [&](const int *i){
	    int a=args.fit_results.M_lat.idx[*i][0];
	    int b=args.fit_results.M_lat.idx[*i][1];
	    return fit_params[a][b];
	  });
      }else if(depth == 1){
	for(int i=0;i<10;i++) assert(args.fit_results.M_lat.idx[i].size() == 1); //single index
	std::vector<superMultiD> fit_params;
	readFromHDF5(fit_params, args.fit_results.M_lat.file, "Main");
	M_lat_sj = NumericTensor<superMultiD,1>({10}, [&](const int *i){
	    return fit_params[args.fit_results.M_lat.idx[*i][0]];
	  });
      }else{
	error_exit(std::cout << "superMultiData constructor: do not support canonical file formats with depths other than 1 or 2\n");
      }
    }
    std::cout << "Unrenormalized lattice matrix elements:\n";
    for(int q=0;q<10;q++)
      std::cout << "Q" << q+1 << " = " << M_lat_sj(&q) << std::endl;
    
    //Load the inputs needed for the Lellouch Luscher factor
    readFromHDF5(mK_sj,args.fit_results.mK,"mK","Main");
    readFromHDF5(Epi_sj,args.fit_results.Epi,"Epi","Main");
    readFromHDF5(Epipi_sj,args.fit_results.Epipi,"Epipi","Main");
        
    //Load the NPR matrix
    if(args.renormalization.stepscale)
      if(args.renormalization.incG1) 
	NPR_sj = loadNPRstepScaleIncG1(args.renormalization.file, args.renormalization.file_ensB_mu1, args.renormalization.file_ensB_mu2);
      else
	NPR_sj = loadNPRstepScale(args.renormalization.file, args.renormalization.file_ensB_mu1, args.renormalization.file_ensB_mu2);
    else
      NPR_sj = loadNPRmulti(args.renormalization.file); //NPR assumed jackknife
    std::cout << "NPR matrix:\n" << NPR_sj << std::endl;
  
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
    std::vector<superMultiD*> sjack_inputs = {&ainv_sj, &omega_expt_sj, &mod_eps_sj, &ImA2_lat_sj, &ReA2_lat_sj, &ReA2_expt_sj, &ReA0_expt_sj, &delta_2_sj};
    layout = mK_sj.getLayout(); //take a copy
    layout = combine(layout, NPR_sj({0,0}).getLayout());
    for(int i=0;i<sjack_inputs.size();i++) layout = combine(layout, sjack_inputs[i]->getLayout());
    
    std::cout << "Created a superMulti layout with size " << layout.nSamplesTotal() << " comprising " << layout.nEnsembles() << " ensembles:\n";
    for(int i=0;i<layout.nEnsembles();i++) std::cout << '\t' << layout.ensTag(i) << " of size " << layout.nSamplesEns(i) << " and type " << toString(layout.ensType(i)) << std::endl;
  
    //Update the layout for the data on the Main ensemble
    for(int i=0;i<10;i++) M_lat_sj(&i).setLayout(layout);
    Epi_sj.setLayout(layout);
    Epipi_sj.setLayout(layout);
    mK_sj.setLayout(layout);

    //Update the layout for the superjackknife data
    for(int i=0;i<7;i++) for(int j=0;j<7;j++) NPR_sj({i,j}).setLayout(layout);
    for(int i=0;i<sjack_inputs.size();i++) sjack_inputs[i]->setLayout(layout);

    //Compute the pion mass  
    superMultiD ppi2_sj(layout, 0.);
    for(int i=0;i<3;i++)
      if(args.twists[i])
	ppi2_sj = ppi2_sj + superMultiD(layout, pow( args.lattice_dispersion_reln ? sin(M_PI/args.L) : M_PI/args.L, 2 ));  
    std::cout << "Pion p^2 = " << ppi2_sj << std::endl;

    mpi_sj = sqrt( Epi_sj*Epi_sj - ppi2_sj );
    std::cout << "m_pi = " << mpi_sj << std::endl;
  }
  

};


#endif

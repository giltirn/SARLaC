#ifndef _ANALYZE_KTOPIPI_GPARITY_MAIN_H_
#define _ANALYZE_KTOPIPI_GPARITY_MAIN_H_

void checkLatticeWilsonCoefficients(const std::pair<NumericTensor<superMultiDistribution<double>,1>, NumericTensor<superMultiDistribution<double>,1> > &lat_Wilson_coeffs,
				    const NumericTensor<superMultiDistribution<double>,1> &M_lat_sj, //lattice matrix elements
				    const superMultiDistribution<double> &ainv_sj, //inverse lattice spacing
				    const superMultiDistribution<double> &F_sj,  //Lellouch-Luscher factor
				    const Args &args){
  std::cout << "Testing lattice Wilson coefficients are computed correctly" << std::endl;
  typedef superMultiDistribution<double> superMultiD;
  superMultiD coeff = ainv_sj*ainv_sj*ainv_sj*F_sj * args.constants.G_F * args.constants.Vud * args.constants.Vus /sqrt(2.0);
  NumericTensor<superMultiD,1> M_unrenorm_phys_std({10}, [&](const int* c){ return coeff * M_lat_sj(c); });
  NumericTensor<superMultiD,1> M_unrenorm_phys_chiral_sj = convertChiralBasis(M_unrenorm_phys_std);

  superMultiD zero(coeff); zeroit(zero);

  superMultiD reA0(zero), imA0(zero);
  for(int i=0;i<7;i++){
    reA0 = reA0 + lat_Wilson_coeffs.first(&i) *  M_unrenorm_phys_chiral_sj(&i);
    imA0 = imA0 + lat_Wilson_coeffs.second(&i) *  M_unrenorm_phys_chiral_sj(&i);
  }

  std::cout << "Got ReA0 = " << reA0 << "  ImA0 = " << imA0 << std::endl;
}


NumericTensor<superMultiDistribution<double>,1> computePhysicalMSbarMatrixElements(const NumericTensor<superMultiDistribution<double>,1> &M_lat_sj, //lattice matrix elements
										       const superMultiDistribution<double> &ainv_sj, //inverse lattice spacing
										       const superMultiDistribution<double> &F_sj,  //Lellouch-Luscher factor
										       const NumericTensor<superMultiDistribution<double>,2> &NPR_sj, //RI-SMOM renormalization factors
										       const MSbarConvert &MSbar, //MSbar matching coefficients
										       const Args &args){
  typedef superMultiDistribution<double> superMultiD;
  
  //Convert M to physical units and infinite volume
  superMultiD coeff = ainv_sj*ainv_sj*ainv_sj*F_sj;
  NumericTensor<superMultiD,1> M_unrenorm_phys_std({10}, [&](const int* c){ return coeff * M_lat_sj(c); });
  std::cout << "Unrenormalized physical matrix elements (physical units, infinite volume) and standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_unrenorm_phys_std(&q) << " GeV^3\n";
  
  //Convert M to chiral basis
  NumericTensor<superMultiD,1> M_unrenorm_phys_chiral_sj = convertChiralBasis(M_unrenorm_phys_std);

  std::cout << "Unrenormalized physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_unrenorm_phys_chiral_sj(&q)<< " GeV^3\n";

  //Apply NPR
  NumericTensor<superMultiD,1> M_RI_chiral_sj = NPR_sj * M_unrenorm_phys_chiral_sj;

  std::cout << "RI-scheme physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_RI_chiral_sj(&q) << " GeV^3\n";

  writeParamsStandard(M_RI_chiral_sj, "matrix_elems_SMOM.hdf5");

  //Convert to MSbar in 7 basis
  NumericTensor<superMultiD,1> M_MSbar_chiral_sj = MSbar.convert(M_RI_chiral_sj, Chiral);

  std::cout << "MSbar-scheme physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q" << q+1 << " = " << M_MSbar_chiral_sj(&q) << " GeV^3\n";
  
  writeParamsStandard(M_MSbar_chiral_sj, "matrix_elems_MSbar_chiral.hdf5");
  
  //Convert to MSbar
  NumericTensor<superMultiD,1> M_MSbar_std_sj = MSbar.convert(M_RI_chiral_sj, Standard);

  std::cout << "MSbar-scheme physical matrix elements in standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_MSbar_std_sj(&q) << " GeV^3\n";
  
  writeParamsStandard(M_MSbar_std_sj, "matrix_elems_MSbar.hdf5");

  return M_MSbar_std_sj;
}


void computeA0(superMultiDistribution<double> &ReA0_sj,
	       superMultiDistribution<double> &ImA0_sj,
	       const NumericTensor<superMultiDistribution<double>,1> &M_MSbar_std_sj,
	       const WilsonCoeffs &wilson_coeffs,
	       const Args &args,
	       const std::string &file_stub = ""){
  
  typedef superMultiDistribution<double> superMultiD;

  const superMultiLayout &layout = M_MSbar_std_sj({0}).getLayout();
  
  double coeff = args.constants.G_F * args.constants.Vud * args.constants.Vus /sqrt(2.0);

  //Apply Wilson coefficients
  NumericTensor<superMultiD,1> ReA0_cpts_sj({10}), ImA0_cpts_sj({10});
  ReA0_sj = ImA0_sj = superMultiD(layout,[](const int s){return 0.;});
  for(int q=0;q<10;q++){
    std::complex<double> w = wilson_coeffs(q);    
    ReA0_cpts_sj(&q) = coeff * w.real() * M_MSbar_std_sj(&q);
    ImA0_cpts_sj(&q) = coeff * w.imag() * M_MSbar_std_sj(&q);
    ReA0_sj = ReA0_sj + ReA0_cpts_sj(&q);
    ImA0_sj = ImA0_sj + ImA0_cpts_sj(&q);
  }
  std::cout << "Contributions to ReA0 and ImA0:\n";
  for(int q=0;q<10;q++)
    std::cout << q+1 << "   " << ReA0_cpts_sj(&q) << "     " << ImA0_cpts_sj(&q) << std::endl;

  writeParamsStandard(ReA0_cpts_sj, file_stub + "ReA0_cpts.hdf5");
  writeParamsStandard(ImA0_cpts_sj, file_stub + "ImA0_cpts.hdf5");
  
  std::cout << "Total:\n";
  std::cout << ReA0_sj << "     " << ImA0_sj << std::endl;

  writeParamsStandard(ReA0_sj, file_stub+"ReA0.hdf5");
  writeParamsStandard(ImA0_sj, file_stub+"ImA0.hdf5");

#define Q(I) ImA0_cpts_sj({I-1})
  superMultiDistribution<double> EWP_d_QCDP = (Q(7)+Q(8)+Q(9)+Q(10))/(Q(3)+Q(4)+Q(5)+Q(6));
  std::cout << "Ratio of EW penguins to QCD penguin contributions to ImA0 = " << EWP_d_QCDP << std::endl;

  writeParamsStandard(EWP_d_QCDP, file_stub+"ImA0_contrib_EWP_div_QCDP.hdf5");
#undef Q
}

//A0 = \sum_{ij} w_i^{MS} C^{MS<-RI (10x7)}_{ij} M^{RI}_j   =  \sum_j (\sum_i w_i^{MS} C^{MS<-RI (10x7)}_{ij}) M^{RI}_j
//Define W_j = (\sum_i w_i^{MS} C^{MS<-RI (10x7)}_{ij})   as the "RI Wilson coefficients"
NumericTensor<double,1> computeRIwilsonCoefficients(const int reIm,
						    const MSbarConvert &MSbar,
						    const WilsonCoeffs &wilson_coeffs){
  const NumericTensor<double,2> &C_MSbar = MSbar.getConversionMatrix(Standard);
  NumericTensor<double,1> RI_Wilson_coeffs({7}, 0.);

  for(int j=0;j<7;j++)
    for(int i=0;i<10;i++){
      std::complex<double> w = wilson_coeffs(i);   
      RI_Wilson_coeffs(&j) = RI_Wilson_coeffs(&j) + (reIm == 1 ? w.imag() : w.real()) * C_MSbar({i,j});
    }
  return RI_Wilson_coeffs;
}

//A0 = \sum_{ij} w_i^{RI} Z^{RI}_{ij} M^lat_{j}
//Define W_j = \sum_i w_i^{RI} Z^{RI}_{ij} as the 'Lattice Wilson coefficients' (chiral basis) in analog to the above
NumericTensor<superMultiDistribution<double>,1> computeLatticeWilsonCoefficients(const NumericTensor<superMultiDistribution<double>,2> &NPR_sj,
										     const NumericTensor<double,1> &RI_Wilson_coeffs){
  typedef superMultiDistribution<double> superMultiD;
  superMultiD zero(NPR_sj({0,0}).getLayout(),[](const int s){return 0.;});
  NumericTensor<superMultiD,1> lat_Wilson_coeffs({7}, zero);
  for(int j=0;j<7;j++)
    for(int i=0;i<7;i++)
      lat_Wilson_coeffs(&j) = lat_Wilson_coeffs(&j) + RI_Wilson_coeffs(&i) * NPR_sj({i,j});
  return lat_Wilson_coeffs;
}

void computeRIandLatticeWilsonCoefficients(std::pair<NumericTensor<double,1>, NumericTensor<double,1> > &RI_Wilson_coeffs,
					   std::pair<NumericTensor<superMultiDistribution<double>,1>, NumericTensor<superMultiDistribution<double>,1> > &lat_Wilson_coeffs,
					   const WilsonCoeffs &wilson_coeffs,
					   const MSbarConvert &MSbar,
					   const NumericTensor<superMultiDistribution<double>,2> &NPR,
					   const std::string &file_stub){
  typedef superMultiDistribution<double> superMultiD;

  std::cout << "Computing RI and lattice Wilson coefficients\n";
  RI_Wilson_coeffs  = std::pair<NumericTensor<double,1>, NumericTensor<double,1> >( computeRIwilsonCoefficients(0,MSbar,wilson_coeffs), computeRIwilsonCoefficients(1,MSbar,wilson_coeffs) );
  for(int reim=0; reim<2; reim++)
    std::cout << (reim==0 ? "Real" : "Imaginary") << " part RI Wilson coefficients:\n" << (reim == 0 ? RI_Wilson_coeffs.first : RI_Wilson_coeffs.second ) << std::endl;
    
  lat_Wilson_coeffs = std::pair<NumericTensor<superMultiD,1>, NumericTensor<superMultiD,1> >( computeLatticeWilsonCoefficients(NPR, RI_Wilson_coeffs.first), computeLatticeWilsonCoefficients(NPR, RI_Wilson_coeffs.second) );
  for(int reim=0; reim<2; reim++)
    std::cout << (reim==0 ? "Real" : "Imaginary") << " part lattice Wilson coefficients:\n" << (reim == 0 ? lat_Wilson_coeffs.first : lat_Wilson_coeffs.second ) << std::endl;

  {
    std::string filename = file_stub + "real_RI_Wilson_coeffs.dat";
    std::ofstream f(filename.c_str());
    f << RI_Wilson_coeffs.first;
  }
  {
    std::string filename = file_stub + "imag_RI_Wilson_coeffs.dat";
    std::ofstream f(filename.c_str());
    f << RI_Wilson_coeffs.second;
  }
  writeParamsStandard(lat_Wilson_coeffs.first, file_stub+"real_lat_Wilson_coeffs.hdf5");
  writeParamsStandard(lat_Wilson_coeffs.second, file_stub+"imag_lat_Wilson_coeffs.hdf5");
}


NumericTensor<superMultiDistribution<double>,1> computePhysicalMSbarMatrixElementsUsingReA0expt(const NumericTensor<superMultiDistribution<double>,1> &M_lat_sj, //lattice matrix elements
												const superMultiDistribution<double> &ainv_sj, //inverse lattice spacing
												const superMultiDistribution<double> &F_sj,  //Lellouch-Luscher factor
												const NumericTensor<superMultiDistribution<double>,2> &NPR_sj, //RI-SMOM renormalization factors
												const MSbarConvert &MSbar, //MSbar matching coefficients
												const NumericTensor<superMultiDistribution<double>,1> &ReA0_lat_Wilson_coeffs,
												const superMultiDistribution<double> &ReA0_expt_sj,
												const int elim_idx,
												const Args &args,
												const std::string &file_stub){
  typedef superMultiDistribution<double> superMultiD;
  std::cout << "Eliminating Q'" << chiralBasisIdx(elim_idx) << " using experimental ReA0\n";
    
  //Convert M to physical units and infinite volume
  superMultiD coeff = ainv_sj*ainv_sj*ainv_sj*F_sj;
  NumericTensor<superMultiD,1> M_unrenorm_phys_std({10}, [&](const int* c){ return coeff * M_lat_sj(c); });
  std::cout << "Unrenormalized physical matrix elements (physical units, infinite volume) and standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_unrenorm_phys_std(&q) << " GeV^3\n";
  
  //Convert M to chiral basis
  NumericTensor<superMultiD,1> M_unrenorm_phys_chiral_sj = convertChiralBasis(M_unrenorm_phys_std);

  std::cout << "Unrenormalized physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_unrenorm_phys_chiral_sj(&q)<< " GeV^3\n";

  //Replace one of the chiral basis lattice matrix elements with the linear combination of Re(A0) from expt and the other matrix elements

  double prefactor = args.constants.G_F * args.constants.Vud * args.constants.Vus /sqrt(2.0);

  const superMultiLayout &layout = NPR_sj({0,0}).getLayout();
  superMultiD sum_other(layout,[](const int s){return 0.;});

  for(int i=0;i<7;i++)
    if(i!=elim_idx)
      sum_other = sum_other + prefactor * ReA0_lat_Wilson_coeffs(&i) * M_unrenorm_phys_chiral_sj(&i);

  //Perform the elimination
  M_unrenorm_phys_chiral_sj(&elim_idx) = ( ReA0_expt_sj - sum_other ) / ReA0_lat_Wilson_coeffs(&elim_idx) / prefactor;

  std::cout << "Unrenormalized physical Q'" << chiralBasisIdx(elim_idx) << " obtained by include Re(A0) from expt: " << M_unrenorm_phys_chiral_sj(&elim_idx) << std::endl;

  //Apply NPR
  NumericTensor<superMultiD,1> M_RI_chiral_sj = NPR_sj * M_unrenorm_phys_chiral_sj;

  std::cout << "RI-scheme physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_RI_chiral_sj(&q) << " GeV^3\n";

  writeParamsStandard(M_RI_chiral_sj, file_stub+"matrix_elems_SMOM.hdf5");

  //Convert to MSbar in 7 basis
  NumericTensor<superMultiD,1> M_MSbar_chiral_sj = MSbar.convert(M_RI_chiral_sj, Chiral);

  std::cout << "MSbar-scheme physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q" << q+1 << " = " << M_MSbar_chiral_sj(&q) << " GeV^3\n";
  
  writeParamsStandard(M_MSbar_chiral_sj, file_stub+"matrix_elems_MSbar_chiral.hdf5");

  //Convert to MSbar in 10 basis
  NumericTensor<superMultiD,1> M_MSbar_std_sj = MSbar.convert(M_RI_chiral_sj, Standard);

  std::cout << "MSbar-scheme physical matrix elements in standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_MSbar_std_sj(&q) << " GeV^3\n";
  
  writeParamsStandard(M_MSbar_std_sj, file_stub+"matrix_elems_MSbar.hdf5");

  return M_MSbar_std_sj;
}



void computeEpsilon(const superMultiDistribution<double> &ReA0_lat_sj,
		    const superMultiDistribution<double> &ImA0_lat_sj,
		    const superMultiDistribution<double> &ReA2_lat_sj,
		    const superMultiDistribution<double> &ImA2_lat_sj,

		    const superMultiDistribution<double> &delta_0_lat_sj,
		    const superMultiDistribution<double> &delta_2_lat_sj,
		    
		    const superMultiDistribution<double> &ReA0_expt_sj,
		    const superMultiDistribution<double> &ReA2_expt_sj,
		    const superMultiDistribution<double> &omega_expt_sj,
		    const superMultiDistribution<double> &mod_eps_expt_sj,
		    
		    const Args &args,
		    const std::string &file_stub = ""){

		    
  typedef superMultiDistribution<double> superMultiD;
		    
  //Compute lattice omega
  superMultiD omega_lat_sj = ReA2_lat_sj / ReA0_lat_sj;
  std::cout << "w = ReA2/ReA0 lattice = " << omega_lat_sj << " vs expt = " << omega_expt_sj << std::endl;

  const superMultiLayout &layout = ReA0_lat_sj.getLayout();
  
  //Compute the epsilon' phase coefficient
  superMultiD re_phase(layout,
		   [&](const int s){ return cos(delta_2_lat_sj.osample(s) - delta_0_lat_sj.osample(s) + M_PI/2 - args.constants.phi_epsilon); });

  std::cout << "cos(delta2-delta0+pi/2-phi_epsilon) = " << re_phase << std::endl;

  superMultiD coeff_exptRe_sj = re_phase * omega_expt_sj /sqrt(2.0)/ mod_eps_expt_sj;
  superMultiD coeff_latRe_sj = re_phase * omega_lat_sj /sqrt(2.0)/ mod_eps_expt_sj;
  
  //Compute epsilon'/epsilon
  superMultiD ep_div_e_EWP_exptRe = coeff_exptRe_sj * ImA2_lat_sj/ReA2_expt_sj;
  superMultiD ep_div_e_QCDP_exptRe = -coeff_exptRe_sj * ImA0_lat_sj/ReA0_expt_sj;
  superMultiD ep_div_e_exptRe = ep_div_e_EWP_exptRe + ep_div_e_QCDP_exptRe;

  std::cout << "Using ReA0, ReA2, omega from experiment:\n";
  std::cout << "Re(eps'/eps)_EWP = " << ep_div_e_EWP_exptRe << std::endl;
  std::cout << "Re(eps'/eps)_QCDP = " << ep_div_e_QCDP_exptRe << std::endl;
  std::cout << "Re(eps'/eps) = " << ep_div_e_exptRe << "\n\n";

  writeParamsStandard(ep_div_e_EWP_exptRe, file_stub+"ep_div_e_EWP_exptRe.hdf5");
  writeParamsStandard(ep_div_e_QCDP_exptRe, file_stub+"ep_div_e_QCDP_exptRe.hdf5");
  writeParamsStandard(ep_div_e_exptRe, file_stub+"ep_div_e_exptRe.hdf5");
  
  superMultiD ep_div_e_EWP_latRe = coeff_latRe_sj * ImA2_lat_sj/ReA2_lat_sj;
  superMultiD ep_div_e_QCDP_latRe = -coeff_latRe_sj * ImA0_lat_sj/ReA0_lat_sj;
  superMultiD ep_div_e_latRe = ep_div_e_EWP_latRe + ep_div_e_QCDP_latRe;

  std::cout << "Using ReA0, ReA2, omega from lattice:\n";
  std::cout << "Re(eps'/eps)_EWP = " << ep_div_e_EWP_latRe << std::endl;
  std::cout << "Re(eps'/eps)_QCDP = " << ep_div_e_QCDP_latRe << std::endl;
  std::cout << "Re(eps'/eps) = " << ep_div_e_latRe << "\n\n";

  writeParamsStandard(ep_div_e_EWP_latRe, file_stub+"ep_div_e_EWP_latRe.hdf5");
  writeParamsStandard(ep_div_e_QCDP_latRe, file_stub+"ep_div_e_QCDP_latRe.hdf5");
  writeParamsStandard(ep_div_e_latRe, file_stub+"ep_div_e_latRe.hdf5");
}


void computeMatrixElementRelations(const NumericTensor<superMultiDistribution<double>,1> &M, const std::string &file_stub = ""){
#define Q(I) M({I-1})
#undef Q
}


#endif

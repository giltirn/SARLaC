#ifndef _ANALYZE_KTOPIPI_GPARITY_MAIN_H_
#define _ANALYZE_KTOPIPI_GPARITY_MAIN_H_


NumericTensor<superJackknifeDistribution<double>,1> computePhysicalMSbarMatrixElements(const NumericTensor<superJackknifeDistribution<double>,1> &M_lat_sj, //lattice matrix elements
							       const superJackknifeDistribution<double> &ainv_sj, //inverse lattice spacing
							       const superJackknifeDistribution<double> &F_sj,  //Lellouch-Luscher factor
							       const NumericTensor<superJackknifeDistribution<double>,2> &NPR_sj, //RI-SMOM renormalization factors
							       const Args &args){
  typedef superJackknifeDistribution<double> superJackD;
  
  //Compute the MSbar matching coefficients
  MSbarConvert MSbar(args.renormalization.mu,args.renormalization.scheme);
  
  //Convert M to physical units and infinite volume
  superJackD coeff = ainv_sj*ainv_sj*ainv_sj*F_sj * args.constants.G_F * args.constants.Vud * args.constants.Vus /sqrt(2.0);
  NumericTensor<superJackD,1> M_unrenorm_phys_std({10}, [&](const int* c){ return coeff * M_lat_sj(c); });
  std::cout << "Unrenormalized physical matrix elements (physical units, infinite volume) and standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_unrenorm_phys_std(&q) << " GeV^3\n";
  
  //Convert M to chiral basis
  NumericTensor<superJackD,1> M_unrenorm_phys_chiral_sj = convertChiralBasis(M_unrenorm_phys_std);

  std::cout << "Unrenormalized physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_unrenorm_phys_chiral_sj(&q)<< " GeV^3\n";

  //Apply NPR
  NumericTensor<superJackD,1> M_RI_chiral_sj = NPR_sj * M_unrenorm_phys_chiral_sj;

  std::cout << "RI-scheme physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_RI_chiral_sj(&q) << " GeV^3\n";
  
  //Convert to MSbar
  NumericTensor<superJackD,1> M_MSbar_std_sj = MSbar.convert(M_RI_chiral_sj);

  std::cout << "MSbar-scheme physical matrix elements in standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_MSbar_std_sj(&q) << " GeV^3\n";
  
  return M_MSbar_std_sj;
}


void computeA0(superJackknifeDistribution<double> &ReA0_sj,
	       superJackknifeDistribution<double> &ImA0_sj,
	       const NumericTensor<superJackknifeDistribution<double>,1> &M_MSbar_std_sj,
	       const Args &args,
	       const std::string &file_stub = ""){
  
  typedef superJackknifeDistribution<double> superJackD;

  //Read the Wilson coefficients
  WilsonCoeffs wilson_coeffs(args.wilson_coeffs_file);

  const superJackknifeLayout &layout = M_MSbar_std_sj({0}).getLayout();
  
  //Apply Wilson coefficients
  NumericTensor<superJackD,1> ReA0_cpts_sj({10}), ImA0_cpts_sj({10});
  ReA0_sj = ImA0_sj = superJackD(layout,[](const int s){return 0.;});
  for(int q=0;q<10;q++){
    std::complex<double> w = wilson_coeffs(q);    
    ReA0_cpts_sj(&q) = w.real() * M_MSbar_std_sj(&q);
    ImA0_cpts_sj(&q) = w.imag() * M_MSbar_std_sj(&q);
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
}


void computeEpsilon(const superJackknifeDistribution<double> &ReA0_lat_sj,
		    const superJackknifeDistribution<double> &ImA0_lat_sj,
		    const superJackknifeDistribution<double> &ReA2_lat_sj,
		    const superJackknifeDistribution<double> &ImA2_lat_sj,

		    const superJackknifeDistribution<double> &delta_0_lat_sj,
		    const superJackknifeDistribution<double> &delta_2_lat_sj,
		    
		    const superJackknifeDistribution<double> &ReA0_expt_sj,
		    const superJackknifeDistribution<double> &ReA2_expt_sj,
		    const superJackknifeDistribution<double> &omega_expt_sj,
		    const superJackknifeDistribution<double> &mod_eps_expt_sj,
		    
		    const Args &args,
		    const std::string &file_stub = ""){

		    
  typedef superJackknifeDistribution<double> superJackD;
		    
  //Compute lattice omega
  superJackD omega_lat_sj = ReA2_lat_sj / ReA0_lat_sj;
  std::cout << "w = ReA2/ReA0 lattice = " << omega_lat_sj << " vs expt = " << omega_expt_sj << std::endl;

  const superJackknifeLayout &layout = ReA0_lat_sj.getLayout();
  
  //Compute the epsilon' phase coefficient
  superJackD re_phase(layout,
		   [&](const int s){ return cos(delta_2_lat_sj.osample(s) - delta_0_lat_sj.osample(s) + M_PI/2 - args.constants.phi_epsilon); });

  std::cout << "cos(delta2-delta0+pi/2-phi_epsilon) = " << re_phase << std::endl;

  superJackD coeff_exptRe_sj = re_phase * omega_expt_sj /sqrt(2.0)/ mod_eps_expt_sj;
  superJackD coeff_latRe_sj = re_phase * omega_lat_sj /sqrt(2.0)/ mod_eps_expt_sj;
  
  //Compute epsilon'/epsilon
  superJackD ep_div_e_EWP_exptRe = coeff_exptRe_sj * ImA2_lat_sj/ReA2_expt_sj;
  superJackD ep_div_e_QCDP_exptRe = -coeff_exptRe_sj * ImA0_lat_sj/ReA0_expt_sj;
  superJackD ep_div_e_exptRe = ep_div_e_EWP_exptRe + ep_div_e_QCDP_exptRe;

  std::cout << "Using ReA0, ReA2, omega from experiment:\n";
  std::cout << "Re(eps'/eps)_EWP = " << ep_div_e_EWP_exptRe << std::endl;
  std::cout << "Re(eps'/eps)_QCDP = " << ep_div_e_QCDP_exptRe << std::endl;
  std::cout << "Re(eps'/eps) = " << ep_div_e_exptRe << "\n\n";

  writeParamsStandard(ep_div_e_EWP_exptRe, file_stub+"ep_div_e_EWP_exptRe.hdf5");
  writeParamsStandard(ep_div_e_QCDP_exptRe, file_stub+"ep_div_e_QCDP_exptRe.hdf5");
  writeParamsStandard(ep_div_e_exptRe, file_stub+"ep_div_e_exptRe.hdf5");
  
  superJackD ep_div_e_EWP_latRe = coeff_latRe_sj * ImA2_lat_sj/ReA2_lat_sj;
  superJackD ep_div_e_QCDP_latRe = -coeff_latRe_sj * ImA0_lat_sj/ReA0_lat_sj;
  superJackD ep_div_e_latRe = ep_div_e_EWP_latRe + ep_div_e_QCDP_latRe;

  std::cout << "Using ReA0, ReA2, omega from lattice:\n";
  std::cout << "Re(eps'/eps)_EWP = " << ep_div_e_EWP_latRe << std::endl;
  std::cout << "Re(eps'/eps)_QCDP = " << ep_div_e_QCDP_latRe << std::endl;
  std::cout << "Re(eps'/eps) = " << ep_div_e_latRe << "\n\n";

  writeParamsStandard(ep_div_e_EWP_latRe, file_stub+"ep_div_e_EWP_latRe.hdf5");
  writeParamsStandard(ep_div_e_QCDP_latRe, file_stub+"ep_div_e_QCDP_latRe.hdf5");
  writeParamsStandard(ep_div_e_latRe, file_stub+"ep_div_e_latRe.hdf5");
}

#endif

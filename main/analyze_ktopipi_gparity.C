#include<distribution.h>
#include<common_defs.h>
#include<superjackknife.h>
#include<alpha_s.h>
#include<parser.h>

#include<analyze_ktopipi_gparity/params.h>
#include<analyze_ktopipi_gparity/LLfactor.h>
#include<analyze_ktopipi_gparity/renormalization.h>
#include<analyze_ktopipi_gparity/WilsonCoeffs.h>
#include<analyze_ktopipi_gparity/utils.h>

int main(const int argc, const char* argv[]){
  typedef superJackknifeDistribution<double> superJackD;

  if(argc != 2){
    std::cout << "To use ./analyze_ktopipi_gparity <params_file>\n";
    std::cout << "Writing template to template.args\n";
    Args tp;
    std::ofstream of("template.args");
    of << tp;
    of.close();    
    exit(0);
  }
  Args args;
  parse(args, argv[1]);

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
  superJackD ainv_sj; readFromXML(ainv_sj,args.other_inputs.ainv,"a^{-1}");
  superJackD omega_expt_sj; readFromXML(omega_expt_sj,args.other_inputs.omega_expt,"omega_expt");
  superJackD mod_eps_sj; readFromXML(mod_eps_sj,args.other_inputs.mod_eps,"mod_eps");
  superJackD ImA2_lat_sj; readFromXML(ImA2_lat_sj,args.other_inputs.ImA2_lat,"ImA2_lat");
  superJackD ReA2_lat_sj; readFromXML(ReA2_lat_sj,args.other_inputs.ReA2_lat,"ReA2_lat");
  superJackD ReA2_expt_sj; readFromXML(ReA2_expt_sj,args.other_inputs.ReA2_expt,"ReA2_expt");
  superJackD ReA0_expt_sj; readFromXML(ReA0_expt_sj,args.other_inputs.ReA0_expt,"ReA0_expt");
  superJackD delta_2_sj; readFromXML(delta_2_sj,args.other_inputs.delta_2_lat,"delta_2");
  //Create a superjackknife layout for treating all these data consistently
  std::vector<superJackD*> sjack_inputs = {&ainv_sj, &omega_expt_sj, &mod_eps_sj, &ImA2_lat_sj, &ReA2_lat_sj, &ReA2_expt_sj, &ReA0_expt_sj, &delta_2_sj};
  superJackknifeLayout layout;
  {
    layout.addEnsemble("Main", M_lat_j({0}).size());
    layout.addEnsemble("NPR", NPR_j({0,0}).size() );
    for(int i=0;i<sjack_inputs.size();i++) layout = combine(layout, sjack_inputs[i]->getLayout());
  }

  std::cout << "Created a superJackknife layout with size " << layout.nSamplesTotal() << " comprising " << layout.nEnsembles() << " ensembles:\n";
  for(int i=0;i<layout.nEnsembles();i++) std::cout << '\t' << layout.ensTag(i) << " of size " << layout.nSamplesEns(i) << std::endl;
  
  //Upcast everything to a superjackknife
  NumericTensor<superJackD,1> M_lat_sj({10}, [&](const int *i){ return superJackD(layout,"Main",M_lat_j(i)); } );
  superJackD Epi_sj(layout,"Main",Epi_j);
  superJackD Epipi_sj(layout,"Main",Epipi_j);
  superJackD mpi_sj(layout,"Main",mpi_j);
  superJackD mK_sj(layout,"Main",mK_j);
  
  NumericTensor<superJackD,2> NPR_sj({7,7}, [&](const int *ij){ return superJackD(layout, "NPR", NPR_j(ij)); } );
  for(int i=0;i<sjack_inputs.size();i++) sjack_inputs[i]->setLayout(layout);
  
  //Compute the MSbar matching coefficients
  MSbarConvert MSbar(args.renormalization.mu,args.renormalization.scheme);
  
  //Read the Wilson coefficients
  WilsonCoeffs wilson_coeffs(args.wilson_coeffs_file);
  
  //Compute the phase shift and it's derivative wrt q
  LuscherZeta zeta(args.twists[0],args.twists[1],args.twists[2]);

  superJackD p_pipi_sj = sqrt( Epipi_sj*Epipi_sj/4 - mpi_sj*mpi_sj );
  superJackD q_pipi_sj = args.L * p_pipi_sj /( 2 * M_PI );
  std::cout << "p = " << p_pipi_sj << std::endl;
  std::cout << "q = " << q_pipi_sj << std::endl;
  
  superJackD delta_0_sj = getPhaseShift(zeta,q_pipi_sj);
  superJackD delta_0_deg_sj = delta_0_sj/M_PI * 180;

  std::cout << "delta_0 = " << delta_0_sj << " rad = " << delta_0_deg_sj << "deg\n";

  superJackD d_delta_by_dq_schenk_sj = getPhaseShiftDerivSchenk(ainv_sj,Epipi_sj,q_pipi_sj,mpi_sj,args.L);
  superJackD d_delta_by_dq_lin_Epipi_sj = getPhaseShiftDerivLinearEpipi(Epipi_sj,q_pipi_sj,mpi_sj,delta_0_sj,args.L);
  superJackD d_delta_by_dq_lin_qpipi_sj = getPhaseShiftDerivLinearQpipi(q_pipi_sj,delta_0_sj);

  std::cout << "ddelta_0/dq (Schenk) = " << d_delta_by_dq_schenk_sj << std::endl;
  std::cout << "ddelta_0/dq (Lin. Epipi) = " << d_delta_by_dq_lin_Epipi_sj << std::endl;
  std::cout << "ddelta_0/dq (Lin. q) = " << d_delta_by_dq_lin_qpipi_sj << std::endl;

  superJackD d_delta_by_dq_sj;
  switch(args.deriv_source){
  case DerivSchenk:
    d_delta_by_dq_sj = d_delta_by_dq_schenk_sj; break;
  case DerivLinearEpipi:
    d_delta_by_dq_sj = d_delta_by_dq_lin_Epipi_sj; break;
  case DerivLinearQpipi:
    d_delta_by_dq_sj = d_delta_by_dq_lin_qpipi_sj; break;
  default:
    error_exit(std::cout << "Unknown phase shift derivative source " << args.deriv_source << std::endl);
  }
  std::cout << "Using ddelta_0/dq source " << args.deriv_source << ", value " << d_delta_by_dq_sj << std::endl;

  superJackD F_sj = computeLLfactor(Epipi_sj,q_pipi_sj,p_pipi_sj,mK_sj,d_delta_by_dq_sj,zeta,args.L);

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

  //Apply Wilson coefficients
  NumericTensor<superJackD,1> ReA0_cpts_sj({10}), ImA0_cpts_sj({10});
  superJackD ReA0_sj(layout,[](const int s){return 0.;}), ImA0_sj(layout,[](const int s){return 0.;});
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

  writeParamsStandard(ReA0_cpts_sj, "ReA0_cpts.hdf5");
  writeParamsStandard(ImA0_cpts_sj, "ImA0_cpts.hdf5");
  
  std::cout << "Total:\n";
  std::cout << ReA0_sj << "     " << ImA0_sj << std::endl;

  writeParamsStandard(ReA0_sj, "ReA0.hdf5");
  writeParamsStandard(ImA0_sj, "ImA0.hdf5");
  
  //Compute lattice omega
  superJackD omega_lat_sj = ReA2_lat_sj / ReA0_sj;
  std::cout << "w = ReA2/ReA0 lattice = " << omega_lat_sj << " vs expt = " << omega_expt_sj << std::endl;

  //Compute the epsilon' phase coefficient
  superJackD re_phase(layout,
		   [&](const int s){ return cos(delta_2_sj.osample(s) - delta_0_sj.osample(s) + M_PI/2 - args.constants.phi_epsilon); });

  std::cout << "cos(delta2-delta0+pi/2-phi_epsilon) = " << re_phase << std::endl;

  superJackD coeff_exptRe_sj = re_phase * omega_expt_sj /sqrt(2.0)/ mod_eps_sj;
  superJackD coeff_latRe_sj = re_phase * omega_lat_sj /sqrt(2.0)/ mod_eps_sj;
  
  //Compute epsilon'/epsilon
  superJackD ep_div_e_EWP_exptRe = coeff_exptRe_sj * ImA2_lat_sj/ReA2_expt_sj;
  superJackD ep_div_e_QCDP_exptRe = -coeff_exptRe_sj * ImA0_sj/ReA0_expt_sj;
  superJackD ep_div_e_exptRe = ep_div_e_EWP_exptRe + ep_div_e_QCDP_exptRe;

  std::cout << "Using ReA0, ReA2, omega from experiment:\n";
  std::cout << "Re(eps'/eps)_EWP = " << ep_div_e_EWP_exptRe << std::endl;
  std::cout << "Re(eps'/eps)_QCDP = " << ep_div_e_QCDP_exptRe << std::endl;
  std::cout << "Re(eps'/eps = " << ep_div_e_exptRe << "\n\n";

  writeParamsStandard(ep_div_e_EWP_exptRe, "ep_div_e_EWP_exptRe.hdf5");
  writeParamsStandard(ep_div_e_QCDP_exptRe, "ep_div_e_QCDP_exptRe.hdf5");
  writeParamsStandard(ep_div_e_exptRe, "ep_div_e_exptRe.hdf5");
  
  superJackD ep_div_e_EWP_latRe = coeff_latRe_sj * ImA2_lat_sj/ReA2_lat_sj;
  superJackD ep_div_e_QCDP_latRe = -coeff_latRe_sj * ImA0_sj/ReA0_sj;
  superJackD ep_div_e_latRe = ep_div_e_EWP_latRe + ep_div_e_QCDP_latRe;

  std::cout << "Using ReA0, ReA2, omega from lattice:\n";
  std::cout << "Re(eps'/eps)_EWP = " << ep_div_e_EWP_latRe << std::endl;
  std::cout << "Re(eps'/eps)_QCDP = " << ep_div_e_QCDP_latRe << std::endl;
  std::cout << "Re(eps'/eps = " << ep_div_e_latRe << "\n\n";

  writeParamsStandard(ep_div_e_EWP_latRe, "ep_div_e_EWP_latRe.hdf5");
  writeParamsStandard(ep_div_e_QCDP_latRe, "ep_div_e_QCDP_latRe.hdf5");
  writeParamsStandard(ep_div_e_latRe, "ep_div_e_latRe.hdf5");
  
  std::cout << "Done\n";
  return 0;
}

#ifndef LL_FACTOR_H
#define LL_FACTOR_H

#include<physics.h>


superMultiDistribution<double> getPhaseShift(const LuscherZeta &lz, const superMultiDistribution<double> &q_pipi){
  std::cout << "Computing phase shift from lattice values\n";
  superMultiDistribution<double> phi(q_pipi.getLayout(), [&](const int i){ return lz.calcPhi(q_pipi.osample(i)); } );

  std::cout << "phi(q) = " << phi << std::endl;
  
  superMultiDistribution<double> delta(q_pipi.getLayout(), 
				       [&](const int i){
					 double d = -phi.osample(i);
					 if(d<0.) d += ceil(-d/M_PI)*M_PI;
					 if(d>M_PI) d -= floor(d/M_PI)*M_PI;
					 return d;
				       });
    
  std::cout << "delta = " << delta << std::endl;
  return delta;
}

superMultiDistribution<double> getPhaseShiftDerivSchenk(const superMultiDistribution<double> &ainv,
							    const superMultiDistribution<double> &Epipi,
							    const superMultiDistribution<double> &q_pipi,
							    const superMultiDistribution<double> &mpi,
							    const int L){					     
  //Pheno curve is in terms of Mandelstam variable s = E_pipi^2,
  //q = pL/2pi     p = ( Epipi^2 /4 - mpi^2 )^{1/2}
  //p = (s/4 - mpi^2)^{1/2}
  //2pi/L q = (s/4 - mpi^2)^{1/2}
  //2pi/L dq/ds = (s/4 - mpi^2)^{-1/2} * 1/2*1/4 = (8*2*pi/L q)^{-1} 
  // dq/ds = (8* [2*pi/L]^2 * q)^{-1} 

  //or
  //2pi/L q = ( Epipi^2 /4 - mpi^2 )^{1/2}
  //2pi/L dq/dEpipi = ( Epipi^2 /4 - mpi^2 )^{-1/2} * 1/2*2*Epipi/4
  //          =  1/(2pi/L q) * Epipi/4
  //dq/dEpipi = Epipi/(4 q) * 1/(2pi/L)^2

  // ds/dq = ds/dEpipi * dEpipi/dq = 2*Epipi / (dq/dEpipi) =   2*Epipi / [ Epipi/(4 q) * 1/(2pi/L)^2 ] = 8 q (2pi/L)^2 

  //For the pheno curve, everything needs to be in units of *MeV*

  std::cout << "Computing phase shift derivative using Schenk formulae\n";
  return superMultiDistribution<double>(ainv.getLayout(),
					    [&](const int b){
					      double ainv_b = ainv.osample(b)* 1000; //in MeV
					      double Epipi_b = Epipi.osample(b) * ainv_b;
					      double qpipi_b = q_pipi.osample(b); //q is dimensionless
					      double mpi_b = mpi.osample(b) * ainv_b;
					      double Lphys = double(L)/ainv_b; // L_latt = L_phys/a, L_phys = a* L_latt = L_latt/a^{-1}
					      double s = Epipi_b*Epipi_b;
					      //double dq_by_dEpipi = Epipi_b /(4*qpipi_b) * pow(2*M_PI/Lphys,-2);
					      //double ds_by_dq = 2*Epipi_b / dq_by_dEpipi;
					      double ds_by_dq = 8 * pow(2*M_PI/Lphys,2) * qpipi_b;
					      double ddelta_by_ds = PhenoCurveSchenk::compute_deriv(s,0,mpi_b,'B',1e-02); //radians/(MeV^2)
					      double ddelta_by_dq = ddelta_by_ds * ds_by_dq;
					      if(b==-1){
						std::cout << "a^{-1} = " << ainv_b << " MeV,  Epipi = " << Epipi_b << " MeV, q =" << qpipi_b << " (dimensionless), mpi = " << mpi_b << " MeV, s = " << s
							  << " MeV,  L = " << Lphys << " MeV\n";
						std::cout << "ds/dq = " << ds_by_dq << " MeV^2\n";
						std::cout << "Computed d(delta)/ds = " << ddelta_by_ds << " MeV^{-2}\n";
						std::cout << "Computed d(delta)/dq = " << ddelta_by_ds * ds_by_dq << std::endl;
					      }
					      return ddelta_by_dq;});
}


superMultiDistribution<double> getPhaseShiftDerivColangelo(const superMultiDistribution<double> &ainv,
							       const superMultiDistribution<double> &Epipi,
							       const superMultiDistribution<double> &q_pipi,
							       const superMultiDistribution<double> &mpi,
							       const int L){					     
  std::cout << "Computing phase shift derivative using Colangelo formulae\n";
  return superMultiDistribution<double>(ainv.getLayout(),
					    [&](const int b){
					      double ainv_b = ainv.osample(b)* 1000; //in MeV
					      double Epipi_b = Epipi.osample(b) * ainv_b;
					      double qpipi_b = q_pipi.osample(b); //q is dimensionless
					      double mpi_b = mpi.osample(b) * ainv_b;
					      double Lphys = double(L)/ainv_b; // L_latt = L_phys/a, L_phys = a* L_latt = L_latt/a^{-1}
					      double s = Epipi_b*Epipi_b;
					      double ds_by_dq = 8 * pow(2*M_PI/Lphys,2) * qpipi_b;
					      double ddelta_by_ds = PhenoCurveColangelo::compute_deriv(s,0,mpi_b,1e-02); //radians/(MeV^2)
					      double ddelta_by_dq = ddelta_by_ds * ds_by_dq;
					      if(b==-1){
						std::cout << "a^{-1} = " << ainv_b << " MeV,  Epipi = " << Epipi_b << " MeV, q =" << qpipi_b << " (dimensionless), mpi = " << mpi_b << " MeV, s = " << s
							  << " MeV,  L = " << Lphys << " MeV\n";
						std::cout << "ds/dq = " << ds_by_dq << " MeV^2\n";
						std::cout << "Computed d(delta)/ds = " << ddelta_by_ds << " MeV^{-2}\n";
						std::cout << "Computed d(delta)/dq = " << ddelta_by_ds * ds_by_dq << std::endl;
					      }
					      return ddelta_by_dq;});
}

//While Colangelo's formula is given in terms of m_pi it is not clear whether it is correct for describing non-physical pion masses
//Edit 1/30/20: It is *NOT* correct to use Colangelo's formular for non-physical pion masses. Use this one!
superMultiDistribution<double> getPhaseShiftDerivColangeloPhysMpi(const superMultiDistribution<double> &ainv,
							       const superMultiDistribution<double> &Epipi,
							       const superMultiDistribution<double> &q_pipi,							       
							       const int L){
  std::cout << "Computing phase shift derivative using Colangelo formulae with physical pion mass as input\n";
  return superMultiDistribution<double>(ainv.getLayout(),
					    [&](const int b){
					      double ainv_b = ainv.osample(b)* 1000; //in MeV
					      double Epipi_b = Epipi.osample(b) * ainv_b;
					      double qpipi_b = q_pipi.osample(b); //q is dimensionless
					      double mpi_b = 135; //MeV!
					      double Lphys = double(L)/ainv_b; // L_latt = L_phys/a, L_phys = a* L_latt = L_latt/a^{-1}
					      double s = Epipi_b*Epipi_b;
					      double ds_by_dq = 8 * pow(2*M_PI/Lphys,2) * qpipi_b;
					      double ddelta_by_ds = PhenoCurveColangelo::compute_deriv(s,0,mpi_b,1e-02); //radians/(MeV^2)
					      double ddelta_by_dq = ddelta_by_ds * ds_by_dq;
					      if(b==-1){
						std::cout << "a^{-1} = " << ainv_b << " MeV,  Epipi = " << Epipi_b << " MeV, q =" << qpipi_b << " (dimensionless), mpi = " << mpi_b << " MeV, s = " << s
							  << " MeV,  L = " << Lphys << " MeV\n";
						std::cout << "ds/dq = " << ds_by_dq << " MeV^2\n";
						std::cout << "Computed d(delta)/ds = " << ddelta_by_ds << " MeV^{-2}\n";
						std::cout << "Computed d(delta)/dq = " << ddelta_by_ds * ds_by_dq << std::endl;
					      }
					      return ddelta_by_dq;});
}




superMultiDistribution<double> getPhaseShiftDerivLinearEpipi(const superMultiDistribution<double> &Epipi,
								 const superMultiDistribution<double> &q_pipi,
								 const superMultiDistribution<double> &mpi,
								 const superMultiDistribution<double> &delta,
								 const int L){
  //We assume a linear function  f( Epipi ) = C*( Epipi - 2*mpi )
  std::cout << "Computing phase shift derivative using linear derivative in Epipi\n";
      
  superMultiDistribution<double> ddelta_by_Epipi = delta / (Epipi - 2.0*mpi);
  //2pi/L q = ( Epipi^2 /4 - mpi^2 )^{1/2}
  //2pi/L dq/dEpipi = ( Epipi^2 /4 - mpi^2 )^{-1/2} * 1/2*2*Epipi/4
  //          =  1/(2pi/L q) * Epipi/4
  //dq/dEpipi = Epipi/(4 q) * 1/(2pi/L)^2
  double Ld(L);
  superMultiDistribution<double> dq_by_dEpipi = Epipi/( 4.*q_pipi ) / (4*M_PI*M_PI/Ld/Ld); 
  return ddelta_by_Epipi / dq_by_dEpipi;
}

superMultiDistribution<double> getPhaseShiftDerivLinearQpipi(const superMultiDistribution<double> &q_pipi,
								 const superMultiDistribution<double> &delta){								 
  std::cout << "Computing phase shift derivative using linear derivative in q\n";
  return delta/q_pipi;
}

superMultiDistribution<double> computedPhiDerivativeQ(const superMultiDistribution<double> &q_pipi,
							  const LuscherZeta &zeta){
  std::cout << "Computing Lellouch-Luscher factor Phi derivative term\n";
  superMultiDistribution<double> dphi_by_dq(q_pipi.getLayout(),
						[&](const int b){ return zeta.calcPhiDeriv(b==-1 ? q_pipi.best() : q_pipi.sample(b)); });

  std::cout << "dphi(q)/dq = " << dphi_by_dq << std::endl;
  return dphi_by_dq;
}

superMultiDistribution<double> computeLLfactor(const superMultiDistribution<double> &Epipi,						   
						   const superMultiDistribution<double> &q_pipi,
						   const superMultiDistribution<double> &p_pipi,
						   const superMultiDistribution<double> &mK,
						   const superMultiDistribution<double> &ddelta_by_dq,
						   const superMultiDistribution<double> &dphi_by_dq,
						   const int L){
  std::cout << "Computing Lellouch-Luscher factor\n";

  superMultiDistribution<double> ddelta_by_dp = (L/2./M_PI) * ddelta_by_dq;   //d/dp = d/d[2pi/L q] = L/2pi d/dq

  std::cout << "d(delta)/dq = " << ddelta_by_dq << std::endl;
  std::cout << "d(delta)/dp = " << ddelta_by_dp << std::endl;

  //Compute F using Eq.23 of arXiv:1106.2714
  superMultiDistribution<double> F2 = 4*M_PI * Epipi*Epipi * mK / p_pipi/p_pipi/p_pipi * (  p_pipi * ddelta_by_dp + q_pipi * dphi_by_dq );

  superMultiDistribution<double> F = sqrt(F2);
  std::cout << "Obtained LL factor F = " << F << std::endl;
  
  return F;
}

void computePhaseShiftAndLLfactor(superMultiDistribution<double> &delta_0, //I=0 phase shift
				  superMultiDistribution<double> &F, //Lellouch-Luscher factor
				  const superMultiDistribution<double> &Epipi, const superMultiDistribution<double> &mpi,
				  const superMultiDistribution<double> &mK, const superMultiDistribution<double> &ainv, const Args &args, const std::string &file_stub = ""){
  typedef superMultiDistribution<double> superMultiD; 
  
  //Compute the phase shift and it's derivative wrt q
  LuscherZeta zeta({args.twists[0],args.twists[1],args.twists[2]},  {0.,0.,0.});

  superMultiD p_pipi = sqrt( Epipi*Epipi/4 - mpi*mpi );
  superMultiD q_pipi = args.L * p_pipi /( 2 * M_PI );
  std::cout << "p = " << p_pipi << std::endl;
  std::cout << "q = " << q_pipi << std::endl;

  superMultiD dphi_by_dq = computedPhiDerivativeQ(q_pipi, zeta);
  
  delta_0 = getPhaseShift(zeta,q_pipi);
  superMultiD delta_0_deg = delta_0/M_PI * 180;

  std::cout << "delta_0 = " << delta_0 << " rad = " << delta_0_deg << "deg\n";

  writeParamsStandard(delta_0_deg, file_stub+"delta0_deg.hdf5");

  superMultiD d_delta_by_dq[5] = { getPhaseShiftDerivSchenk(ainv,Epipi,q_pipi,mpi,args.L),
				   getPhaseShiftDerivColangelo(ainv,Epipi,q_pipi,mpi,args.L),
				   getPhaseShiftDerivLinearEpipi(Epipi,q_pipi,mpi,delta_0,args.L),
				   getPhaseShiftDerivLinearQpipi(q_pipi,delta_0),
				   getPhaseShiftDerivColangeloPhysMpi(ainv,Epipi,q_pipi,args.L) };
  
  PhaseShiftDerivativeSource d_delta_by_dq_src[5] = { PhaseShiftDerivativeSource::DerivSchenk, PhaseShiftDerivativeSource::DerivColangelo, 
						      PhaseShiftDerivativeSource::DerivLinearEpipi, PhaseShiftDerivativeSource::DerivLinearQpipi,
						      PhaseShiftDerivativeSource::DerivColangeloPhysMpi };

  std::string d_delta_by_dq_descr[5] = { "Schenk", "Colangelo", "Lin. Epipi", "Lin. q", "Colangelo phys. mpi." };
  
  int i_use = -1;
  for(int i=0;i<5;i++){
    superMultiD F_i = computeLLfactor(Epipi,q_pipi,p_pipi,mK,d_delta_by_dq[i],dphi_by_dq,args.L);
    std::cout << "ddelta_0/dq (" << d_delta_by_dq_descr[i] << ") = " << d_delta_by_dq[i] << std::endl;
    std::cout << "F (" << d_delta_by_dq_descr[i] << ") = " << F_i << std::endl;

    if(d_delta_by_dq_src[i] == args.deriv_source){
      F = F_i; i_use = i;
    }
  }
  if(i_use == -1) error_exit(std::cout << "computePhaseShiftAndLLfactor: Unknown phase shift derivative source\n");

  std::cout << "Using F=" << F << " from source " << d_delta_by_dq_src[i_use] << " for final result\n";
}



#endif

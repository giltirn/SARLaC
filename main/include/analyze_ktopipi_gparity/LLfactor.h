#ifndef LL_FACTOR_H
#define LL_FACTOR_H

#include<luscherzeta.h>

superJackknifeDistribution<double> getPhaseShift(const LuscherZeta &lz, const superJackknifeDistribution<double> &q_pipi){
  std::cout << "Computing phase shift from lattice values\n";
  superJackknifeDistribution<double> phi(q_pipi.getLayout(), [&](const int i){ return lz.calcPhi(i==-1 ? q_pipi.best() : q_pipi.sample(i)); } );

  std::cout << "phi(q) = " << phi << std::endl;
  
  superJackknifeDistribution<double> delta(q_pipi.getLayout(), [&](const int i){
      double d = -phi.osample(i);
      if(d<0.) d += ceil(-d/M_PI)*M_PI;
      if(d>M_PI) d -= floor(d/M_PI)*M_PI;
      return d;
    });
    
  std::cout << "delta = " << delta << std::endl;
  return delta;
}

superJackknifeDistribution<double> getPhaseShiftDerivSchenk(const superJackknifeDistribution<double> &ainv,
							    const superJackknifeDistribution<double> &Epipi,
							    const superJackknifeDistribution<double> &q_pipi,
							    const superJackknifeDistribution<double> &mpi,
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
  return superJackknifeDistribution<double>(ainv.getLayout(),
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
					      double ddelta_by_ds = PhenoCurve::compute_deriv(s,0,mpi_b,'B',1e-02); //radians/(MeV^2)
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


superJackknifeDistribution<double> getPhaseShiftDerivLinearEpipi(const superJackknifeDistribution<double> &Epipi,
								 const superJackknifeDistribution<double> &q_pipi,
								 const superJackknifeDistribution<double> &mpi,
								 const superJackknifeDistribution<double> &delta,
								 const int L){
  //We assume a linear function  f( Epipi ) = C*( Epipi - 2*mpi )
  std::cout << "Computing phase shift derivative using linear derivative in Epipi\n";
      
  superJackknifeDistribution<double> ddelta_by_Epipi = delta / (Epipi - 2.0*mpi);
  //2pi/L q = ( Epipi^2 /4 - mpi^2 )^{1/2}
  //2pi/L dq/dEpipi = ( Epipi^2 /4 - mpi^2 )^{-1/2} * 1/2*2*Epipi/4
  //          =  1/(2pi/L q) * Epipi/4
  //dq/dEpipi = Epipi/(4 q) * 1/(2pi/L)^2
  double Ld(L);
  superJackknifeDistribution<double> dq_by_dEpipi = Epipi/( 4.*q_pipi ) / (4*M_PI*M_PI/Ld/Ld); 
  return ddelta_by_Epipi / dq_by_dEpipi;
}

superJackknifeDistribution<double> getPhaseShiftDerivLinearQpipi(const superJackknifeDistribution<double> &q_pipi,
								 const superJackknifeDistribution<double> &delta){								 
  std::cout << "Computing phase shift derivative using linear derivative in q\n";
  return delta/q_pipi;
}

superJackknifeDistribution<double> computeLLfactor(const superJackknifeDistribution<double> &Epipi,						   
						   const superJackknifeDistribution<double> &q_pipi,
						   const superJackknifeDistribution<double> &p_pipi,
						   const superJackknifeDistribution<double> &mK,
						   const superJackknifeDistribution<double> &ddelta_by_dq,
						   const LuscherZeta &zeta,
						   const int L){
  std::cout << "Computing Lellouch-Luscher factor\n";
  superJackknifeDistribution<double> dphi_by_dq(q_pipi.getLayout(),
						[&](const int b){ return zeta.calcPhiDeriv(b==-1 ? q_pipi.best() : q_pipi.sample(b)); });

  std::cout << "dphi(q)/dq = " << dphi_by_dq << std::endl;
  superJackknifeDistribution<double> ddelta_by_dp = (L/2./M_PI) * ddelta_by_dq;   //d/dp = d/d[2pi/L q] = L/2pi d/dq

  std::cout << "d(delta)/dq = " << ddelta_by_dq << std::endl;
  std::cout << "d(delta)/dp = " << ddelta_by_dp << std::endl;

  //Compute F using Eq.23 of arXiv:1106.2714
  superJackknifeDistribution<double> F2 = 4*M_PI * Epipi*Epipi * mK / p_pipi/p_pipi/p_pipi * (  p_pipi * ddelta_by_dp + q_pipi * dphi_by_dq );

  superJackknifeDistribution<double> F = sqrt(F2);
  std::cout << "Obtained LL factor F = " << F << std::endl;
  
  return F;
}



#endif

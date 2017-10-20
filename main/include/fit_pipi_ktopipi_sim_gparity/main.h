#ifndef FIT_PIPI_KTOPIPI_SIM_GPARITY_MAIN_H_
#define FIT_PIPI_KTOPIPI_SIM_GPARITY_MAIN_H_

typedef correlationFunction<amplitudeDataCoord, jackknifeDistributionD> jackAmplitudeCorrelationFunction;
typedef correlationFunction<amplitudeDataCoord, doubleJackknifeDistributionD> doubleJackAmplitudeCorrelationFunction;
typedef correlationFunction<amplitudeDataCoordSim, jackknifeDistributionD> jackAmplitudeSimCorrelationFunction;
typedef correlationFunction<amplitudeDataCoordSim, doubleJackknifeDistributionD> doubleJackAmplitudeSimCorrelationFunction;


void convertChiralBasisAndInsertKtoPipi(jackAmplitudeSimCorrelationFunction &data_combined_j,
					doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,				 
					const std::vector<jackAmplitudeCorrelationFunction> &A0_all_j,
					const std::vector<doubleJackAmplitudeCorrelationFunction> &A0_all_dj){
  const int ndata_q = A0_all_j[0].size();
  for(int i=1;i<10;i++)
    if(A0_all_j[i].size() != ndata_q)
      error_exit(std::cout << "convertChiralBasisAndInsert expected all matrix elements to have the same number of data points\n");

  if(ndata_q == 0) return;
  
  const int nsample = A0_all_j[0].value(0).size();
  
  //Convert the data to the chiral basis and gather into single containers    
  //Convert Q123 -> Q'123
  static const double Q123rot[3][3] = {  { 3    ,  2,    -1     },
					 { 2./5 , -2./5,  1./5  },
					 {-3./5,   3./5,  1./5  } };
  
  jackknifeDistributionD zero_j(nsample,0.);
  doubleJackknifeDistributionD zero_dj(nsample,0.);
      
  for(int d=0;d<ndata_q;d++){
    for(int i=1;i<10;i++) assert(A0_all_j[i].coord(d) == A0_all_j[0].coord(d)); //should all be at the same coordinate, just different q
      
    std::vector<jackknifeDistributionD> Qprime_j(7, zero_j);
    std::vector<doubleJackknifeDistributionD> Qprime_dj(7, zero_dj);
#define Q(i,j) Q123rot[i-1][j-1]
      
#define MOj(i) Qprime_j[i-1]
#define MIj(i) A0_all_j[i-1].value(d)
#define MOdj(i) Qprime_dj[i-1]
#define MIdj(i) A0_all_dj[i-1].value(d)      
      
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
      amplitudeDataCoordSim cc(A0_all_j[0].coord(d), q);	
      data_combined_j.push_back(jackAmplitudeSimCorrelationFunction::ElementType(cc,Qprime_j[q]));
      data_combined_dj.push_back(doubleJackAmplitudeSimCorrelationFunction::ElementType(cc,Qprime_dj[q]));
    }
  }
}


void insertPipi(jackAmplitudeSimCorrelationFunction &data_combined_j,
		doubleJackAmplitudeSimCorrelationFunction &data_combined_dj,
		const doubleJackCorrelationFunction &pipi_data_dj){
  const int dummy = 0;
  const int idx = -1;
  for(int i=0;i<pipi_data_dj.size();i++){
    amplitudeDataCoordSim c(pipi_data_dj.coord(i),dummy,idx);
    data_combined_j.push_back(jackAmplitudeSimCorrelationFunction::ElementType(c,pipi_data_dj.value(i).toJackknife()));
    data_combined_dj.push_back(doubleJackAmplitudeSimCorrelationFunction::ElementType(c,pipi_data_dj.value(i)));
  }
}

void vectorizeAndConvert10basis(std::vector<jackknifeDistributionD> &into,		    
				const jackknifeDistribution<TwoPointThreePointSimFitParams> &params){
  const int nsample = params.size();
  
  into.resize(15,jackknifeDistributionD(nsample,0.));
#pragma omp parallel for
  for(int s=0;s<nsample;s++){
    for(int i=0;i<5;i++)
      into[i].sample(s) = params.sample(s)(i);
  
    //Convert Q'123 -> Q123
    static const double Q123invrot[3][3] = {  {1./5,   1,   0},
					      {1./5,   0,   1},
					      {  0 ,   3,   2} };
#define MO(i) into[i+5-1].sample(s)
#define MI(i) params.sample(s)(i+5-1)
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
}

#endif

#ifndef _FIT_KTOPIPI_KTOSIGMA_GPARITY_PROJFIT_H
#define _FIT_KTOPIPI_KTOSIGMA_GPARITY_PROJFIT_H

#include<config.h>
#include<utils/macros.h>

CPSFIT_START_NAMESPACE

void projectionFit(const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktopipi_A0_all_j,
		   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &ktopipi_A0_all_dj,
		   const std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > &ktosigma_A0_all_j,
		   const std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > &ktosigma_A0_all_dj,
		   const jackknifeDistributionD &mK, const jackknifeDistributionD &cK, 
		   const jackknifeDistributionD &E0, const jackknifeDistributionD &E1,
		   const NumericSquareMatrix<jackknifeDistributionD> &coeffs, //row = (0=pipi, 1=sigma)  col = (0=gnd state, 1=exc state)
		   const int Lt, const int tmin_k_op, const int tmin_op_snk, bool correlated){
  int nsample = mK.size();
  /*
    Define linear combination a O_0 + b O_1   s.t.   
    a <O_0|0> + b<O_1|0> = 1
    a <O_0|1> + b<O_1|1> = 0

    Thus
    a = -b<O_1|1>/<O_0|1> 
    
    -b<O_1|1><O_0|0>/<O_0|1> + b<O_1|0> = 1

  */
  jackknifeDistributionD one(nsample,1.);
  jackknifeDistributionD b = one/( coeffs(1,0) - coeffs(1,1)*coeffs(0,0)/coeffs(0,1) );
  jackknifeDistributionD a = -b * coeffs(1,1)/coeffs(0,1);
    
  std::cout << "Combination coefficients " << a << " " << b << std::endl;

  //Compute linear combination
  std::vector<correlationFunction<amplitudeDataCoord, jackknifeDistributionD> > A0_all_j(10);
  std::vector<correlationFunction<amplitudeDataCoord, doubleJackknifeA0StorageType> > A0_all_dj(10);  

  for(int q=0;q<10;q++){
    assert( ktosigma_A0_all_j[q].size() == ktopipi_A0_all_j[q].size() );

    int N = ktopipi_A0_all_j[q].size();

    A0_all_j[q].resize(N);
    A0_all_dj[q].resize(N);

    for(int i=0; i< N; i++){
      assert( ktosigma_A0_all_j[q].coord(i) == ktopipi_A0_all_j[q].coord(i) );
      A0_all_j[q].coord(i) = A0_all_dj[q].coord(i) = ktosigma_A0_all_j[q].coord(i);
      A0_all_j[q].value(i) = a * ktopipi_A0_all_j[q].value(i) + b * ktosigma_A0_all_j[q].value(i);

      //Coeffs are single-jackknife. We don't want to do a triple jackknife to get double-jackknife coefficients so perhaps the best we can do is use the same
      //coefficient for all samples of a given outer index

      A0_all_dj[q].value(i).resize(nsample);
      for(int s=0;s<nsample;s++){
	A0_all_dj[q].value(i).sample(s) = a.sample(s) * ktopipi_A0_all_dj[q].value(i).sample(s) + b.sample(s) * ktosigma_A0_all_dj[q].value(i).sample(s);
      }
    }
  }


  KtoPiPiFitFunc fitfunc = KtoPiPiFitFunc::FitSeparate;
  FitKtoPiPiFreezeFixed freeze;
  freeze.addFreeze(0, cK);   
  freeze.addFreeze(1, mK);
  freeze.addFreeze(2, one); //unit-normalized the combined pipi/sigma operator
  freeze.addFreeze(3, E0);


  FitKtoPiPiOptions opt;
  opt.import_freeze_data = true;
  opt.freeze_data_mgr = &freeze;

  fitAndPlot(A0_all_j,A0_all_dj, Lt, tmin_k_op, tmin_op_snk, fitfunc, correlated, opt);
}

CPSFIT_END_NAMESPACE

#endif

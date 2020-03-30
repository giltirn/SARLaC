#ifndef ANALYZE_KTOPIPI_CONVERT_CHIRAL_H_
#define ANALYZE_KTOPIPI_CONVERT_CHIRAL_H_

//Covert index 0...6  to chiral basis index  1,2,3,5,6,7,8
inline int chiralBasisIdx(const int q){ return q<=2 ? q+1 : q+2; }

//Convert from Q -> Q'' basis (following Qi's convention, pg 66 of his thesis), or the inverse if 'reverse' is true
template<typename DistributionType>
NumericTensor<DistributionType,1> convertChiralBasis(const NumericTensor<DistributionType,1> &Min, const bool reverse = false){
  //Convert Q123 -> Q'123
  static const double Q123rot[3][3] = {  { 3    ,  2,    -1     },
					 { 2./5 , -2./5,  1./5  },
					 {-3./5,   3./5,  1./5  } };
  //Convert Q'123 -> Q123
  static const double Q123invrot[3][3] = {  {1./5,   1,   0},
					    {1./5,   0,   1},
					    {  0 ,   3,   2} };

#define MO(i) Mout({i-1})
#define MI(i) Min({i-1})
#define Q(i,j) Q123rot[i-1][j-1]
#define Qinv(i,j) Q123invrot[i-1][j-1]
    
  if(!reverse){
    //We throw away Q4, Q9 and Q10 but in practise they might differ from their linear combinations due to the random numbers used
    //Can we find some way to include an average to take advantage of these data?
    NumericTensor<DistributionType,1> Mout({7});
    for(int i=1;i<=3;i++)
      MO(i) = Q(i,1)*MI(1) + Q(i,2)*MI(2) + Q(i,3)*MI(3);
    for(int i=4;i<=7;i++)
      MO(i) = MI(i+1); //5->4  6->5 etc
    return Mout;
  }else{
    NumericTensor<DistributionType,1> Mout({10});
      
    for(int i=1;i<=3;i++)
      MO(i) = Qinv(i,1)*MI(1) + Qinv(i,2)*MI(2) + Qinv(i,3)*MI(3);

    MO(4) = MO(2) + MO(3) - MO(1); //Q4 = Q2 + Q3 - Q1    [Lehner, Sturm, arXiv:1104.4948 eq 9]
    for(int i=5;i<=8;i++) MO(i) = MI(i-1); //4->5 5->6 etc
	
    MO(9) = 3./2*MO(1)  -1./2*MO(3); //Q9 = 3/2 Q1 - 1/2 Q3
    MO(10) = 1./2*MO(1)  -1./2*MO(3) + MO(2); //Q10 = 1/2(Q1 - Q3) + Q2
    return Mout;
  }
#undef MO
#undef MI
#undef Q
#undef Qinv
}


//Convert to chiral basis by choosing optimal combination of data that minimizes
//q^2 = ( T_aj Q'_j - Q_a ) cov^-1_ab ( T_bj Q'_j - Q_b )
//where cov is the 10x10 covariance matrix


struct ConvChiralFitFunc{
  typedef NumericTensor<double,2> MatrixType;

  typedef NumericVector<double> Params; //7 coeffs
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef int GeneralizedCoordinate; //index in 10 basis

  MatrixType T7to10;

  ConvChiralFitFunc(){
    //T matrix maps between chiral basis and conventional basis
    const static double tin[] =
      {  1./5, 1, 0, 0, 0, 0, 0,
	 1./5, 0, 1, 0, 0, 0, 0,
	 0   , 3, 2, 0, 0, 0, 0,
	 0   , 2, 3, 0, 0, 0, 0,
	 0   , 0, 0, 1, 0, 0, 0,
	 0   , 0, 0, 0, 1, 0, 0,
	 0   , 0, 0, 0, 0, 1, 0,
	 0   , 0, 0, 0, 0, 0, 1,
	 3./10,  0, -1, 0, 0, 0, 0,
	 3./10, -1,  0, 0, 0, 0, 0 };
    T7to10 = MatrixType({10,7}, [](const int *c){ return tin[c[1]+7*c[0]]; });
  }


  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate i, const ParameterType &p, const Boost &b) const{
    T out(0.);
    for(int j=0;j<7;j++)
      out = out + T7to10({i, j})*b(p(j), j);
    return out;
  }


  inline ValueType value(const GeneralizedCoordinate x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<7;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  int Nparams() const{ return 7; }
};



template<typename DistributionType>
NumericTensor<DistributionType,1> convertChiralFit(const NumericTensor<DistributionType,1> &Min){
  assert(Min.size(0) == 10);

  auto binit = Min({0}).getInitializer();

  NumericVector<double> guess(7, 1.0);
  ConvChiralFitFunc fitfunc;
  genericFitFuncWrapper<ConvChiralFitFunc> fwrap(fitfunc, guess);

  simpleFitWrapper<DistributionType> fitter(fwrap, MinimizerType::MarquardtLevenberg);
  
  DistributionType chisq(binit,0.), chisq_per_dof(binit,0.);
  int dof;
  typename DistributionType::template rebase< NumericVector<double> > params(binit, guess);

  correlationFunction<int, DistributionType> Q10(10);
  for(int i=0;i<10;i++){
    Q10.value(i) = Min(&i);      
    Q10.coord(i) = i;
  }

  NumericSquareMatrix<DistributionType> cov(10, [&](const int i, const int j){ return DistributionType(binit, DistributionType::covariance(Min(&i), Min(&j)) ); });
  std::cout << "Covariance matrix of 10 operator matrix elements:" << std::endl;
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++)
      std::cout << cov(i,j).best() << " ";
    std::cout << std::endl;
  }


  fitter.importCovarianceMatrix(cov);
  
  std::cout << "Correlation matrix of 10 operator matrix elements:" << std::endl;
  for(int i=0;i<10;i++){
    for(int j=0;j<10;j++)
      std::cout << fitter.getCorrelationMatrix()(i,j).best() << " ";
    std::cout << std::endl;
  }

  fitter.fit(params, chisq, chisq_per_dof, dof, Q10);
  
  int qmap[7] = {1,2,3,5,6,7,8};
  
  std::vector<DistributionType> Q7_fit(7);
  for(int i=0;i<7;i++)
    Q7_fit[i] = distributionStructPeek(params, i);
  
  //Test we reproduce the 10-basis results
  std::vector<DistributionType> Q10_repro(10, DistributionType(binit,0.));
  for(int i=0;i<10;i++)
    for(int j=0;j<7;j++) 
      Q10_repro[i] = Q10_repro[i] + fitfunc.T7to10({i,j}) * Q7_fit[j];
  
  writeParamsStandard(Q10_repro, "matrix_elems_unrenorm_std_repro_from_7basis.hdf5");

  std::cout << "Q in 10 basis reproduced from 7-basis result" << std::endl;
  for(int i=0;i<10;i++){
    DistributionType diff = Q10_repro[i] - Q10.value(i);
    std::cout << i+1 << " " << Q10_repro[i] << " diff " << diff << std::endl;
  }

  std::cout << "Chisq: " << chisq << std::endl;
  std::cout << "Dof: " << dof << std::endl;
  std::cout << "Chisq/dof: " << chisq_per_dof << std::endl;

  std::cout << "Results:" << std::endl;
  for(int i=0;i<7;i++) std::cout << qmap[i] << " " << Q7_fit[i] << std::endl;

  NumericTensor<DistributionType,1> out({7});
  for(int i=0;i<7;i++) out(&i) = std::move( Q7_fit[i] );
  return out;
}



NumericTensor<superMultiDistribution<double>,1> computePhysicalUnrenormalizedMatrixElementsChiral(const NumericTensor<superMultiDistribution<double>,1> &M_lat_sj, //lattice matrix elements
												  const superMultiDistribution<double> &ainv_sj, //inverse lattice spacing
												  const superMultiDistribution<double> &F_sj,  //Lellouch-Luscher factor
												  const Args &args){
  assert(M_lat_sj.size(0) == 10);
  typedef superMultiDistribution<double> superMultiD;
  
  //Convert M to physical units and infinite volume
  superMultiD coeff = ainv_sj*ainv_sj*ainv_sj*F_sj;
  NumericTensor<superMultiD,1> M_unrenorm_phys_std({10}, [&](const int* c){ return coeff * M_lat_sj(c); });
  std::cout << "Unrenormalized physical matrix elements (physical units, infinite volume) and standard basis:\n";
  for(int q=0;q<10;q++)
    std::cout << "Q" << q+1 << " = " << M_unrenorm_phys_std(&q) << " GeV^3\n";
  
  writeParamsStandard(M_unrenorm_phys_std, "matrix_elems_unrenorm_std.hdf5");

  //Convert M to chiral basis
  NumericTensor<superMultiD,1> M_unrenorm_phys_chiral_sj;
  if(args.chiral_conv_method == ChiralConvertMethod::Original){
    M_unrenorm_phys_chiral_sj = convertChiralBasis(M_unrenorm_phys_std);
  }else{
    M_unrenorm_phys_chiral_sj = convertChiralFit(M_unrenorm_phys_std);
  }
  
  std::cout << "Unrenormalized physical matrix elements in chiral basis:\n";
  for(int q=0;q<7;q++)
    std::cout << "Q'" << chiralBasisIdx(q) << " = " << M_unrenorm_phys_chiral_sj(&q)<< " GeV^3\n";
  
  writeParamsStandard(M_unrenorm_phys_chiral_sj, "matrix_elems_unrenorm_chiral.hdf5");
  
  return M_unrenorm_phys_chiral_sj;
}


#endif

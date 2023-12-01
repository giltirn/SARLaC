#ifndef _ANALYZE_KTOPIPI_RENORMALIZATION_H_
#define _ANALYZE_KTOPIPI_RENORMALIZATION_H_

NumericTensor<jackknifeDistributionD,2> loadNPR(const std::string &file, const int size = 7){
  XMLreader rd(file);
  UKvalenceDistributionContainer<jackknifeDistributionD> con;
  read(rd, con, "data_in_file");
  
  NumericTensor<jackknifeDistributionD,2> npr({size,size});
  assert(con.Nentries == size*size);
  
  for(int ij=0;ij<size*size;ij++){
    int j = ij %size;
    int i = ij /size;
    npr({i,j}) = std::move(con.list[ij]);
  }
  return npr;
}

NumericTensor<superMultiDistribution<double>,2> loadNPRmulti(const std::string &file, const std::string &ens_tag = "NPR"){
  NumericTensor<jackknifeDistributionD,2> mat = loadNPR(file);
  superMultiLayout* layout = new superMultiLayout;
  layout->addEnsemble(ens_tag, MultiType::Jackknife, mat({0,0}).size());
  return NumericTensor<superMultiDistribution<double>,2>({7,7}, [&](const int* c){ return superMultiDistribution<double>(*layout, 0, mat(c)); });
}
 
//Compute step-scaling matrix between mu1 and mu2 on ensB and apply to matrix at mu1 on ensA under superjackknife
//If lambda_out != NULL  the step-scaling matrix will be copied to that location
NumericTensor<superMultiDistribution<double>,2> loadNPRstepScale(const std::string &file_mu1_ensA,
								 const std::string &file_mu1_ensB,
								 const std::string &file_mu2_ensB,
								 NumericTensor<superMultiDistribution<double>,2> *lambda_out = NULL,
								 const std::string &ens_tag_stub = "NPR"){
  std::cout << "Computing NPR via step-scaling" << std::endl;
  
  NumericTensor<jackknifeDistributionD,2> mat_mu1_ensA = loadNPR(file_mu1_ensA);
  NumericTensor<jackknifeDistributionD,2> mat_mu1_ensB = loadNPR(file_mu1_ensB);
  NumericTensor<jackknifeDistributionD,2> mat_mu2_ensB = loadNPR(file_mu2_ensB);

  std::cout << "NPR matrix at mu1 on base ensemble:\n" << mat_mu1_ensA << std::endl;
  std::cout << "NPR matrix at mu2 on fine ensemble:\n" << mat_mu2_ensB << std::endl;
  std::cout << "NPR matrix at mu1 on fine ensemble:\n" << mat_mu1_ensB << std::endl;
  
  NumericTensor<jackknifeDistributionD,2> mat_mu1_ensB_inv(mat_mu1_ensB);
  svd_inverse(mat_mu1_ensB_inv, mat_mu1_ensB);
  
  std::cout << "Inverse of NPR matrix at mu1 on fine ensemble:\n" << mat_mu1_ensB_inv << std::endl;

  NumericTensor<jackknifeDistributionD,2> Lambda_mu2_mu1_ensB = mat_mu2_ensB * mat_mu1_ensB_inv;
  
  std::cout << "Non-perturbative running between mu1 and mu2:\n" << Lambda_mu2_mu1_ensB << std::endl;

  superMultiLayout* layout = new superMultiLayout;
  layout->addEnsemble(ens_tag_stub + "_ensA", MultiType::Jackknife, mat_mu1_ensA({0,0}).size());
  layout->addEnsemble(ens_tag_stub + "_ensB", MultiType::Jackknife, mat_mu1_ensB({0,0}).size());

  NumericTensor<superMultiDistribution<double>,2> mat_mu1_ensA_sj({7,7}, [&](const int* c){ return superMultiDistribution<double>(*layout, 0, mat_mu1_ensA(c)); });
  NumericTensor<superMultiDistribution<double>,2> Lambda_mu2_mu1_ensB_sj({7,7}, [&](const int* c){ return superMultiDistribution<double>(*layout, 1, Lambda_mu2_mu1_ensB(c)); });

  if(lambda_out != NULL) *lambda_out = Lambda_mu2_mu1_ensB_sj;

  return NumericTensor<superMultiDistribution<double>,2>( Lambda_mu2_mu1_ensB_sj * mat_mu1_ensA_sj );
}



//Compute step-scaling matrix between mu1 and mu2 on ensB and apply to matrix at mu1 on ensA under superjackknife
//This version includes the G1 operator such that the internal matrices are 8x8 and a conversion is required
//If lambda_out != NULL  the step-scaling matrix will be copied to that location
NumericTensor<superMultiDistribution<double>,2> loadNPRstepScaleIncG1(const std::string &file_mu1_ensA,
								 const std::string &file_mu1_ensB,
								 const std::string &file_mu2_ensB,
								 NumericTensor<superMultiDistribution<double>,2> *lambda_out = NULL,
								 const std::string &ens_tag_stub = "NPR"){
  std::cout << "Computing NPR via step-scaling" << std::endl;
  
  NumericTensor<jackknifeDistributionD,2> mat_mu1_ensA = loadNPR(file_mu1_ensA,8);
  NumericTensor<jackknifeDistributionD,2> mat_mu1_ensB = loadNPR(file_mu1_ensB,8);
  NumericTensor<jackknifeDistributionD,2> mat_mu2_ensB = loadNPR(file_mu2_ensB,8);

  std::cout << "8x8 NPR matrix at mu1 on base ensemble:\n" << mat_mu1_ensA << std::endl;
  std::cout << "8x8 NPR matrix at mu2 on fine ensemble:\n" << mat_mu2_ensB << std::endl;
  std::cout << "8x8 NPR matrix at mu1 on fine ensemble:\n" << mat_mu1_ensB << std::endl;
  
  NumericTensor<jackknifeDistributionD,2> mat_mu1_ensB_inv(mat_mu1_ensB);
  svd_inverse(mat_mu1_ensB_inv, mat_mu1_ensB);
  
  std::cout << "Inverse of NPR matrix at mu1 on fine ensemble:\n" << mat_mu1_ensB_inv << std::endl;

  NumericTensor<jackknifeDistributionD,2> Lambda_mu2_mu1_ensB = mat_mu2_ensB * mat_mu1_ensB_inv;
  
  std::cout << "Non-perturbative running between mu1 and mu2:\n" << Lambda_mu2_mu1_ensB << std::endl;

  superMultiLayout* layout = new superMultiLayout;
  layout->addEnsemble(ens_tag_stub + "_ensA", MultiType::Jackknife, mat_mu1_ensA({0,0}).size());
  layout->addEnsemble(ens_tag_stub + "_ensB", MultiType::Jackknife, mat_mu1_ensB({0,0}).size());

  NumericTensor<superMultiDistribution<double>,2> mat_mu1_ensA_sj({8,8}, [&](const int* c){ return superMultiDistribution<double>(*layout, 0, mat_mu1_ensA(c)); });
  NumericTensor<superMultiDistribution<double>,2> Lambda_mu2_mu1_ensB_sj({8,8}, [&](const int* c){ return superMultiDistribution<double>(*layout, 1, Lambda_mu2_mu1_ensB(c)); });

  if(lambda_out != NULL) *lambda_out = Lambda_mu2_mu1_ensB_sj;

  NumericTensor<superMultiDistribution<double>,2> mat_mu2_stepscale_8x8( Lambda_mu2_mu1_ensB_sj * mat_mu1_ensA_sj );
  
  std::cout << "8x8 NPR matrix including G1:\n" << mat_mu2_stepscale_8x8 << std::endl;

  //To convert to 7x7 use Greg's thesis Eq. 7.22 and 7.24
  //c_i are the first 7 elements of the last column
  //d_i are the first 7 elements of the last row
  //Z_G1 is the (8,8) element  (indexing from 1!)
  const static std::vector<double> s = {0, 1, 7./3, -1./3, 1, 0, 0};
  
  std::vector<superMultiDistribution<double> > k(7, superMultiDistribution<double>(*layout));

  for(int j=0;j<7;j++){
    superMultiDistribution<double> siZij(*layout, 0.);
    for(int i=0;i<7;i++)
      siZij = siZij + s[i]*mat_mu2_stepscale_8x8({i,j});

    superMultiDistribution<double> skck(*layout, 0.);
    for(int k=0;k<7;k++){
      const auto & ck = mat_mu2_stepscale_8x8({k,7});
      skck = skck + s[k]*ck;
    }
    const auto & dj = mat_mu2_stepscale_8x8({7,j});
    const auto & ZG1 = mat_mu2_stepscale_8x8({7,7});

    k[j] = ( siZij - dj )/( ZG1 - skck );
    std::cout << "Computed k_" << j << " = " << k[j] << std::endl;
  }

  NumericTensor<superMultiDistribution<double>,2> mat_mu2_stepscale_7x7({7,7});
  for(int i=0;i<7;i++){
    const auto & ci = mat_mu2_stepscale_8x8({i,7});
    for(int j=0;j<7;j++)
      mat_mu2_stepscale_7x7({i,j}) = mat_mu2_stepscale_8x8({i,j}) + ci*k[j];
  }
  std::cout << "7x7 NPR matrix including G1:\n" << mat_mu2_stepscale_7x7 << std::endl;

  return mat_mu2_stepscale_7x7;
}



 


enum OperatorBasis { Standard, Chiral };

//Compute the Delta S=1 RI->MSbar conversion factors using Lehner, Sturm, Phys.Rev. D84 (2011) 014001 
class MSbarConvert{
  typedef NumericTensor<double,2> MatrixType;
  
  MatrixType conv_std; //10x7 matrix
  MatrixType conv_chiral; //7x7 matrix

  MatrixType dr_norm_gamma_gamma; //7x7
  MatrixType dr_norm_qslash_qslash; //7x7
  MatrixType dr_norm_qslash_gamma; //7x7
  MatrixType dr_norm_gamma_qslash; //7x7

  MatrixType T; //10x7
  MatrixType unit_chiral; //7x7
 
public:

  //O(alpha_s)
  const MatrixType &getdrNorm(const RIscheme scheme) const{ 
    switch(scheme){
    case RIscheme::GammaGamma:
      return dr_norm_gamma_gamma;
    case RIscheme::QslashQslash:
      return dr_norm_qslash_qslash;
    case RIscheme::GammaQslash:
      return dr_norm_gamma_qslash;
    case RIscheme::QslashGamma:
      return dr_norm_qslash_gamma;
    default:
      assert(0);
    }
  }

  const MatrixType &getT() const{ return T; }


  //The matrix Delta T_1^MSbar defined in Eq 65. Requires precomputation of alpha_s   O(alpha_s)
  static inline MatrixType computedT(const double alpha_s_over_4pi){
    const int Nc=3;
    MatrixType dT({10,7}, 0.);
    dT({3,1}) = alpha_s_over_4pi * (3./Nc - 2);
    dT({3,2}) = alpha_s_over_4pi * (2./Nc - 3);
    dT({3,3}) = alpha_s_over_4pi * 1./Nc;
    dT({3,4}) = -alpha_s_over_4pi;
    return dT;
  }

  MSbarConvert(): dr_norm_gamma_gamma({7,7}), dr_norm_qslash_qslash({7,7}), dr_norm_gamma_qslash({7,7}), dr_norm_qslash_gamma({7,7}){
    //Input matrices for compute. These are mu-independent and so can be stored
    unit_chiral = MatrixType({7,7}, [](const int *c){ return c[0]==c[1] ? 1. : 0.; });
  
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
    T = MatrixType({10,7}, [](const int *c){ return tin[c[1]+7*c[0]]; });
    
    //Map from chiral basis indices 1,2,3,5,6,7,8  
    //to standard indices           0,1,2,3,4,5,6
    static const int chiral_basis_map[9] = {-1, 0, 1, 2, -1, 3, 4, 5, 6 }; //God I hate these conventions

    //\Delta_r*/(\alpha_s/(4\pi)) in GammaGamma scheme with Nc=3 and \xi=0 (Landau gauge) from table IX
    {
#define r(i,j) dr_norm_gamma_gamma({chiral_basis_map[i],chiral_basis_map[j]})
      r(1,1) = 0.21184;
      r(2,2) = -0.10592;
      r(2,3) = 0.31777;
      r(3,2) = 0.09554;
      r(3,3) = -0.62444;
      r(3,5) = 0.07407;
      r(3,6) = -0.22222;
      r(5,5) = 0.04319;
      r(5,6) = -0.12957;
      r(6,2) = -1.66667;
      r(6,3) = -3.88889;
      r(6,5) = -1.05815;
      r(6,6) = 2.82894;
      r(7,7) = 0.04319;
      r(7,8) = -0.12957;
      r(8,7) = -1.61371;
      r(8,8) = 4.49561;
#undef r
    }
  
    //\Delta_r*/(\alpha_s/(4\pi)) in QslashQslash scheme with Nc=3 and \xi=0 (Landau gauge) from table IX
    {
#define r(i,j) dr_norm_qslash_qslash({chiral_basis_map[i],chiral_basis_map[j]})
      r(1,1) = -0.45482;
      r(2,2) = 2.62741;
      r(2,3) = 4.91777;
      r(3,2) = -4.50446;
      r(3,3) = -8.69111;
      r(3,5) = 0.07407;
      r(3,6) = -0.22222;
      r(5,5) = 0.04319;
      r(5,6) = -0.12957;
      r(6,2) = -1.66667;
      r(6,3) = -3.88889;
      r(6,5) = -0.05815;
      r(6,6) = -0.17106;
      r(7,7) = 0.04319;
      r(7,8) = -0.12957;
      r(8,7) = -0.61371;
      r(8,8) = 1.49561;
#undef r
    }

    { //Table X
#define r(i,j) dr_norm_qslash_gamma({chiral_basis_map[i],chiral_basis_map[j]})
      r(1,1) = 2.21184;
      r(2,2) = 5.29408;
      r(2,3) = 4.91777;
      r(3,2) = -4.50446;
      r(3,3) = -6.02444;
      r(3,5) = 0.07407;
      r(3,6) = -0.22222;
      r(5,5) = 2.70986;
      r(5,6) = -0.12957;
      r(6,2) = -1.66667;
      r(6,3) = -3.88889;
      r(6,5) = -0.05815;
      r(6,6) = 2.49561;
      r(7,7) = 2.70986;
      r(7,8) = -0.12957;
      r(8,7) = -0.61371;
      r(8,8) = 4.16227;
#undef r
    }

    { //Table VIII
#define r(i,j) dr_norm_gamma_qslash({chiral_basis_map[i],chiral_basis_map[j]})
      r(1,1) = -2.45482;
      r(2,2) = -2.77259;
      r(2,3) = 0.31777;
      r(3,2) = 0.09554;
      r(3,3) = -3.29111;
      r(3,5) = 0.07407;
      r(3,6) = -0.22222;
      r(5,5) = -2.62348;
      r(5,6) = -0.12957;
      r(6,2) = -1.66667;
      r(6,3) = -3.88889;
      r(6,5) = -1.05815;
      r(6,6) = 0.16227;
      r(7,7) = -2.62348;
      r(7,8) = -0.12957;
      r(8,7) = -1.61371;
      r(8,8) = 1.82894;
#undef r
    }
  }
  MSbarConvert(const double mu_GeV, const RIscheme scheme, const PerturbativeVariables &pv): MSbarConvert(){ compute(mu_GeV,scheme,pv); }
  MSbarConvert(const std::string &file, const OperatorBasis basis): MSbarConvert(){ read(file, basis); }

  //Compute the 10x7 matrix that converts from the chiral to 10-basis in MSbar (Nc=3)
  MatrixType computeMSbar7to10basisConv(const double alpha_s_over_4pi) const{
    return MatrixType(T + computedT(alpha_s_over_4pi));
  }


  //Read the MSbar conversion matrix
  void read(const std::string &file, const OperatorBasis basis){
    std::ifstream ff(file.c_str());
    if (!ff.good()) error_exit(std::cout << "MSbarConvert::read failed to open file " << file << std::endl);

    int bsz = basis == Standard ? 10 : 7;
    MatrixType &into = basis == Standard ? conv_std : conv_chiral;

    int c[2];
    for(c[0]=0;c[0]<bsz;c[0]++)
      for(c[1]=0;c[1]<7;c[1]++){
	ff >> into(c);
      }
    assert(!ff.bad());    
    ff.close();

    std::cout << "Read MSbar conversion matrix:\n" << into << std::endl;
  }


  //Compute the MSbar conversion matrix
  void compute(const double alpha_s, const RIscheme scheme, const bool vrb = true){
    double alpha_s_over_4pi = alpha_s / 4. / M_PI;
    const int Nc = 3;
      
    //NOTE: For MSbar you must also include \Delta T in the conversion (here for basis I, the usual 10-basis)
    MatrixType dT = computedT(alpha_s_over_4pi);

    if(vrb) std::cout << "Computing RI->MSbar matrix with RI scheme " << scheme << " and \\alpha_s(mu) = " << alpha_s << '\n';

    const MatrixType &dr_norm = scheme == RIscheme::GammaGamma ? dr_norm_gamma_gamma : dr_norm_qslash_qslash;
    MatrixType R = unit_chiral + dr_norm * alpha_s_over_4pi;
    conv_chiral = R;
    conv_std = (T + dT)*R;

    if(vrb){
      std::cout << "Using \\Delta r/(as/4pi) = \n" << dr_norm << std::endl;
      std::cout << "\n\nR^RI->MSbar(7x7) =\n" << conv_chiral << std::endl;
      std::cout << "\n\nR^RI->MSbar(10x7) =\n" << conv_std << std::endl;
    }
  }

  void compute(const double mu_GeV, const RIscheme scheme, const PerturbativeVariables &pv, const bool vrb = true){
    if(vrb) std::cout << "Computing 2-loop alpha_s at mu=" << mu_GeV << " GeV in MSbar scheme\n";
    double alpha_s = ComputeAlphaS::alpha_s(mu_GeV, pv, 3);
    if(vrb) printf("alpha_s at mu=%f: %.8f\n",mu_GeV,alpha_s);
    compute(alpha_s, scheme, vrb);    
  }

  //Convert the matrix elements in the RI scheme and chiral basis to MSbar matrix elements in the 10-basis or 7-basis
  NumericTensor<superMultiDistribution<double>,1> convert(const NumericTensor<superMultiDistribution<double>,1> &Min, const OperatorBasis basis) const{
    superMultiDistribution<double> zro(Min({0}).getLayout());
    zro.zero();

    int bsz = basis == Standard ? 10 : 7;
    const MatrixType &cmat = basis == Standard ? conv_std : conv_chiral;

    NumericTensor<superMultiDistribution<double>,1> Mout({bsz},zro);
    int c[2];
    for(c[0]=0;c[0]<bsz;c[0]++)
      for(c[1]=0;c[1]<7;c[1]++)
	Mout(&c[0]) = Mout(&c[0]) + cmat(c) * Min(&c[1]);
    return Mout;
  }

  const MatrixType & getConversionMatrix(const OperatorBasis basis) const{ return basis == Standard ? conv_std : conv_chiral; }
};


#endif

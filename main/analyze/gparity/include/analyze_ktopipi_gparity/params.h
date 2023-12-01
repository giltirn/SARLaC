#ifndef _ANALYZE_KTOPIPI_PARAMS_H_
#define _ANALYZE_KTOPIPI_PARAMS_H_

#include<parser.h>

GENERATE_ENUM_AND_PARSER(PhaseShiftDerivativeSource, (DerivSchenk)(DerivLinearEpipi)(DerivLinearQpipi)(DerivColangelo)(DerivColangeloPhysMpi)(DerivColangeloPhysMpiCharged) );
GENERATE_ENUM_AND_PARSER(RIscheme, (QslashQslash)(GammaGamma)(QslashGamma)(GammaQslash) );
GENERATE_ENUM_AND_PARSER(ChiralConvertMethod, (Original)(Fit) );

#define FILE_IDX_PAIR_MEMBERS (std::string, file)(int, idx)(std::string, operation)
struct FileIdxPair{
  GENERATE_MEMBERS(FILE_IDX_PAIR_MEMBERS);
  FileIdxPair(): operation(""){}
  FileIdxPair(const std::string &f, const int i, const std::string &operation = ""): file(f),idx(i), operation(operation){}
};
GENERATE_PARSER(FileIdxPair, FILE_IDX_PAIR_MEMBERS);

#define M_FILE_INFO_MEMBERS (std::string, file)(std::vector<std::vector<int> >, idx)
struct MfileInfo{
  GENERATE_MEMBERS(M_FILE_INFO_MEMBERS);
  MfileInfo(){}
  MfileInfo(const std::string &f, const std::vector<std::vector<int> > &i): file(f),idx(i){}
  MfileInfo(const std::string &f): file(f), idx(10, std::vector<int>(2)){
    for(int i=0;i<10;i++){ //default layout for separate fits to the matrix elements (stored in vector<vector<dist> > format)
      idx[i][0] = i;
      idx[i][1] = 4;
    }
  }
    
};
GENERATE_PARSER(MfileInfo, M_FILE_INFO_MEMBERS);


#define FIT_RESULTS_MEMBERS (MfileInfo, M_lat)(FileIdxPair, mK)(FileIdxPair, Epi)(FileIdxPair, Epipi)
struct FitResults{ //these are all expected to be hdf5 files currently
  GENERATE_MEMBERS(FIT_RESULTS_MEMBERS);

  FitResults(): M_lat("32c_216cfgs_results/ktopipi.hdf5"), mK("32c_216cfgs_results/mk.hdf5", 1), Epi("32c_216cfgs_results/Epi.hdf5", 1), Epipi("32c_216cfgs_results/pipi.hdf5", 1){}
};
GENERATE_PARSER(FitResults, FIT_RESULTS_MEMBERS);

#define OTHER_INPUTS_MEMBERS (FileIdxPair, ainv)(FileIdxPair, omega_expt)(FileIdxPair, mod_eps)(FileIdxPair, ReA0_expt)(FileIdxPair, ReA2_expt)(FileIdxPair, ReA2_lat)(FileIdxPair, ImA2_lat)(FileIdxPair, delta_2_lat)
struct OtherInputs{ //these are all expected to be xml-format superjackknife files currently
  GENERATE_MEMBERS(OTHER_INPUTS_MEMBERS);
  OtherInputs(): ainv("32c_216cfgs_results/ainv_gaussian.bootxml",0),
		 omega_expt("32c_216cfgs_results/fit_inputs/omega.bootxml",0),
		 mod_eps("32c_216cfgs_results/fit_inputs/mod_eps.bootxml",0),
		 ReA0_expt("32c_216cfgs_results/fit_inputs/reA0_expt.bootxml",0),
		 ReA2_expt("32c_216cfgs_results/fit_inputs/reA2_expt.bootxml",0),
		 ReA2_lat("32c_216cfgs_results/fit_inputs/reA2_lat_gaussian.bootxml",0), ImA2_lat("32c_216cfgs_results/fit_inputs/imA2_lat.bootxml",0),
		 delta_2_lat("32c_216cfgs_results/fit_inputs/delta2_gaussian.bootxml",0){}
};
GENERATE_PARSER(OtherInputs,OTHER_INPUTS_MEMBERS);

//If stepscale = true,  mu must be set to the *high scale*, mu2. "File" should be the matrix on the original ensemble (ensA) at mu1,  and filenames for the matrix
//at mu1 and mu2 should be provided for the data on the finer ensemble ensB
//Option incG1 is valid only for step-scaling with the G1 operator included; here the internal matrices are 8x8 and a conversion to the 7x7 matrix is required
#define RENORMALIZATION_MEMBERS (RIscheme, scheme)(double, mu)(std::string, file)(bool, stepscale)(bool, incG1)(std::string, file_ensB_mu1)(std::string, file_ensB_mu2)
struct Renormalization{
  GENERATE_MEMBERS(RENORMALIZATION_MEMBERS);
Renormalization(): scheme(RIscheme::QslashQslash), mu(1.531), file("32c_216cfgs_results/qslash_1.53GeV_noG1.xml"),stepscale(false),incG1(false),file_ensB_mu1(""),file_ensB_mu2(""){}
};
GENERATE_PARSER(Renormalization, RENORMALIZATION_MEMBERS);

#define CONSTANT_MEMBERS (double, G_F)(double, Vud)(double, Vus)(double, phi_epsilon)(std::complex<double>, tau)
struct Constants{ //Currently treated as constants but some could be upgraded to jackknifes to add in errors
  GENERATE_MEMBERS(CONSTANT_MEMBERS);
Constants(): G_F(1.1663787e-05), Vud(0.97425), Vus(0.2253), phi_epsilon(0.75956729046793223188), tau(0.001543, -0.000635){}
};
GENERATE_PARSER(Constants, CONSTANT_MEMBERS);



#define ARGS_MEMBERS \
  (FitResults, fit_results)			\
  (OtherInputs, other_inputs)			\
  (Renormalization, renormalization)		\
  (Constants, constants)			\
  (PerturbativeInputs, perturbative_inputs)	\
  (std::string, wilson_coeffs_file)		\
  (std::vector<int>, twists)			\
  (int, L)					\
  (PhaseShiftDerivativeSource, deriv_source)	\
  (bool, lattice_dispersion_reln)		\
  (ChiralConvertMethod, chiral_conv_method)

struct Args{
  GENERATE_MEMBERS(ARGS_MEMBERS);
  Args(): wilson_coeffs_file("32c_216cfgs_results/WilsonCoeffs.dat"), //cf WilsonCoeffs.h for format
    twists({1,1,1}),
    L(32),
    deriv_source(PhaseShiftDerivativeSource::DerivLinearQpipi),
    lattice_dispersion_reln(false), chiral_conv_method(ChiralConvertMethod::Original){}
};
GENERATE_PARSER(Args, ARGS_MEMBERS);


#endif

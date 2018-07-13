#ifndef _PIPI_SIGMA_SIM_FIT_FITFUNC_H_
#define _PIPI_SIGMA_SIM_FIT_FITFUNC_H_

GENERATE_ENUM_AND_PARSER(SimFitType, (PiPiToSigma)(Sigma2pt)(PiPi2pt) );

#define SIM_FIT_COORD_ARGS \
  (SimFitType, type)	   \
  (double, t)

struct SimFitCoord{
  GENERATE_MEMBERS(SIM_FIT_COORD_ARGS);  

  SimFitCoord() = default;
  SimFitCoord(const SimFitType type, const double t): type(type), t(t){}
};

GENERATE_PARSER(SimFitCoord, SIM_FIT_COORD_ARGS);

typedef correlationFunction<SimFitCoord, doubleJackknifeDistributionD> simFitCorrFuncDJ;
typedef correlationFunction<SimFitCoord, jackknifeDistributionD> simFitCorrFuncJ;

#define SIM_FIT_PARAMS_MEMBERS (cpipi0)(cpipi1)(csigma0)(csigma1)(E0)(E1)(Cpipipipi)(Cpipisigma)(Csigmasigma)
DEF_ENUMERATED_STRUCT( (SimFitParams, double, SIM_FIT_PARAMS_MEMBERS , (1)(1)(1)(1)(0.5)(0.8)(0)(0)(0) ) );

#define _ENUMERATED_STRUCT_BOOST_MEMBER(R,DUMMY,IDX,ELEM) T ELEM = b(p. ELEM, IDX);
#define ENUMERATED_STRUCT_BOOST_MEMBERS(DEF) BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_BOOST_MEMBER, , DEF)

class FitSim{
  int Lt;
  double tsep_pipi;
  double Ascale;
  double Cscale;

public:
  typedef SimFitParams Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoord GeneralizedCoordinate;

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    ENUMERATED_STRUCT_BOOST_MEMBERS(SIM_FIT_PARAMS_MEMBERS);

    switch(x.type){
    case SimFitType::PiPi2pt:
      return 
	cpipi0 * cpipi0 * Ascale * ( exp(-E0*x.t) + exp(-E0*(Lt-2.*tsep_pipi-x.t)) ) +
	cpipi1 * cpipi1 * Ascale * ( exp(-E1*x.t) + exp(-E1*(Lt-2.*tsep_pipi-x.t)) ) +
	Cpipipipi * Cscale;
    case SimFitType::PiPiToSigma:
      return 
	cpipi0 * csigma0 * Ascale * ( exp(-E0*x.t) + exp(-E0*(Lt-tsep_pipi-x.t)) ) +
	cpipi1 * csigma1 * Ascale * ( exp(-E1*x.t) + exp(-E1*(Lt-tsep_pipi-x.t)) ) +
	Cpipisigma * Cscale;
    case SimFitType::Sigma2pt:
      return 
	csigma0 * csigma0 * Ascale * ( exp(-E0*x.t) + exp(-E0*(Lt-x.t)) ) +
	csigma1 * csigma1 * Ascale * ( exp(-E1*x.t) + exp(-E1*(Lt-x.t)) ) +
	Csigmasigma * Cscale;
    default:
      assert(0);
    }
  }
  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d;
    for(int i=0;i<9;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSim(const int Lt, const double tsep_pipi, const double Ascale=1e13, const double Cscale=1e13): Lt(Lt), tsep_pipi(tsep_pipi), Ascale(Ascale), Cscale(Cscale){}

  static inline int Nparams(){ return 9; }
};

#define _ENUMERATED_STRUCT_JACKPRINT_MEMBER(R,DUMMY,IDX,ELEM) os << BOOST_PP_STRINGIZE(ELEM) << " = (" << cen. ELEM << " +- " << err. ELEM << ")\n";
#define ENUMERATED_STRUCT_JACKPRINT_MEMBERS(DEF) BOOST_PP_SEQ_FOR_EACH_I(_ENUMERATED_STRUCT_JACKPRINT_MEMBER, , DEF)

template<>
struct pipiParamsPrinter<FitSim>: public distributionPrinter<jackknifeDistribution< FitSim::Params> >{
  void print(std::ostream &os, const jackknifeDistribution< FitSim::Params> &dist) const{
    FitSim::Params cen = dist.best();
    FitSim::Params err = dist.standardError();
    ENUMERATED_STRUCT_JACKPRINT_MEMBERS(SIM_FIT_PARAMS_MEMBERS);
  }
};

#endif

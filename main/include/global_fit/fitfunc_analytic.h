#ifndef _GLOBALFIT_FITFUNC_ANALYTIC_H_
#define _GLOBALFIT_FITFUNC_ANALYTIC_H_

#define ANALYTIC_PQ_FIT_MPI2_PARAMS (AnalyticPQfitMpi2Params, double, (Zl)(Zh)(Ra)(c0)(ca)(cl)(cvl)(ch), (0.)(0.)(0.)(0.)(0.)(0.)(0.)(0.) )
DEF_ENUMERATED_STRUCT( ANALYTIC_PQ_FIT_MPI2_PARAMS );

class AnalyticPQfitMpi2{
public:
  typedef double ValueType;
  typedef AnalyticPQfitMpi2Params ParameterType;
  typedef AnalyticPQfitMpi2Params ValueDerivativeType; //derivative wrt parameters
  typedef DataParams GeneralizedCoordinate;

  bool a2_term;
  double mh0; //expansion heavy mass in normalization of reference ensemble

  AnalyticPQfitMpi2(bool a2_term, double mh0): a2_term(a2_term), mh0(mh0){}

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    ENUMERATED_STRUCT_BOOST_MEMBERS(ANALYTIC_PQ_FIT_MPI2_PARAMS);

    T mvl_r = Zl*Ra*(x.mx + x.my)/2.0;
    T ml_r = Zl*Ra*x.ml;
    T mh_r = (Zh*Ra*x.mh - mh0);
    T a2cpt = a2_term ? ca/Ra/Ra : T(0.);

    return ( c0*(1. + a2cpt) + cvl*mvl_r + cl*ml_r + ch*mh_r )/Ra/Ra;
  }
  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d;
    for(int i=0;i<8;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }
    
  inline int Nparams() const{ return 8; }

  ParameterType guess(){ return ParameterType(1,1,1,1,1,1,1,1); }
};


#define ANALYTIC_PQ_FIT_MK2_PARAMS (AnalyticPQfitMK2Params, double, (Zl)(Zh)(Ra)(c0)(ca)(cl)(cvl)(ch)(cvh), (0.)(0.)(0.)(0.)(0.)(0.)(0.)(0.)(0.) )
DEF_ENUMERATED_STRUCT( ANALYTIC_PQ_FIT_MK2_PARAMS );

class AnalyticPQfitMK2{
public:
  typedef double ValueType;
  typedef AnalyticPQfitMK2Params ParameterType;
  typedef AnalyticPQfitMK2Params ValueDerivativeType; //derivative wrt parameters
  typedef DataParams GeneralizedCoordinate;

  bool a2_term;
  double mh0; //expansion heavy mass in normalization of reference ensemble

  AnalyticPQfitMK2(bool a2_term, double mh0): a2_term(a2_term), mh0(mh0){}

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    ENUMERATED_STRUCT_BOOST_MEMBERS(ANALYTIC_PQ_FIT_MK2_PARAMS);

    T mvl_r = Zl*Ra*x.mx;
    T mvh_r = (Zh*Ra*x.my-mh0);
    T ml_r = Zl*Ra*x.ml;
    T mh_r = (Zh*Ra*x.mh-mh0);
    T a2cpt = a2_term ? ca/Ra/Ra : T(0.);

    return ( c0*(1. + a2cpt) + cvl*mvl_r + cl*ml_r + cvh*mvh_r + ch*mh_r )/Ra/Ra;
  }
  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d;
    for(int i=0;i<9;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }
    
  inline int Nparams() const{ return 9; }

  ParameterType guess(){ return ParameterType(1,1,1,1,1,1,1,1,1); }
};




#define ANALYTIC_PQ_FIT_MOMEGA_PARAMS (AnalyticPQfitMOmegaParams, double, (Zl)(Zh)(Ra)(c0)(ca)(cl)(ch)(cvh), (0.)(0.)(0.)(0.)(0.)(0.)(0.)(0.) )
DEF_ENUMERATED_STRUCT( ANALYTIC_PQ_FIT_MOMEGA_PARAMS );

class AnalyticPQfitMOmega{
public:
  typedef double ValueType;
  typedef AnalyticPQfitMOmegaParams ParameterType;
  typedef AnalyticPQfitMOmegaParams ValueDerivativeType; //derivative wrt parameters
  typedef DataParams GeneralizedCoordinate;

  bool a2_term;
  double mh0; //expansion heavy mass in normalization of reference ensemble

  AnalyticPQfitMOmega(bool a2_term, double mh0): a2_term(a2_term), mh0(mh0){}

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    ENUMERATED_STRUCT_BOOST_MEMBERS(ANALYTIC_PQ_FIT_MOMEGA_PARAMS);

    T mvh_r = (Zh*Ra*x.my-mh0);
    T ml_r = Zl*Ra*x.ml;
    T mh_r = (Zh*Ra*x.mh-mh0);
    T a2cpt = a2_term ? ca/Ra/Ra : T(0.);

    return ( c0*(1. + a2cpt) + cl*ml_r + cvh*mvh_r + ch*mh_r )/Ra;
  }
  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d;
    for(int i=0;i<8;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }
    
  inline int Nparams() const{ return 8; }

  ParameterType guess(){ return ParameterType(1,1,1,1,1,1,1,1); }
};


DEF_ENUMERATED_STRUCT( (						\
			AnalyticGlobalFitParams_base, double,		\
			(mpi2_c0)(mpi2_ca)(mpi2_cl)(mpi2_cvl)(mpi2_ch)	\
			(mK2_c0)(mK2_ca)(mK2_cl)(mK2_cvl)(mK2_ch)(mK2_cvh)	\
			(mOmega_c0)(mOmega_ca)(mOmega_cl)(mOmega_ch)(mOmega_cvh) \
			,						\
			(0.1*0.1)(1e-3)(1e-3)(1e-3)(1e-3)		\
			(0.3*0.3)(1e-3)(1e-3)(1e-3)(1e-3)(1e-3)		\
			(0.8*0.8)(1e-3)(1e-3)(1e-3)(1e-3)		\
			)						\
		       );

typedef GlobalFitParams<AnalyticGlobalFitParams_base> AnalyticGlobalFitParams; //includes Zl, Zh, Ra vectors

GENERATE_PARSER( AnalyticGlobalFitParams, (std::vector<double>, Zl)(std::vector<double>, Zh)(std::vector<double>, Ra)(AnalyticGlobalFitParams_base, params) );

//Bindings for parameters other than scaling parameters
DEF_GLOBALFIT_BINDING(AnalyticGlobalFitBindings_mpi2, AnalyticPQfitMpi2Params, AnalyticGlobalFitParams_base, (mpi2_c0,c0)(mpi2_ca,ca)(mpi2_cl,cl)(mpi2_cvl,cvl)(mpi2_ch,ch) );
DEF_GLOBALFIT_BINDING(AnalyticGlobalFitBindings_mK2, AnalyticPQfitMK2Params, AnalyticGlobalFitParams_base, (mK2_c0,c0)(mK2_ca,ca)(mK2_cl,cl)(mK2_cvl,cvl)(mK2_ch,ch)(mK2_cvh,cvh) );
DEF_GLOBALFIT_BINDING(AnalyticGlobalFitBindings_mOmega, AnalyticPQfitMOmegaParams, AnalyticGlobalFitParams_base, (mOmega_c0,c0)(mOmega_ca,ca)(mOmega_cl,cl)(mOmega_ch,ch)(mOmega_cvh,cvh) );

struct AnalyticGlobalFitPolicies{
  AnalyticPQfitMpi2 fit_mpi2;
  AnalyticPQfitMK2 fit_mK2;
  AnalyticPQfitMOmega fit_mOmega;
  
  int nlat;

  typedef AnalyticGlobalFitParams ParameterType;

  typedef AnalyticGlobalFitBindings_mpi2 Bindings_mpi2;
  typedef AnalyticGlobalFitBindings_mK2 Bindings_mK2;
  typedef AnalyticGlobalFitBindings_mOmega Bindings_mOmega;

  struct fitFuncSetupType{
    int nlat;
    bool a2_term;
    double mh0;
  };

  inline int Nparams() const{ ParameterType p(nlat); return p.size(); }

  inline ParameterType guess() const{ ParameterType p(nlat); return p; }

  AnalyticGlobalFitPolicies(const fitFuncSetupType &params): nlat(params.nlat), 
							     fit_mpi2(params.a2_term, params.mh0),
							     fit_mK2(params.a2_term, params.mh0),
							     fit_mOmega(params.a2_term, params.mh0){}
};

typedef GlobalFit<AnalyticGlobalFitPolicies> AnalyticGlobalFit;


#endif

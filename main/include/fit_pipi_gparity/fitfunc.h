#ifndef _PIPI_FITFUNC_H_
#define _PIPI_FITFUNC_H_

#include<parser.h>

class FitCoshPlusConstant{
public:
  struct Params{
    double A;
    double E;
    double C;
    Params(): A(1), E(0.3), C(0){}
    Params(const double _A, const double _E, const double _C): A(_A), E(_E), C(_C){}
    inline int size() const{ return 3; }
    inline std::string print() const{ std::ostringstream os; os << "A=" << A << " E=" << E << " C=" << C; return os.str(); }
    inline void zero(){ A=E=C=0.; }
    inline double &operator()(const int i){
      switch(i){
      case 0:
	return A;
      case 1:
	return E;
      case 2:
	return C;
      }
    }
    inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<Params*>(this)->operator()(i) ); }

    typedef FitCoshPlusConstant::Params ET_tag;
    template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,FitCoshPlusConstant::Params>::value, int>::type = 0>
											      Params(U&& expr);

    inline double pipiEnergy() const{ return E; }
    inline double constant() const{ return C; }
  };

  struct Derivs{
    double dA;
    double dE;
    double dC;
    Derivs(){}
    Derivs(const double _dA, const double _dE, const double _dC): dA(_dA), dE(_dE), dC(_dC){}
    inline int size() const{ return 3; }
    inline std::string print() const{ std::ostringstream os; os << "dA=" << dA << " dE=" << dE << " dC=" << dC; return os.str(); }
    inline void resize(const int n){};
    inline void zero(){ dA=dE=dC=0.; }
    inline double &operator()(const int i){
      switch(i){
      case 0:
	return dA;
      case 1:
	return dE;
      case 2:
	return dC;
      }
    }
    inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<Derivs*>(this)->operator()(i) ); }      
  };
private:
    
  const double Lt;
  const double tsep_pipi;
  const double Ascale;
  const double Cscale;
public:
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Derivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitCoshPlusConstant(const double _Lt, const double _tsep_pipi, const double _Ascale = 1e13, const double _Cscale = 1e13):
    Lt(_Lt), tsep_pipi(_tsep_pipi), Ascale(_Ascale), Cscale(_Cscale){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
    return p.A * Ascale * ( exp(-p.E*t) + exp(-p.E*(Lt-2.*tsep_pipi-t)) ) + p.C * Cscale;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    yderivs.dA = Ascale * ( exp(-p.E*t) + exp(-p.E*(Lt-2.*tsep_pipi-t)) );
    yderivs.dE = p.A * Ascale * ( (-t)*exp(-p.E*t) - (Lt-2.*tsep_pipi-t)*exp(-p.E*(Lt-2.*tsep_pipi-t)) );
    yderivs.dC = Cscale;    
    return yderivs;
  }

  inline int Nparams() const{ return 3; }  
};

  
template<>
struct getElem<FitCoshPlusConstant::Params>{
  static inline auto elem(const FitCoshPlusConstant::Params &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const FitCoshPlusConstant::Params &v){ return 0; }
};

template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, FitCoshPlusConstant::Params::ET_tag>::value && !std::is_same<U,FitCoshPlusConstant::Params>::value, int>::type>
FitCoshPlusConstant::Params::Params(U&& expr){
  this->A = expr[0];
  this->E = expr[1];
  this->C = expr[2];  
}

GENERATE_PARSER_GM(FitCoshPlusConstant::Params, FitCoshPlusConstant_Params_grammar, (double, A)(double, E)(double, C) )

template<>
struct printStats<jackknifeDistribution<FitCoshPlusConstant::Params> >{
  inline static std::string centralValue(const jackknifeDistribution<FitCoshPlusConstant::Params> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.A << ", " << best.E << ", " << best.C << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<FitCoshPlusConstant::Params> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.A << ", " << err.E << ", " << err.C << ")";
    return os.str();
  }
};

template<typename FitFunc>
struct pipiParamsPrinter{};

template<>
struct pipiParamsPrinter<FitCoshPlusConstant>: public distributionPrinter<jackknifeDistribution<FitCoshPlusConstant::Params> >{
  void print(std::ostream &os, const jackknifeDistribution<FitCoshPlusConstant::Params> &dist) const{
    FitCoshPlusConstant::Params cen = dist.best();
    FitCoshPlusConstant::Params err = dist.standardError();
    os << "A = (" << cen.A << " +- " << err.A << ") E = (" << cen.E << " +- " << err.E << ") C = (" << cen.C << " +- " << err.C << ")";
  }
};



struct WindowBoundedMapping{
  static inline double map(const double pint, const double min, const double max){
    return min + (sin(pint)+1)*(max-min)/2.;
  }
  //dmap/dpint
  static inline double deriv(const double pint, const double min, const double max){
    return cos(pint)*(max-min)/2.;
  }
};
struct MinBoundedMapping{
  static inline double map(const double pint, const double min){
    return min -1.+ sqrt( pint*pint + 1 );
  }
  //dmap/dpint
  static inline double deriv(const double pint, const double min){
    return pint/sqrt( pint*pint + 1 );
  }
};
struct MaxBoundedMapping{
  static inline double map(const double pint, const double max){
    return max +1.- sqrt( pint*pint + 1 );
  }
  //dmap/dpint
  static inline double deriv(const double pint, const double max){
    return -pint/sqrt( pint*pint + 1 );
  }
};


//#define PIPI_2EXP_E1BOUNDED

class FitCoshPlusConstantDoubleExp{
public:
  struct Params{
    double A0;
    double E0;
    double A1;
    double E1;
    double C;
    Params(): A0(1), E0(0.3), A1(1), E1(1), C(0){}
    Params(const double _A0, const double _E0, const double _A1, const double _E1, const double _C): A0(_A0), E0(_E0), A1(_A1), E1(_E1), C(_C){}
    inline int size() const{ return 5; }
    inline std::string print() const{ std::ostringstream os; os << "A0=" << A0 << " E0=" << E0 << "A1=" << A1 << " E1=" << E1 << " C=" << C; return os.str(); }

    inline void zero(){ A0=E0=A1=E1=C=0.; }
    inline double &operator()(const int i){
      switch(i){
      case 0:
	return A0;
      case 1:
	return E0;
      case 2:
	return A1;
      case 3:
	return E1;
      case 4:
	return C;
      }
    }
    inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<Params*>(this)->operator()(i) ); }

    typedef FitCoshPlusConstantDoubleExp::Params ET_tag;
    template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,FitCoshPlusConstantDoubleExp::Params>::value, int>::type = 0>
											      Params(U&& expr);
    inline double pipiEnergy() const{ return E0; }
    inline double constant() const{ return C; }
  };

  struct Derivs{
    double dA0;
    double dE0;
    double dA1;
    double dE1;    
    double dC;
    Derivs(){}
    Derivs(const double _dA0, const double _dE0, const double _dA1, const double _dE1, const double _dC): dA0(_dA0), dE0(_dE0), dA1(_dA1), dE1(_dE1), dC(_dC){}
    inline int size() const{ return 5; }
    inline std::string print() const{ std::ostringstream os; os << "dA0=" << dA0 << " dE0=" << dE0 << "dA1=" << dA1 << " dE1=" << dE1 << " dC=" << dC; return os.str(); }
    inline void resize(const int n){};
    inline void zero(){ dA0=dE0=dA1=dE1=dC=0.; }
    inline double &operator()(const int i){
      switch(i){
      case 0:
	return dA0;
      case 1:
	return dE0;
      case 2:
	return dA1;
      case 3:
	return dE1;	
      case 4:
	return dC;
      }
    }
    inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<Derivs*>(this)->operator()(i) ); }      
  };
private:
    
  const double Lt;
  const double tsep_pipi;
  const double Ascale;
  const double Cscale;
public:
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Derivs ValueDerivativeType; //derivative wrt parameters
  typedef double GeneralizedCoordinate; //time coord

  FitCoshPlusConstantDoubleExp(const double _Lt, const double _tsep_pipi, const double _Ascale = 1e13, const double _Cscale = 1e13):
    Lt(_Lt), tsep_pipi(_tsep_pipi), Ascale(_Ascale), Cscale(_Cscale){}
  
  //Params are A, m  
  ValueType value(const GeneralizedCoordinate &t, const ParameterType &p) const{
#ifdef PIPI_2EXP_E1BOUNDED
    double E1_true = MinBoundedMapping::map(p.E1, p.E0); //force E1 to be > E0
#else 
    double E1_true = p.E1;
#endif

    return
      p.A0 * Ascale * ( exp(-p.E0*t) + exp(-p.E0*(Lt-2.*tsep_pipi-t)) ) +
      p.A1 * Ascale * ( exp(-E1_true*t) + exp(-E1_true*(Lt-2.*tsep_pipi-t)) ) +
      p.C * Cscale;
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &t, const ParameterType &p) const{
#ifdef PIPI_2EXP_E1BOUNDED
    double E1_true = MinBoundedMapping::map(p.E1, p.E0);
    double dE1_true_by_dE1 = MinBoundedMapping::deriv(p.E1, p.E0);
#else
    double E1_true = p.E1;
    double dE1_true_by_dE1 = 1;
#endif

    ValueDerivativeType yderivs;
    yderivs.dA0 = Ascale * ( exp(-p.E0*t) + exp(-p.E0*(Lt-2.*tsep_pipi-t)) );
    yderivs.dE0 = p.A0 * Ascale * ( (-t)*exp(-p.E0*t) - (Lt-2.*tsep_pipi-t)*exp(-p.E0*(Lt-2.*tsep_pipi-t)) );

    yderivs.dA1 = Ascale * ( exp(-E1_true*t) + exp(-E1_true*(Lt-2.*tsep_pipi-t)) );
    yderivs.dE1 = p.A1 * Ascale * ( (-t)*exp(-E1_true*t) - (Lt-2.*tsep_pipi-t)*exp(-E1_true*(Lt-2.*tsep_pipi-t)) ) * dE1_true_by_dE1;
	    
    yderivs.dC = Cscale;    
    return yderivs;
  }

  inline int Nparams() const{ return 5; }
};

  
template<>
struct getElem<FitCoshPlusConstantDoubleExp::Params>{
  static inline auto elem(const FitCoshPlusConstantDoubleExp::Params &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const FitCoshPlusConstantDoubleExp::Params &v){ return 0; }
};

template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, FitCoshPlusConstantDoubleExp::Params::ET_tag>::value && !std::is_same<U,FitCoshPlusConstantDoubleExp::Params>::value, int>::type>
FitCoshPlusConstantDoubleExp::Params::Params(U&& expr){
  this->A0 = expr[0];
  this->E0 = expr[1];
  this->A1 = expr[2];
  this->E1 = expr[3];
  this->C = expr[4];  
}

GENERATE_PARSER_GM(FitCoshPlusConstantDoubleExp::Params, FitCoshPlusConstantDoubleExp_Params_grammar, (double, A0)(double, E0)(double, A1)(double, E1)(double, C) )

template<>
struct printStats<jackknifeDistribution<FitCoshPlusConstantDoubleExp::Params> >{
  inline static std::string centralValue(const jackknifeDistribution<FitCoshPlusConstantDoubleExp::Params> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.A0 << ", " << best.E0 << ", " << best.A1 << ", " << best.E1 << ", " << best.C << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<FitCoshPlusConstantDoubleExp::Params> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.A0 << ", " << err.E0 << ", " << err.A1 << ", " << err.E1 << ", " << err.C << ")";
    return os.str();
  }
};

template<>
struct pipiParamsPrinter<FitCoshPlusConstantDoubleExp>: public distributionPrinter<jackknifeDistribution<FitCoshPlusConstantDoubleExp::Params> >{
  void print(std::ostream &os, const jackknifeDistribution<FitCoshPlusConstantDoubleExp::Params> &dist) const{
    int nsample = dist.size();
    jackknifeDistribution<double> E1(dist.size());
#ifdef PIPI_2EXP_E1BOUNDED
    for(int s=0;s<nsample;s++) E1.sample(s) = MinBoundedMapping::map(dist.sample(s).E1, dist.sample(s).E0);
#else
    for(int s=0;s<nsample;s++) E1.sample(s) = dist.sample(s).E1;
#endif

    FitCoshPlusConstantDoubleExp::Params cen = dist.best();
    FitCoshPlusConstantDoubleExp::Params err = dist.standardError();

    double cen_E1 = E1.best();
    double err_E1 = E1.standardError();
    
    os << "A0 = (" << cen.A0 << " +- " << err.A0 << ") E0 = (" << cen.E0 << " +- " << err.E0 << ") A1 = (" << cen.A1 << " +- " << err.A1 << ") E1 = (" << cen_E1 << " +- " << err_E1 << ") C = (" << cen.C << " +- " << err.C << ")";
  }
};


#endif

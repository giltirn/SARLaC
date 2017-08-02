#ifndef _PIPI_FITFUNC_H_
#define _PIPI_FITFUNC_H_

#include<parser.h>

class FitCoshPlusConstant{
public:
  struct Params{
    double A;
    double E;
    double C;
    Params(){}
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

struct pipiParamsPrinter: public distributionPrinter<jackknifeDistribution<FitCoshPlusConstant::Params> >{
  void print(std::ostream &os, const jackknifeDistribution<FitCoshPlusConstant::Params> &dist) const{
    FitCoshPlusConstant::Params cen = dist.best();
    FitCoshPlusConstant::Params err = dist.standardError();
    os << "A = (" << cen.A << " +- " << err.A << ") E = (" << cen.E << " +- " << err.E << ") C = (" << cen.C << " +- " << err.C << ")";
  }
};
#endif

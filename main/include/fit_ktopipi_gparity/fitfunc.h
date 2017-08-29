#ifndef _FIT_KTOPIPI_GPARITY_FITFUNC_H
#define _FIT_KTOPIPI_GPARITY_FITFUNC_H

#include<parser.h>

struct amplitudeDataCoord{
  double t;
  int tsep_k_pi;
  inline amplitudeDataCoord(double _t, int _tsep_k_pi): t(_t), tsep_k_pi(_tsep_k_pi){}
};
std::ostream & operator<<(std::ostream &os, const amplitudeDataCoord &c){
  os << "(t=" << c.t << ", tsep_k_pi=" << c.tsep_k_pi << ")";
  return os;
}

class FitKtoPiPi{
public:
  struct Params{
    double AK;
    double mK;
    double Apipi;
    double Epipi;
    double M;

    Params(): AK(1), mK(0.5), Apipi(1), Epipi(0.5), M(1){}
    Params(const double _AK, const double _mK, const double _Apipi, const double _Epipi, const double _M): AK(_AK), mK(_mK), Apipi(_Apipi), Epipi(_Epipi), M(_M){}

    inline int size() const{ return 5; }
    inline std::string print() const{ std::ostringstream os; os << "AK=" << AK << " mK=" << mK << " Apipi=" << Apipi << " Epipi=" << Epipi << " M=" << M; return os.str(); }
    inline void zero(){ AK=mK=Apipi=Epipi=M=0.; }
    inline double &operator()(const int i){
      switch(i){
      case 0:
	return AK;
      case 1:
	return mK;
      case 2:
	return Apipi;
      case 3:
	return Epipi;
      case 4:
	return M;
      }
    }  
    inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<Params*>(this)->operator()(i) ); }

    typedef FitKtoPiPi::Params ET_tag;
    template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,FitKtoPiPi::Params>::value, int>::type = 0>
    Params(U&& expr);
  };
public:
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType; //derivative wrt parameters
  typedef amplitudeDataCoord GeneralizedCoordinate; 

  ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    return p.AK*p.Apipi*p.M*exp(-p.Epipi * c.tsep_k_pi)*exp( -(p.mK - p.Epipi)*c.t )/sqrt(2.);
  }
  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p) const{
    ValueDerivativeType yderivs;
    double v = value(c,p);    
    yderivs.AK = v/p.AK;
    yderivs.Apipi = v/p.Apipi;
    yderivs.M = v/p.M;
    yderivs.Epipi = -c.tsep_k_pi*v +c.t*v;
    yderivs.mK = -c.t*v;
    return yderivs;
  }

  inline int Nparams() const{ return 5; }
};

  
template<>
struct getElem<FitKtoPiPi::Params>{
  static inline auto elem(const FitKtoPiPi::Params &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const FitKtoPiPi::Params &v){ return 0; }
};

template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, FitKtoPiPi::Params::ET_tag>::value && !std::is_same<U,FitKtoPiPi::Params>::value, int>::type>
FitKtoPiPi::Params::Params(U&& expr){
  this->AK = expr[0];
  this->mK = expr[1];
  this->Apipi = expr[2];
  this->Epipi = expr[3];
  this->M = expr[4];
}

GENERATE_PARSER_GM(FitKtoPiPi::Params, FitKtoPiPi_Params_grammar, (double, AK)(double, mK)(double, Apipi)(double, Epipi)(double, M) )

template<>
struct printStats<jackknifeDistribution<FitKtoPiPi::Params> >{
  inline static std::string centralValue(const jackknifeDistribution<FitKtoPiPi::Params> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " << best.Apipi << ", " << best.Epipi << ", " << best.M << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<FitKtoPiPi::Params> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " << err.Apipi << ", " << err.Epipi << ", " << err.M << ")";
    return os.str();
  }
};

template<typename FitFunc>
struct ktopipiParamsPrinter{};

template<>
struct ktopipiParamsPrinter<FitKtoPiPi>: public distributionPrinter<jackknifeDistribution<FitKtoPiPi::Params> >{
  void print(std::ostream &os, const jackknifeDistribution<FitKtoPiPi::Params> &dist) const{
    FitKtoPiPi::Params cen = dist.best();
    FitKtoPiPi::Params err = dist.standardError();
    os << "AK = (" << cen.AK << " +- " << err.AK << ") mK = (" << cen.mK << " +- " << err.mK << ")"
       << "Apipi = (" << cen.Apipi << " +- " << err.Apipi << ") Epipi = (" << cen.Epipi << " +- " << err.Epipi << ")"      
       << "M = (" << cen.M << " +- " << err.M << ")";
  }
};


#endif

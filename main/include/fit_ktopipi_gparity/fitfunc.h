#ifndef _FIT_KTOPIPI_GPARITY_FITFUNC_H
#define _FIT_KTOPIPI_GPARITY_FITFUNC_H

#include<parser.h>

struct amplitudeDataCoord{
  double t; //K->op separation
  int tsep_k_pi;
  inline amplitudeDataCoord(){}
  inline amplitudeDataCoord(double _t, int _tsep_k_pi): t(_t), tsep_k_pi(_tsep_k_pi){}
};
std::ostream & operator<<(std::ostream &os, const amplitudeDataCoord &c){
  os << "(t=" << c.t << ", tsep_k_pi=" << c.tsep_k_pi << ")";
  return os;
}
#ifdef HAVE_HDF5

void write(HDF5writer &writer, const amplitudeDataCoord &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.t, "t");
  write(writer,value.tsep_k_pi, "tsep_k_pi");
  writer.leave();
}
void read(HDF5reader &reader, amplitudeDataCoord &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.t, "t");
  read(reader,value.tsep_k_pi, "tsep_k_pi");
  reader.leave();
}

#endif



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

  static ValueType value(const GeneralizedCoordinate &c, const ParameterType &p){
    return p.AK*p.Apipi*p.M*exp(-p.Epipi * c.tsep_k_pi)*exp( -(p.mK - p.Epipi)*c.t )/sqrt(2.);
  }
  static ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p){
    ValueDerivativeType yderivs;
    double v = value(c,p);    
    yderivs.AK = v/p.AK;
    yderivs.Apipi = v/p.Apipi;
    yderivs.M = v/p.M;
    yderivs.Epipi = -c.tsep_k_pi*v +c.t*v;
    yderivs.mK = -c.t*v;
    return yderivs;
  }

  static inline int Nparams(){ return 5; }
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





//Simultaneous fit
struct amplitudeDataCoordSim{
  double t; //K->op separation
  int tsep_k_pi;
  int idx;
  inline amplitudeDataCoordSim(){}
  inline amplitudeDataCoordSim(double _t, int _tsep_k_pi, int _idx): t(_t), tsep_k_pi(_tsep_k_pi), idx(_idx){}
  inline amplitudeDataCoordSim(const amplitudeDataCoord &c, const int _idx): t(c.t), tsep_k_pi(c.tsep_k_pi), idx(_idx){}
};
std::ostream & operator<<(std::ostream &os, const amplitudeDataCoordSim &c){
  os << "(t=" << c.t << ", tsep_k_pi=" << c.tsep_k_pi << ", idx=" << c.idx << ")";
  return os;
}
#ifdef HAVE_HDF5

void write(HDF5writer &writer, const amplitudeDataCoordSim &value, const std::string &tag){
  writer.enter(tag);
  write(writer,value.t, "t");
  write(writer,value.tsep_k_pi, "tsep_k_pi");
  write(writer,value.idx, "idx");
  writer.leave();
}
void read(HDF5reader &reader, amplitudeDataCoordSim &value, const std::string &tag){
  reader.enter(tag);
  read(reader,value.t, "t");
  read(reader,value.tsep_k_pi, "tsep_k_pi");
  read(reader,value.idx, "idx");
  reader.leave();
}

#endif


class FitKtoPiPiSim{
public:
  struct Params{
    double AK;
    double mK;
    double Apipi;
    double Epipi;
    std::vector<double> M;

    Params(): AK(1), mK(0.5), Apipi(1), Epipi(0.5), M(10,1){}
    Params(const double _AK, const double _mK, const double _Apipi, const double _Epipi, const double _M): AK(_AK), mK(_mK), Apipi(_Apipi), Epipi(_Epipi), M(10,_M){}
    Params(const Params &r) = default;
    Params(Params &&r) = default;

    Params& operator=(const Params &r) = default;
    Params& operator=(Params &&r) = default;
    
    inline int size() const{ return 14; }

    inline std::string print() const{
      std::ostringstream os; os << "AK=" << AK << " mK=" << mK << " Apipi=" << Apipi << " Epipi=" << Epipi;
      for(int i=0;i<10;i++) os << " M[" << i << "]=" << M[i];
      return os.str();
    }
    inline void zero(){
      AK=mK=Apipi=Epipi=0.;
      for(int i=0;i<10;i++) M[i] = 0.;
    }
    inline double &operator()(const int i){
      if(i<4){
	switch(i){
	case 0:
	  return AK;
	case 1:
	  return mK;
	case 2:
	  return Apipi;
	case 3:
	  return Epipi;
	}
      }else return M[i-4];
    }
    inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<Params*>(this)->operator()(i) ); }

    typedef FitKtoPiPiSim::Params ET_tag;
    template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,FitKtoPiPiSim::Params>::value, int>::type = 0>
    Params(U&& expr);
  };
public:
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType; //derivative wrt parameters
  typedef amplitudeDataCoordSim GeneralizedCoordinate; 

  static ValueType value(const GeneralizedCoordinate &c, const ParameterType &p){
    return p.AK*p.Apipi*p.M[c.idx]*exp(-p.Epipi * c.tsep_k_pi)*exp( -(p.mK - p.Epipi)*c.t )/sqrt(2.);
  }
  static ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p){
    ValueDerivativeType yderivs;  yderivs.zero();
    double v = value(c,p);
    yderivs.AK = v/p.AK;
    yderivs.Apipi = v/p.Apipi;
    yderivs.Epipi = -c.tsep_k_pi*v +c.t*v;
    yderivs.mK = -c.t*v;
    yderivs.M[c.idx] = v/p.M[c.idx];
    return yderivs;
  }

  static inline int Nparams(){ return 14; }
};

  
template<>
struct getElem<FitKtoPiPiSim::Params>{
  static inline auto elem(const FitKtoPiPiSim::Params &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const FitKtoPiPiSim::Params &v){ return 0; }
};

template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, FitKtoPiPiSim::Params::ET_tag>::value && !std::is_same<U,FitKtoPiPiSim::Params>::value, int>::type>
FitKtoPiPiSim::Params::Params(U&& expr): M(10){
  this->AK = expr[0];
  this->mK = expr[1];
  this->Apipi = expr[2];
  this->Epipi = expr[3];
  for(int i=4;i<14;i++) this->M[i-4] = expr[i];
}

GENERATE_PARSER_GM(FitKtoPiPiSim::Params, FitKtoPiPiSim_Params_grammar, (double, AK)(double, mK)(double, Apipi)(double, Epipi)(std::vector<double>, M) )

template<>
struct printStats<jackknifeDistribution<FitKtoPiPiSim::Params> >{
  inline static std::string centralValue(const jackknifeDistribution<FitKtoPiPiSim::Params> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " << best.Apipi << ", " << best.Epipi;
    for(int i=0;i<10;i++) os << ", " << best.M[i];
    os << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<FitKtoPiPiSim::Params> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " << err.Apipi << ", " << err.Epipi;
    for(int i=0;i<10;i++) os << ", " << err.M[i];
    os << ")";
    return os.str();
  }
};


template<>
struct ktopipiParamsPrinter<FitKtoPiPiSim>: public distributionPrinter<jackknifeDistribution<FitKtoPiPiSim::Params> >{
  void print(std::ostream &os, const jackknifeDistribution<FitKtoPiPiSim::Params> &dist) const{
    FitKtoPiPiSim::Params cen = dist.best();
    FitKtoPiPiSim::Params err = dist.standardError();
    os << "AK = (" << cen.AK << " +- " << err.AK << ") mK = (" << cen.mK << " +- " << err.mK << ") "
       << "Apipi = (" << cen.Apipi << " +- " << err.Apipi << ") Epipi = (" << cen.Epipi << " +- " << err.Epipi << ")";
    for(int i=0;i<10;i++) os << " M[" << i << "] = (" << cen.M[i] << " +- " << err.M[i] << ")";
  }
};





#endif

#ifndef _FIT_KTOPIPI_GPARITY_FITFUNC_H
#define _FIT_KTOPIPI_GPARITY_FITFUNC_H

#include<parser.h>
#include<containers/enumerated_struct.h>

struct amplitudeDataCoord{
  double t; //K->op separation
  int tsep_k_pi;
  inline amplitudeDataCoord(){}
  inline amplitudeDataCoord(double _t, int _tsep_k_pi): t(_t), tsep_k_pi(_tsep_k_pi){}

  inline bool operator==(const amplitudeDataCoord &r) const{ return t==r.t && tsep_k_pi==r.tsep_k_pi; }
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
#define FIT_KTOPIPI_PARAMS (Params, double, (AK)(mK)(Apipi)(Epipi)(M), (1.0)(0.5)(1.0)(0.5)(1.0) )
  DEF_ENUMERATED_STRUCT_MEMBER(FitKtoPiPi, FIT_KTOPIPI_PARAMS);

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
DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(FitKtoPiPi, FIT_KTOPIPI_PARAMS);

template<>
struct CPSfit::printStats<jackknifeDistribution<FitKtoPiPi::Params> >{
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


//Same as above but with extra constant
class FitKtoPiPiWithConstant{
public:
#define FIT_KTOPIPI_WITH_CONSTANT_PARAMS (Params, double, (AK)(mK)(Apipi)(Epipi)(M)(C), (1.0)(0.5)(1.0)(0.5)(1.0)(0.) )
  DEF_ENUMERATED_STRUCT_MEMBER(FitKtoPiPiWithConstant, FIT_KTOPIPI_WITH_CONSTANT_PARAMS);

public:
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType; //derivative wrt parameters
  typedef amplitudeDataCoord GeneralizedCoordinate; 

  static ValueType value(const GeneralizedCoordinate &c, const ParameterType &p){
    return p.AK*p.Apipi*p.M*exp(-p.Epipi * c.tsep_k_pi)*exp( -(p.mK - p.Epipi)*c.t )/sqrt(2.) + p.C;
  }
  static ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p){
    ValueDerivativeType yderivs;
    double v = value(c,p);    
    yderivs.AK = v/p.AK;
    yderivs.Apipi = v/p.Apipi;
    yderivs.M = v/p.M;
    yderivs.Epipi = -c.tsep_k_pi*v +c.t*v;
    yderivs.mK = -c.t*v;
    yderivs.C = 1.;
    return yderivs;
  }

  static inline int Nparams(){ return 6; }
};
DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(FitKtoPiPiWithConstant, FIT_KTOPIPI_WITH_CONSTANT_PARAMS);

template<>
struct CPSfit::printStats<jackknifeDistribution<FitKtoPiPiWithConstant::Params> >{
  inline static std::string centralValue(const jackknifeDistribution<FitKtoPiPiWithConstant::Params> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " << best.Apipi << ", " << best.Epipi << ", " << best.M << ", " << best.C << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<FitKtoPiPiWithConstant::Params> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " << err.Apipi << ", " << err.Epipi << ", " << err.M << ", " << err.C << ")";
    return os.str();
  }
};



//Same as above but with extra constant
class FitKtoPiPiTwoExp{
public:
#define FIT_KTOPIPI_TWO_EXP_PARAMS (Params, double, (AK)(mK)(A0)(E0)(A1)(E1)(M0)(M1), (1.0)(0.5)(1.0)(0.5)(1.0)(0.5)(1.0)(1.0) )
  DEF_ENUMERATED_STRUCT_MEMBER(FitKtoPiPiTwoExp, FIT_KTOPIPI_TWO_EXP_PARAMS);

public:
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType; //derivative wrt parameters
  typedef amplitudeDataCoord GeneralizedCoordinate; 

  static ValueType value(const GeneralizedCoordinate &c, const ParameterType &p){
    return 
      p.AK*p.A0*p.M0*exp(-p.E0 * c.tsep_k_pi)*exp( -(p.mK - p.E0)*c.t )/sqrt(2.) + 
      p.AK*p.A1*p.M1*exp(-p.E1 * c.tsep_k_pi)*exp( -(p.mK - p.E1)*c.t )/sqrt(2.);
  }
  static ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p){
    double v0 = p.AK*p.A0*p.M0*exp(-p.E0 * c.tsep_k_pi)*exp( -(p.mK - p.E0)*c.t )/sqrt(2.);
    double v1 = p.AK*p.A1*p.M1*exp(-p.E1 * c.tsep_k_pi)*exp( -(p.mK - p.E1)*c.t )/sqrt(2.);
    double v = v0 + v1;

    ValueDerivativeType yderivs;
    yderivs.AK = v/p.AK;
    yderivs.mK = -c.t*v;
    yderivs.A0 = v0/p.A0;
    yderivs.E0 = -c.tsep_k_pi*v0 +c.t*v0;
    yderivs.A1 = v1/p.A1;
    yderivs.E1 = -c.tsep_k_pi*v1 +c.t*v1;
    yderivs.M0 = v0/p.M0;
    yderivs.M1 = v1/p.M1;
    return yderivs;
  }

  static inline int Nparams(){ return 8; }

  static double getM0data(const double y, const GeneralizedCoordinate &c, const ParameterType &p){
    double v2 = p.AK*p.A1*p.M1*exp(-p.E1 * c.tsep_k_pi)*exp( -(p.mK - p.E1)*c.t )/sqrt(2.);
    double v1woM = p.AK*p.A0*exp(-p.E0 * c.tsep_k_pi)*exp( -(p.mK - p.E0)*c.t )/sqrt(2.);
    return (y - v2)/v1woM;
  }
};
DEF_ENUMERATED_STRUCT_MEMBER_EXTERNAL(FitKtoPiPiTwoExp, FIT_KTOPIPI_TWO_EXP_PARAMS);

template<>
struct CPSfit::printStats<jackknifeDistribution<FitKtoPiPiTwoExp::Params> >{
  inline static std::string centralValue(const jackknifeDistribution<FitKtoPiPiTwoExp::Params> &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " 
			      << best.A0 << ", " << best.E0 << ", " 
			      << best.A1 << ", " << best.E1 << ", "
			      << best.M0 << ", " << best.M1 << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<FitKtoPiPiTwoExp::Params> &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " 
			      << err.A0 << ", " << err.E0 << ", " 
			      << err.A1 << ", " << err.E1 << ", "
			      << err.M0 << ", " << err.M1 << ")";
    return os.str();
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



template<int N>
struct SimFitParams{
  double AK;
  double mK;
  double Apipi;
  double Epipi;
  std::vector<double> M;

  SimFitParams(): AK(1), mK(0.5), Apipi(1), Epipi(0.5), M(N,1){}
  SimFitParams(const double _AK, const double _mK, const double _Apipi, const double _Epipi, const double _M): AK(_AK), mK(_mK), Apipi(_Apipi), Epipi(_Epipi), M(N,_M){}
  SimFitParams(const SimFitParams<N> &r) = default;
  SimFitParams(SimFitParams<N> &&r) = default;

  SimFitParams<N>& operator=(const SimFitParams<N> &r) = default;
  SimFitParams<N>& operator=(SimFitParams<N> &&r) = default;
    
  inline int size() const{ return N+4; }

  inline std::string print() const{
    std::ostringstream os; os << "AK=" << AK << " mK=" << mK << " Apipi=" << Apipi << " Epipi=" << Epipi;
    for(int i=0;i<N;i++) os << " M[" << i << "]=" << M[i];
    return os.str();
  }
  inline void zero(){
    AK=mK=Apipi=Epipi=0.;
    for(int i=0;i<N;i++) M[i] = 0.;
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
  inline const double &operator()(const int i) const{ return const_cast<const double &>( const_cast<SimFitParams*>(this)->operator()(i) ); }

  typedef SimFitParams<N> ET_tag;
  
  template<typename U, typename std::enable_if<std::is_same<typename U::ET_tag, ET_tag>::value && !std::is_same<U,SimFitParams<N> >::value, int>::type = 0>
  SimFitParams(U&& expr): M(N){
    this->AK = expr[0];
    this->mK = expr[1];
    this->Apipi = expr[2];
    this->Epipi = expr[3];
    for(int i=4;i<N+4;i++) this->M[i-4] = expr[i];
  }
};

template<int N>
struct CPSfit::getElem<SimFitParams<N> >{
  static inline auto elem(const SimFitParams<N> &v, const int i)->decltype(v(i)){ return v(i); }
  static inline int common_properties(const SimFitParams<N> &v){ return 0; }
};

template<int N>
struct CPSfit::printStats<jackknifeDistribution<SimFitParams<N> > >{
  inline static std::string centralValue(const jackknifeDistribution<SimFitParams<N> > &d){
    auto best = d.best();
    std::ostringstream os; os << "(" << best.AK << ", " << best.mK << ", " << best.Apipi << ", " << best.Epipi;
    for(int i=0;i<N;i++) os << ", " << best.M[i];
    os << ")";
    return os.str();
  }
  inline static std::string error(const jackknifeDistribution<SimFitParams<N> > &d){
    auto err = d.standardError();
    std::ostringstream os; os << "(" << err.AK << ", " << err.mK << ", " << err.Apipi << ", " << err.Epipi;
    for(int i=0;i<N;i++) os << ", " << err.M[i];
    os << ")";
    return os.str();
  }
};

template<int N>
class FitKtoPiPiSim{
public:
  typedef SimFitParams<N> Params;
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

  static inline int Nparams(){ return N+4; }
};

  
GENERATE_PARSER_GM(FitKtoPiPiSim<10>::Params, FitKtoPiPiSim_10_Params_grammar, (double, AK)(double, mK)(double, Apipi)(double, Epipi)(std::vector<double>, M) )
GENERATE_PARSER_GM(FitKtoPiPiSim<7>::Params, FitKtoPiPiSim_7_Params_grammar, (double, AK)(double, mK)(double, Apipi)(double, Epipi)(std::vector<double>, M) )

template<int N>
struct ktopipiParamsPrinter<FitKtoPiPiSim<N> >: public distributionPrinter<jackknifeDistribution<typename FitKtoPiPiSim<N>::Params> >{
  void print(std::ostream &os, const jackknifeDistribution<typename FitKtoPiPiSim<N>::Params> &dist) const{
    typename FitKtoPiPiSim<N>::Params cen = dist.best();
    typename FitKtoPiPiSim<N>::Params err = dist.standardError();
    os << "AK = (" << cen.AK << " +- " << err.AK << ") mK = (" << cen.mK << " +- " << err.mK << ") "
       << "Apipi = (" << cen.Apipi << " +- " << err.Apipi << ") Epipi = (" << cen.Epipi << " +- " << err.Epipi << ")";
    for(int i=0;i<N;i++) os << " M[" << i << "] = (" << cen.M[i] << " +- " << err.M[i] << ")";
  }
};






template<typename FitFunc>
struct fitReturnType{};

template<>
struct fitReturnType<FitKtoPiPi>{ typedef std::vector<jackknifeDistribution<FitKtoPiPi::Params> > type; };

template<>
struct fitReturnType<FitKtoPiPiWithConstant>{ typedef std::vector<jackknifeDistribution<FitKtoPiPiWithConstant::Params> > type; };

template<>
struct fitReturnType<FitKtoPiPiTwoExp>{ typedef std::vector<jackknifeDistribution<FitKtoPiPiTwoExp::Params> > type; };

template<int N>
struct fitReturnType<FitKtoPiPiSim<N> >{ typedef jackknifeDistribution<typename FitKtoPiPiSim<N>::Params> type; };
  
template<>
struct fitReturnType<FitKtoPiPiSim<7> >{ typedef jackknifeDistribution<typename FitKtoPiPiSim<10>::Params> type; }; //return parameters converted to 10-basis

//Call static method 'call' on object with shared inputs type Inputs 
template<typename Inputs, template<typename> class F>
inline void fitfuncCall(const KtoPiPiFitFunc fitfunc, const Inputs &inputs){
  switch(fitfunc){
  case FitSeparate:
    return F<FitKtoPiPi>::call(inputs);
  case FitSimultaneous:
    return F<FitKtoPiPiSim<10> >::call(inputs);
  case FitSimultaneousChiralBasis:
    return F<FitKtoPiPiSim<7> >::call(inputs);
  case FitSeparateWithConstant:
    return F<FitKtoPiPiWithConstant>::call(inputs);
  case FitSeparateTwoExp:
    return F<FitKtoPiPiTwoExp>::call(inputs);
  default:
    error_exit(std::cout << "fitfuncCall(..) Unknown fit function " << fitfunc << std::endl);
  }
}


#endif

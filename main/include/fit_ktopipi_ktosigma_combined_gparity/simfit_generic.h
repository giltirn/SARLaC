#ifndef KTOPIPI_SIMFIT_GENERIC_H_
#define KTOPIPI_SIMFIT_GENERIC_H_

#include<config.h>
#include<utils/macros.h>

#include<containers/tagged_value_container.h>
#include<tensors/dual_number.h>

//Generic and flexible one and two state fits

CPSFIT_START_NAMESPACE

//Different fit functions require different numbers of parameters
struct SimFitCoordGen{
  double t; //time coordinate
  double tsep_k_snk;

  /*For the data for which this is the coordinate, provide the mapping of the inner parameter name to that in the complete set of parameters
    Here "inner parameter names" are the hard-coded namings for the parameters of the common fit form; eg for FitSimGenTwiState below these are Asnk0, Asnk1, E0, E1, etc
    The names in the full set are arbitrary

    As an example, imagine you are fitting K->pipi and K->sigma
    - there are 2 sink operators for which we may decide to name their overlaps with the ground state as  Apipi0 and Asigma0, and with the excited state as Apipi1 and Asigma1
    - there are a number of common parameters including E0, E1, AK, etc

    For a piece of K->pipi data, the param map should then be: Asnk0 -> Apipi0   Asnk1 -> Apipi1     AK -> AK, E0 -> E0....
    For a piece of K->sigma data, it should be               : Asnk0 -> Asigma0, Asnk1 -> Asigma1,   ....
    etc
  */
  std::unordered_map<std::string, std::string> const* param_map;  

  SimFitCoordGen() = default;
  SimFitCoordGen(const double t, const double tsep_k_snk, std::unordered_map<std::string, std::string> const* param_map): t(t), tsep_k_snk(tsep_k_snk),param_map(param_map){}
};

class FitSimGenTwoState{
  int nparams;
public:
  typedef taggedValueContainer<double,std::string> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoordGen GeneralizedCoordinate;

  inline const std::string &paramTag(const std::string &pname, const GeneralizedCoordinate &coord) const{ 
    auto it = coord.param_map->find(pname);
    if(it == coord.param_map->end()) error_exit(std::cout << "FitSimGenTwoState::paramTag could not find parameter name " << pname << std::endl);
    return it->second;
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    size_t idx; double v;

#define BOOSTIT(NM)  v = p(idx, paramTag(#NM,x)); auto NM = b(v, idx)
    BOOSTIT(AK);
    BOOSTIT(mK);
    BOOSTIT(Asnk0); //snk op with state 0
    BOOSTIT(Asnk1); //snk op with state 1
    BOOSTIT(E0);
    BOOSTIT(E1);
    BOOSTIT(M0);
    BOOSTIT(M1);
#undef BOOSTIT
    
    return AK*Asnk0*M0*exp(-E0 * x.tsep_k_snk)*exp( -(mK - E0)*x.t )/sqrt(2.) +
      AK*Asnk1*M1*exp(-E1 * x.tsep_k_snk)*exp( -(mK - E1)*x.t )/sqrt(2.);
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenTwoState(const int nparams): nparams(nparams){}

  inline int Nparams() const{ return nparams; }
};


class FitSimGenThreeState{
  int nparams;
public:
  typedef taggedValueContainer<double,std::string> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoordGen GeneralizedCoordinate;

  inline const std::string &paramTag(const std::string &pname, const GeneralizedCoordinate &coord) const{ 
    auto it = coord.param_map->find(pname);
    if(it == coord.param_map->end()) error_exit(std::cout << "FitSimGenThreeState::paramTag could not find parameter name " << pname << std::endl);
    return it->second;
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    size_t idx; double v;

#define BOOSTIT(NM)  v = p(idx, paramTag(#NM,x)); auto NM = b(v, idx)
    BOOSTIT(AK);
    BOOSTIT(mK);
    BOOSTIT(Asnk0); //snk op with state 0
    BOOSTIT(Asnk1); //snk op with state 1
    BOOSTIT(Asnk2); //snk op with state 2
    BOOSTIT(E0);
    BOOSTIT(E1);
    BOOSTIT(E2);
    BOOSTIT(M0);
    BOOSTIT(M1);
    BOOSTIT(M2);
#undef BOOSTIT
    
    return AK*Asnk0*M0*exp(-E0 * x.tsep_k_snk)*exp( -(mK - E0)*x.t )/sqrt(2.) +
      AK*Asnk1*M1*exp(-E1 * x.tsep_k_snk)*exp( -(mK - E1)*x.t )/sqrt(2.) +
      AK*Asnk2*M2*exp(-E2 * x.tsep_k_snk)*exp( -(mK - E2)*x.t )/sqrt(2.);
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenThreeState(const int nparams): nparams(nparams){}

  inline int Nparams() const{ return nparams; }
};




class FitSimGenMultiState{
  int nparams;
  int nstate;
public:
  typedef taggedValueContainer<double,std::string> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoordGen GeneralizedCoordinate;

  static inline const std::string &paramTag(const std::string &pname, const GeneralizedCoordinate &coord){ 
    auto it = coord.param_map->find(pname);
    if(it == coord.param_map->end()) error_exit(std::cout << "FitSimGenMultiState::paramTag could not find parameter name " << pname << std::endl);
    return it->second;
  }
  template<typename T, typename Boost>
  static inline T boostit(const std::string &pname, const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b){
    double v; size_t idx; 
    v = p(idx, paramTag(pname,x)); //get idx and value
    return b(v,idx);
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    size_t idx; double v;

    auto AK = boostit<T,Boost>("AK", x,p,b);
    auto mK = boostit<T,Boost>("mK", x,p,b);
    
    T out(0.);
    
    for(int s=0;s<nstate;s++){
      auto Asnk = boostit<T,Boost>(stringize("Asnk%d",s), x,p,b);
      auto E = boostit<T,Boost>(stringize("E%d",s), x,p,b);
      auto M = boostit<T,Boost>(stringize("M%d",s), x,p,b);
	
      out = out + AK*Asnk*M*exp(-E * x.tsep_k_snk)*exp( -(mK - E)*x.t )/sqrt(2.);
    }
    return out;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenMultiState(const int nparams, const int nstate): nparams(nparams), nstate(nstate){}

  inline int Nparams() const{ return nparams; }
};



//Fit the weighted average of the data with the kaon time dependence and operator normalization removed
//Here the SimFitCoordGen::t is interpreted as the operator->sink separation, and tsep_k_snk is not used
class FitSimGenMultiStateWavg{
  int nparams;
  int nstate;
public:
  typedef taggedValueContainer<double,std::string> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoordGen GeneralizedCoordinate;

  static inline const std::string &paramTag(const std::string &pname, const GeneralizedCoordinate &coord){ 
    auto it = coord.param_map->find(pname);
    if(it == coord.param_map->end()) error_exit(std::cout << "FitSimGenMultiState::paramTag could not find parameter name " << pname << std::endl);
    return it->second;
  }
  template<typename T, typename Boost>
  static inline T boostit(const std::string &pname, const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b){
    double v; size_t idx; 
    v = p(idx, paramTag(pname,x)); //get idx and value
    return b(v,idx);
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    size_t idx; double v;
   
    T out(0.);
    
    for(int s=0;s<nstate;s++){
      auto Asnk = boostit<T,Boost>(stringize("Asnk%d",s), x,p,b);
      auto E = boostit<T,Boost>(stringize("E%d",s), x,p,b);
      auto M = boostit<T,Boost>(stringize("M%d",s), x,p,b);
	
      out = out + Asnk*M*exp( -E*x.t )/sqrt(2.);
    }
    return out;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenMultiStateWavg(const int nparams, const int nstate): nparams(nparams), nstate(nstate){}

  inline int Nparams() const{ return nparams; }
};



CPSFIT_END_NAMESPACE

#endif

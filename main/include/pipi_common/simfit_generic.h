#ifndef SIMFIT_GENERIC_H_
#define SIMFIT_GENERIC_H_

#include<config.h>
#include<utils/macros.h>

#include<containers/tagged_value_container.h>
#include<tensors/dual_number.h>

//Generic and flexible one and two state fits

CPSFIT_START_NAMESPACE

//Different fit functions require different numbers of parameters
struct SimFitCoordGen{
  double t; //time coordinate

  /*The data typically is symmetric as   C(Lt - fold_offset - t) ~ C(t) 
    fold_offset should be  2*tsep_pipi  for pipi2pt,   tsep_pipi  for pipi->sigma and 0 for sigma 2pt
  */
  int fold_offset; 

  /*For the data for which this is the coordinate, provide the mapping of the inner parameter name to that in the complete set of parameters
    Here "inner parameter names" are the hard-coded namings for the parameters of the common fit form; eg for FitSimGenOneState below these are Asrc, Asnk, E, Csys
    The names in the full set are arbitrary

    As an example, imagine you are fitting pipi 2pt, sigma 2pt and pipi->sigma data to a single shared state
    - there are 2 operators for which we may decide to name their overlaps with the ground state as  Apipi and Asigma
    - there are 3 constants associated with the 3 correlators, which we may name   Cpipi2pt, Cpipitosigma, Csigma2pt
    - the energy of the state we might call Egnd

    For a piece of pipi 2pt data, the param map should then be:  Asrc -> Apipi,  Asnk -> Apipi, E -> Egnd, Csys -> Cpipi2pt
    For a piece of pipi->sigma data, it should be             :  Asrc -> Apipi,  Asnk -> Asigma, E-> Egnd, Csys -> Cpipitosigma
    etc
  */
  std::unordered_map<std::string, std::string> const* param_map;  

  SimFitCoordGen() = default;
  SimFitCoordGen(const double t, std::unordered_map<std::string, std::string> const* param_map, int fold_offset): t(t), param_map(param_map), fold_offset(fold_offset){}
};

class FitSimGenOneState{
  int Lt;
  double Ascale;
  double Cscale;
  int nparams;
public:
  typedef taggedValueContainer<double,std::string> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoordGen GeneralizedCoordinate;

  inline const std::string &paramTag(const std::string &pname, const GeneralizedCoordinate &coord) const{ 
    auto it = coord.param_map->find(pname);
    if(it == coord.param_map->end()) error_exit(std::cout << "FitSimGenOneState::paramTag could not find parameter name " << pname << std::endl);
    return it->second;
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    const std::string &Asrctag = paramTag("Asrc",x); //<0|O_src(tsrc)|state>
    const std::string &Asnktag = paramTag("Asnk",x); //<state|O_snk(tsnk)|0>
    const std::string &Etag = paramTag("E",x); //Energy of state
    const std::string &Csystag = paramTag("Csys",x);   //constant term associated with the system of states produced by the operators
   
    size_t idx; double v;
#define BOOSTIT(NM) v = p(idx, NM ## tag); auto NM = b(v, idx)
    BOOSTIT(Asrc);
    BOOSTIT(Asnk);
    BOOSTIT(E);
    BOOSTIT(Csys);
#undef BOOSTIT

    return Asrc * Asnk * Ascale * ( exp(-E*x.t) + exp(-E*(Lt-x.fold_offset-x.t)) ) + Csys * Cscale;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenOneState(const int Lt, const int nparams, const double Ascale=1e13, const double Cscale=1e13): Lt(Lt), nparams(nparams), Ascale(Ascale), Cscale(Cscale){}

  inline int Nparams() const{ return nparams; }
};


class FitSimGenTwoState{
  int Lt;
  double Ascale;
  double Cscale;
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
    const std::string &Asrc0tag = paramTag("Asrc0",x); //src op with state 0
    const std::string &Asrc1tag = paramTag("Asrc1",x); //src op with state 1
    const std::string &Asnk0tag = paramTag("Asnk0",x); //snk op with state 0
    const std::string &Asnk1tag = paramTag("Asnk1",x); //snk op with state 1
    const std::string &E0tag = paramTag("E0",x);
    const std::string &E1tag = paramTag("E1",x);
    const std::string &Csystag = paramTag("Csys",x);
   
    size_t idx; double v;
#define BOOSTIT(NM) v = p(idx, NM ## tag); auto NM = b(v, idx)
    BOOSTIT(Asrc0);
    BOOSTIT(Asrc1);
    BOOSTIT(Asnk0);
    BOOSTIT(Asnk1);
    BOOSTIT(E0);
    BOOSTIT(E1);
    BOOSTIT(Csys);
#undef BOOSTIT

    return Asrc0 * Asnk0 * Ascale * ( exp(-E0*x.t) + exp(-E0*(Lt-x.fold_offset-x.t)) ) 
      + Asrc1 * Asnk1 * Ascale * ( exp(-E1*x.t) + exp(-E1*(Lt-x.fold_offset-x.t)) ) 
      + Csys * Cscale;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenTwoState(const int Lt, const int nparams, const double Ascale=1e13, const double Cscale=1e13): Lt(Lt), nparams(nparams), Ascale(Ascale), Cscale(Cscale){}

  inline int Nparams() const{ return nparams; }
};



class FitSimGenThreeState{
  int Lt;
  double Ascale;
  double Cscale;
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
#define BOOSTIT(NM) v = p(idx, paramTag(#NM,x)); auto NM = b(v, idx)
    BOOSTIT(Asrc0); //src op with state 0
    BOOSTIT(Asrc1);
    BOOSTIT(Asrc2);
    BOOSTIT(Asnk0);
    BOOSTIT(Asnk1);
    BOOSTIT(Asnk2);
    BOOSTIT(E0);
    BOOSTIT(E1);
    BOOSTIT(E2);
    BOOSTIT(Csys);
#undef BOOSTIT

    return Asrc0 * Asnk0 * Ascale * ( exp(-E0*x.t) + exp(-E0*(Lt-x.fold_offset-x.t)) ) 
      + Asrc1 * Asnk1 * Ascale * ( exp(-E1*x.t) + exp(-E1*(Lt-x.fold_offset-x.t)) ) 
      + Asrc2 * Asnk2 * Ascale * ( exp(-E2*x.t) + exp(-E2*(Lt-x.fold_offset-x.t)) ) 
      + Csys * Cscale;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenThreeState(const int Lt, const int nparams, const double Ascale=1e13, const double Cscale=1e13): Lt(Lt), nparams(nparams), Ascale(Ascale), Cscale(Cscale){}

  inline int Nparams() const{ return nparams; }
};




//We can stabilize the fit by enforcing that E2 > E1 > E0 > 0 by parameterizing in terms of log(E0)  ,  log(E1 - E0), log(E2 - E1)
class FitSimGenThreeStateLogEdiff{
  int Lt;
  double Ascale;
  double Cscale;
  int nparams;
public:
  typedef taggedValueContainer<double,std::string> Params;
  typedef double ValueType;
  typedef Params ParameterType;
  typedef Params ValueDerivativeType;
  typedef SimFitCoordGen GeneralizedCoordinate;

  inline const std::string &paramTag(const std::string &pname, const GeneralizedCoordinate &coord) const{ 
    auto it = coord.param_map->find(pname);
    if(it == coord.param_map->end()) error_exit(std::cout << "FitSimGenThreeStateLogEdiff::paramTag could not find parameter name " << pname << std::endl);
    return it->second;
  }

  template<typename T, typename Boost>
  T eval(const GeneralizedCoordinate &x, const ParameterType &p, const Boost &b) const{
    size_t idx; double v;
#define BOOSTIT(NM) v = p(idx, paramTag(#NM,x)); auto NM = b(v, idx)
    BOOSTIT(Asrc0); //src op with state 0
    BOOSTIT(Asrc1);
    BOOSTIT(Asrc2);
    BOOSTIT(Asnk0);
    BOOSTIT(Asnk1);
    BOOSTIT(Asnk2);
    BOOSTIT(logE0);
    BOOSTIT(logE1mE0);
    BOOSTIT(logE2mE1);
    BOOSTIT(Csys);
#undef BOOSTIT

    T E0 = exp(logE0);
    T E1 = exp(logE1mE0) + E0;
    T E2 = exp(logE2mE1) + E1;    

    return Asrc0 * Asnk0 * Ascale * ( exp(-E0*x.t) + exp(-E0*(Lt-x.fold_offset-x.t)) ) 
      + Asrc1 * Asnk1 * Ascale * ( exp(-E1*x.t) + exp(-E1*(Lt-x.fold_offset-x.t)) ) 
      + Asrc2 * Asnk2 * Ascale * ( exp(-E2*x.t) + exp(-E2*(Lt-x.fold_offset-x.t)) ) 
      + Csys * Cscale;
  }

  inline ValueType value(const GeneralizedCoordinate &x, const ParameterType &p) const{
    return eval<double>(x,p,[&](const double a, const int i){ return a; });
  }
  inline ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &x, const ParameterType &p) const{
    ValueDerivativeType d = p;
    for(int i=0;i<nparams;i++) d(i) = eval<dual>(x,p,[&](const double v, const int j){ return dual(v, j==i ? 1.:0.); }).xp;
    return d;
  }

  FitSimGenThreeStateLogEdiff(const int Lt, const int nparams, const double Ascale=1e13, const double Cscale=1e13): Lt(Lt), nparams(nparams), Ascale(Ascale), Cscale(Cscale){}

  inline int Nparams() const{ return nparams; }
};





CPSFIT_END_NAMESPACE

#endif

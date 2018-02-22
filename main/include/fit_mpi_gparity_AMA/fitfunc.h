#ifndef _FIT_MPI_GPARITY_AMA_FITFUNC_H
#define _FIT_MPI_GPARITY_AMA_FITFUNC_H

#include<fit.h>
#include<containers.h>

class TypeInfo{
  std::map<std::string, int> pmap;
  std::map<int, std::string> upmap;
  std::map<int, FitFuncType> fitfuncs;
  int cur;
public:
  TypeInfo(): cur(0){}

  void registerType(const std::string &type, const FitFuncType ff){
    std::map<std::string, int>::const_iterator it = pmap.find(type);
    if(it != pmap.end())
      assert(fitfuncs[it->second] == ff);
    else{
      pmap[type] = cur;
      upmap[cur] = type;
      fitfuncs[cur] = ff;
      ++cur;
    }
  }

  FitFuncType fitFunc(const int idx) const{
    std::map<int, FitFuncType>::const_iterator it = fitfuncs.find(idx);
    if(it == fitfuncs.end()) error_exit(std::cout << "TypeInfo::fitFunc could not find index " << idx << " in type map\n");
    return it->second;
  }
  
  int typeIdx(const std::string &type) const{
    std::map<std::string, int>::const_iterator it = pmap.find(type);
    if(it == pmap.end()) error_exit(std::cout << "TypeInfo::typeIdx could not find string " << type << " in type map\n");
    return it->second;
  }

  inline bool containsType(const std::string &type) const{ return pmap.count(type) != 0; }
  
  const std::string & typeName(const int idx) const{
    std::map<int, std::string>::const_iterator it = upmap.find(idx);
    if(it == upmap.end()) error_exit(std::cout << "TypeInfo::typeName could not find index " << idx << " in type map\n");
    return it->second;
  }

  inline int nTypes() const{ return pmap.size(); }
};
    
  

//Coordinate, parameters and param derivatives for aggregate fit func
struct Coord{
  double t; //time
  int type;
  TypeInfo const* pmap;
  
  Coord(){}
  Coord(const double _t, const std::string _type, const TypeInfo &_pmap): t(_t), pmap(&_pmap){
    type = _pmap.typeIdx(_type);
  }
};
std::ostream & operator<<(std::ostream &os, const Coord &c){
  os << "(" << c.t << "," << c.pmap->typeName(c.type) << ")";
  return os;
}

class filterCoordTrange{
  double t_min;
  double t_max;
  bool _invert;
public:
  filterCoordTrange(const double _t_min, const double _t_max): t_min(_t_min), t_max(_t_max), _invert(false){}

  void invert(){ _invert = !_invert; } 
  
  template<typename T>
  bool accept(const Coord &x, const T &y) const{
    bool cond = x.t >= t_min && x.t <= t_max;
    return _invert ? !cond : cond;
  }
};

class FitMpi{
public:
  typedef double ValueType;
  typedef parameterVector<double> ParameterType;
  typedef parameterVector<double> ValueDerivativeType;
  typedef Coord GeneralizedCoordinate;
private:
  FitCosh fcosh;
  FitSinh fsinh;
  FitExp fexp;
  
  inline StandardFitParams reduce(const ParameterType &p, const int type) const{
    const double &m = p(0);
    const double &A = p(type+1);
    return StandardFitParams(A,m);
  }
  int nparams;
  int ntypes;
  std::vector<StandardFitFuncBase const*> fitfuncs;
public:

  FitMpi(const double _Lt, const TypeInfo &pmap): fcosh(_Lt), fsinh(_Lt), ntypes(pmap.nTypes()){
    nparams = ntypes+1; //mass, A0, A1, ....
    fitfuncs.resize(ntypes);
    for(int i=0;i<ntypes;i++){
      switch(pmap.fitFunc(i)){
      case FCosh:
	fitfuncs[i] = &fcosh; break;
      case FSinh:
	fitfuncs[i] = &fsinh; break;
      case FExp:
	fitfuncs[i] = &fexp; break;
      default:
	assert(0);
      }
    }    
  }
  
  inline ValueType value(const GeneralizedCoordinate &c, const ParameterType &p) const{
    return fitfuncs[c.type]->value(c.t,reduce(p,c.type));
  }

  ValueDerivativeType parameterDerivatives(const GeneralizedCoordinate &c, const ParameterType &p) const{
    ValueDerivativeType yderivs(nparams);
    yderivs.zero();
    StandardFitParamDerivs subderivs = fitfuncs[c.type]->parameterDerivatives(c.t,reduce(p,c.type));
    yderivs(0) = subderivs.dm;
    yderivs(c.type+1) = subderivs.dA;
    return yderivs;
  }

  inline int Nparams() const{ return nparams; }
};


class FitMpiEffectiveMass{
  FitCosh fcosh;
  FitSinh fsinh;
  FitExp fexp;
  StandardFitFuncBase* fitfunc;
public:
  typedef double ValueType;
  typedef singleValueContainer<double> ParameterType;
  typedef singleValueContainer<double> ValueDerivativeType;
  typedef double GeneralizedCoordinate;

  FitMpiEffectiveMass(const double _Lt, const std::string &_type, const TypeInfo &pmap): fcosh(_Lt), fsinh(_Lt){
    int idx = pmap.typeIdx(_type);
    switch(pmap.fitFunc(idx)){
    case FCosh:
      fitfunc = new FitCosh(_Lt); break;
    case FSinh:
      fitfunc = new FitSinh(_Lt); break;
    case FExp:
      fitfunc = new FitExp; break;
    default:
      assert(0);
    }
  }    

  inline double value(const double t, const ParameterType &params) const{    
    StandardFitParams p(1000, *params);
    return fitfunc->value(t,p)/fitfunc->value(t+1,p);
  }
  
  ValueDerivativeType parameterDerivatives(const double t, const ParameterType &params) const{
    ValueDerivativeType yderivs;
    StandardFitParams p(1000, *params);
    
    double value_t = fitfunc->value(t,p);
    StandardFitParamDerivs subderivs_t = fitfunc->parameterDerivatives(t,p);

    double value_tp1 = fitfunc->value(t+1,p);
    StandardFitParamDerivs subderivs_tp1 = fitfunc->parameterDerivatives(t+1,p);
 
    *yderivs = subderivs_t.dm/value_tp1 - value_t/value_tp1/value_tp1 * subderivs_tp1.dm;
    return yderivs;
  }

  inline int Nparams() const{ return 1; }

  ~FitMpiEffectiveMass(){
    delete fitfunc;
  }
  
};





#endif
